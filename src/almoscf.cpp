/*
 *
 *   PURPOSE: To implement class SCF, which carries out HF self-consistent field calculations.
 *
 *   DATE        AUTHOR         CHANGES
 *   ===============================================================
 *   15/09/15    Robert Shaw    Original code.
 *
 */

#include "almoscf.hpp"
#include "scf.hpp"
#include "logger.hpp"
#include "integrals.hpp"
#include <cmath>
#include <Eigen/Dense>
#include "mvector.hpp"
#include <libint2.hpp>

// Constructor
ALMOSCF::ALMOSCF(Molecule& m, Fock& f) : molecule(m), focker(f) 
{
	// Build Sinv
	Matrix S = focker.getIntegrals().getOverlap();
	Eigen::MatrixXd s_(S.nrows(), S.ncols());
	for (int i = 0; i < S.nrows(); i++)
		for (int j = 0; j<= i; j++)
			s_(i, j) = s_(j, i) = S(i, j);
	s_ = s_.inverse();
	Sinv.assign(S.nrows(), S.ncols(), 0.0);
	for (int i = 0; i < S.nrows(); i++) 
		for (int j = 0; j <= i; j++)
			Sinv(i, j) = Sinv(j, i) = s_(i, j);
	
	// Zero out energies
	dimer_energy = e_frz = e_pol = e_ct = e_int = 0.0;
	
	// Add the fragment Fock matrices and perform monomer calculations
	std::vector<Fragment>& frags = molecule.getLog().getFragments(); 
	nfrags = frags.size();
	int start = 0;
	int shell_offset = 0; 
	Matrix T(S.nrows(), S.nrows(), 0.0);
	Matrix tx;
	Matrix ty;
	for (auto f : frags) {
		
		f.buildShellBasis();
		f.calcEnuc();
		std::vector<libint2::Shell>& shells = f.getBasis().getIntShells();
		
		int f_nbfs = focker.getIntegrals().nbasis(shells);
		fragments.push_back(FockFragment(start, start+f_nbfs, shell_offset, S, focker.getIntegrals(), f));
		shell_offset += shells.size(); 
		
		SCF hf(f, fragments[fragments.size()-1]);
		hf.rhf();
		monomer_energies.push_back(hf.getEnergy()); 

		Matrix& cp = fragments[fragments.size()-1].getCP();
		if (start == 0) tx = fragments[fragments.size()-1].getCP();
		else ty = fragments[fragments.size()-1].getCP();
		for (int i = 0; i < cp.nrows(); i++)
			for (int j = 0; j < cp.ncols(); j++)
				T(i+start, j+start) = cp(i, j);
		start += f_nbfs;
	}

	T = T.transpose() * S * T;
	Eigen::MatrixXd tst(T.nrows(), T.ncols());
	for (int i = 0; i < T.nrows(); i++)
		for (int j = 0; j < T.nrows(); j++)
			tst(i, j) = T(i, j);
	tst = tst.inverse();

	T.assign(tx.ncols(), ty.ncols(), 0.0);
	for (int i = 0; i < tx.ncols(); i++)
		for (int j = tx.ncols(); j < tx.ncols() + ty.ncols(); j++)
			T(i, j-tx.ncols()) = tst(i, j);
	(tx * T * ty.transpose()).print(); std::cout << std::endl;
	for (int i = 0; i < tx.ncols(); i++)
		for (int j = tx.ncols(); j < tx.ncols() + ty.ncols(); j++)
			T(i, j - tx.ncols()) = Sinv(i, j);
	T.print(); std::cout << std::endl;
	
}

// Routines

// Calculate the scf energy
void ALMOSCF::calcE()
{
	
}

void ALMOSCF::rscf()
{

}

