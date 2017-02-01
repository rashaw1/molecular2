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
#include <cmath>
#include "mvector.hpp"
#include <libint2.hpp>

// Constructor
ALMOSCF::ALMOSCF(Molecule& m, Fock& f) : molecule(m), focker(f) 
{	
	// Zero out energies
	dimer_energy = e_frz = e_pol = e_ct = e_int = 0.0;
	
	// Add the fragment Fock matrices and perform monomer calculations
	std::vector<Fragment>& frags = molecule.getLog().getFragments(); 
	nfrags = frags.size();
	int start = 0;
	int nf = 0;
	for (int i = 0; i < frags.size(); i++) {
		
		frags[i].buildShellBasis();
		frags[i].calcEnuc();
		std::vector<libint2::Shell>& shells = frags[i].getBasis().getIntShells();

		int f_nbfs = ints[nf].nbasis(shells);
		ints.push_back(IntegralEngine(frags[i], focker.getIntegrals(), start, start+f_nbfs));
		start += f_nbfs;
		fragments.push_back(FockFragment(ints[nf], frags[i]));
		
		SCF hf(frags[i], fragments[fragments.size()-1]);
		hf.rhf();
		monomer_energies.push_back(hf.getEnergy()); 

		nf++;
	}

	molecule.getLog().title("ALMO Calculation");
	molecule.getLog().print("Making initial density...\n");
	makeDens(); 
	molecule.getLog().localTime();
	
}

// Routines
void ALMOSCF::makeDens() {

	// Make inverse overlap metric
	Matrix sigma = focker.getS();
	Matrix& hcore = focker.getHCore();
	
	//Build T matrix
	int nbfs = sigma.nrows();
	int nocc = focker.getMolecule().getNel() / 2;
	
	Eigen::MatrixXd S = Eigen::MatrixXd::Zero(nbfs, nbfs);
	Eigen::MatrixXd H = Eigen::MatrixXd::Zero(nbfs, nbfs);
	for (int i = 0; i < nbfs; i++) {
		for (int j = 0; j <= i; j++) {
			S(i, j) = S(j, i) = sigma(i, j); 
			H(i, j) = H(j, i) = hcore(i, j);
		}
	}
	
	Eigen::MatrixXd Tocc = Eigen::MatrixXd::Zero(nbfs, nocc); 
	int row_offset = 0; int col_offset = 0; 
	for (auto f : fragments) {
		Matrix& f_cp = f.getCP(); 
		int f_nocc = f.getMolecule().getNel() / 2;
		int f_nbfs = f.getHCore().nrows(); 
		
		for (int col = 0; col < f_nocc; col++)
			for (int row = 0; row < f_nbfs; row++)
				Tocc(row+row_offset, col+col_offset) = f_cp(row, col); 
		
		row_offset += f_nbfs; 
		col_offset += f_nocc; 
	}
	
	P = Tocc.transpose() * S * Tocc;
	P = Tocc * P.inverse() * Tocc.transpose(); 
	
	// Build Fock matrix
	Eigen::MatrixXd F(nbfs, nbfs);
	IntegralEngine& integrals = focker.getIntegrals();
	for (int i = 0; i < nbfs; i++) {
		for (int j = 0; j < nbfs; j++) {
			F(i, j) = H(i, j);
			
			for (int r = 0; r < nbfs; r++)
				for (int s = 0; s < nbfs; s++)
					F(i, j) += (2*integrals.getERI(i, j, r, s) - integrals.getERI(i, s, r, j)) * P(s, r); 
		}
	}
	
	double energy = (P * (H + F)).trace() + focker.getMolecule().getEnuc();
	//for (auto en : monomer_energies) energy -= en; 
	std::cout << energy << std::endl;  
	
	
	
}

// Calculate the scf energy
void ALMOSCF::calcE()
{
	
}

void ALMOSCF::rscf()
{

}

