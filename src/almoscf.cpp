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
#include <Eigen/Cholesky>
#include <Eigen/Eigenvalues>

// Constructor
ALMOSCF::ALMOSCF(Molecule& m, Fock& f) : molecule(m), focker(f) 
{	
	// Zero out energies
	dimer_energy = e_frz = e_pol = e_ct = e_int = 0.0;
	MAX = 6;
	diis.init(MAX, m.getLog().diis());
	P = Eigen::MatrixXd::Zero(focker.getHCore().nrows(), focker.getHCore().nrows());
}

void ALMOSCF::setFragments(bool unrestricted)
{
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
		IntegralEngine new_engine(frags[i], focker.getIntegrals(), start, start+f_nbfs);
		ints.push_back(new_engine);
		fragments.push_back(FockFragment(ints[nf], frags[i], start, start+f_nbfs));
		start += f_nbfs;
		
		SCF hf(frags[i], fragments[fragments.size()-1]);
		if (unrestricted) hf.uhf();
		else hf.rhf(); 
		monomer_energies.push_back(hf.getEnergy()); 

		nf++;
	}
}

// Routines
void ALMOSCF::compute() {

	// Make inverse overlap metric
	Matrix sigma = focker.getS();
	Matrix& hcore = focker.getHCore();
	
	//Build T matrix
	int nbfs = sigma.nrows();
	int nocc = focker.getMolecule().getNel() / 2;
	
	Eigen::MatrixXd P_old = P; 
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
	
	delta_d = (P - P_old).norm(); 
	
	// Build Fock matrix
	F = Eigen::MatrixXd::Zero(nbfs, nbfs);
	IntegralEngine& integrals = focker.getIntegrals();
	for (int i = 0; i < nbfs; i++) {
		for (int j = 0; j < nbfs; j++) {
			F(i, j) = H(i, j);
			
			for (int r = 0; r < nbfs; r++)
				for (int s = 0; s < nbfs; s++)
					F(i, j) += (2*integrals.getERI(i, j, r, s) - integrals.getERI(i, s, r, j)) * P(s, r); 
		}
	}
	
	if (molecule.getLog().diis()) { // Archive for averaging
		if (focks.size() == MAX) {
			focks.erase(focks.begin());
		}
		focks.push_back(F);
	}		
	
	// Calculate errors
	Eigen::MatrixXd Q = - S * P; 
	for (int i = 0; i < nbfs; i++) Q(i, i) += 1.0; 
	Eigen::MatrixXd QFP = Q * F; 
	if (molecule.getLog().diis()) {
		QFP = QFP * P; 
		
		std::vector<Vector> errs; 
		int offset = 0;
		for (int i = 0; i < fragments.size(); i++) {
			int f_nbfs = fragments[i].getHCore().nrows();
			Eigen::MatrixXd e = QFP.block(offset, offset, f_nbfs, f_nbfs) * fragments[i].Sxx;
			e = (e - e.transpose()).eval(); 
			
			Matrix new_err(f_nbfs, f_nbfs, 0.0);
			for (int i = 0; i < f_nbfs; i++)  
				for (int j = 0; j <= i; j++)
					new_err(i, j) = new_err(j, i) = e(i, j);
			errs.push_back(Vector(new_err, true)); 
			
			offset += f_nbfs; 
		}
		
		Vector weights = diis.compute(errs);
		if (focks.size() > 2) {
			// Average the fock matrices according to the weights
			offset = focks.size() - weights.size();
			F = Eigen::MatrixXd::Zero(nbfs, nbfs);
			for (int i = offset; i < focks.size(); i++)
				 F = F + weights[i-offset]*focks[i]; 
		}
		 
		QFP = Q * F; 
	} 
	
	double new_dimer_energy = (P * (H + F)).trace() + focker.getMolecule().getEnuc();
	delta_e = new_dimer_energy - dimer_energy;
	dimer_energy = new_dimer_energy;
	
	Eigen::MatrixXd QFQ = QFP * Q.transpose();
	QFP = QFP * P; 
	Eigen::MatrixXd PFP = P * F * P; 
	
	for (int i = 0; i < fragments.size(); i++)
		fragments[i].buildFock(QFQ, QFP, PFP);
	
}

// Calculate the perturbative correction
void ALMOSCF::perturb(bool order4)
{
	int nbfs = focker.getHCore().nrows(); 
	int nocc = focker.getMolecule().getNel() / 2;
	int nvirt = nbfs - nocc;
	
	Eigen::MatrixXd T = Eigen::MatrixXd::Zero(nbfs, nocc); 
	Eigen::MatrixXd V = Eigen::MatrixXd::Zero(nbfs, nvirt);
	int row_offset = 0; int occ_col_offset = 0; int virt_col_offset = 0;
	for (auto f : fragments) {
		Matrix& f_cp = f.getCP(); 
		int f_nocc = f.getMolecule().getNel() / 2;
		int f_nbfs = f.getHCore().nrows(); 
		int f_nvirt = f_nbfs - f_nocc; 
		
		for (int col = 0; col < f_nocc; col++)
			for (int row = 0; row < f_nbfs; row++)
				T(row+row_offset, col+occ_col_offset) = f_cp(row, col); 
		
		for (int col = f_nocc; col < f_nbfs; col++)
			for (int row = 0; row < f_nbfs; row++)
				V(row+row_offset, col+virt_col_offset-f_nocc) = f_cp(row, col);
		
		row_offset += f_nbfs; 
		occ_col_offset += f_nocc; 
		virt_col_offset += f_nvirt; 
	}
	
	// Form MO overlap matrix
	Eigen::MatrixXd sigma(nbfs, nbfs);
	Eigen::MatrixXd CMI = Eigen::MatrixXd::Zero(nbfs, nbfs);
	Matrix S = focker.getS();
	for (int i = 0; i < nbfs; i++) 
		for (int j = 0; j <= i; j++) 
			CMI(i, j) = CMI(j, i) = S(i, j);
	sigma.block(0, 0, nocc, nocc) = T.transpose() * CMI * T; 
	sigma.block(nocc, 0, nvirt, nocc) = V.transpose() * CMI * T;
	sigma.block(nocc, nocc, nvirt, nvirt) = V.transpose() * CMI * V; 
	Eigen::LLT<Eigen::MatrixXd> llt(sigma); 
	Eigen::MatrixXd C0 = llt.matrixL();
	C0 = C0.inverse();
	CMI.block(0, 0, nbfs, nocc) = T;
	CMI.block(0, nocc, nbfs, nvirt) = V;
	C0 = CMI * C0.transpose(); 
	
	sigma = C0.transpose() * F * C0; 
	Eigen::MatrixXd F0occ = sigma.block(0, 0, nocc, nocc);
	Eigen::MatrixXd F0virt = sigma.block(nocc, nocc, nvirt, nvirt); 
	Eigen::MatrixXd Fov = sigma.block(0, nocc, nocc, nvirt); 
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es_occ(F0occ);
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es_virt(F0virt);
	T = es_occ.eigenvectors();
	V = es_virt.eigenvectors(); 
	Fov = T.transpose() * Fov * V; 
	Eigen::VectorXd eps_i = es_occ.eigenvalues(); 
	Eigen::VectorXd eps_a = es_virt.eigenvalues();
	
	if (order4) {
		e_pert_2 = 0.0; 
		e_pert_4 = 0.0;
		
		Eigen::MatrixXd Fvo = V.transpose() * sigma.block(nocc, 0, nvirt, nocc) * T;
		for (int i = 0; i < nocc; i++) {
			double sum1 = 0.0, sum2 = 0.0;
			for (int a = 0; a < nvirt; a++) {
				auto delta = eps_i(i) - eps_a(a);
				auto val = Fov(i, a) * Fov(i, a) / delta; 
				e_pert_2 += val;
				sum1 += val; 
				sum2 += val / delta;
				
				for (int b = 0; b < nvirt; b++) {
					for (int j = 0; j < nocc; j++) {
						if (i == j) continue;
						e_pert_4 += Fov(i, a) * Fvo(a, j) * Fov(j, b) * Fvo(b, i) / ( delta * (eps_i(i) - eps_i(j)) * (eps_i(i) - eps_a(b)) ); 
					}
				}
			}
			e_pert_4 -= sum1 * sum2; 
		}
	} else {
		e_pert_2 = 0.0;
		for (int i = 0; i < nocc; i++)
			for (int a = 0; a < nvirt; a++)
				e_pert_2 += Fov(i, a) * Fov(i, a) / (eps_i(i) - eps_a(a)); 
	}
	
}

void ALMOSCF::rscf()
{
	molecule.getLog().title("ALMO Calculation");
	
	setFragments();
	
	molecule.getLog().initIteration(); 
	delta_d = 1.0; delta_e = 1.0; 
	bool converged = false;
	int iter = 0;
	while (!converged && iter < molecule.getLog().maxiter()) {
		compute(); 
		molecule.getLog().iteration(iter++, dimer_energy, delta_e, delta_d);
		converged = (fabs(delta_e) < molecule.getLog().converge()) && (fabs(delta_d) < molecule.getLog().converge());
	}
	
	if (converged) {
		double energy = dimer_energy;
		for (auto en : monomer_energies) energy -= en; 
		molecule.getLog().result("ALMO Interaction Energy", energy * Logger::TOKCAL, "kcal / mol"); 
		perturb(true);
		e_pert_2 *= 2.0;
		e_pert_4 *= 2.0;
		molecule.getLog().result("E(2)", e_pert_2 * Logger::TOKCAL, "kcal /mol"); 
		molecule.getLog().result("E(4)", e_pert_4 * Logger::TOKCAL, "kcal / mol"); 
		molecule.getLog().result("Total ALMO+CT interaction energy", (energy + e_pert_2 + e_pert_4) * Logger::TOKCAL, "kcal / mol");
	} else {
		molecule.getLog().result("ALMO SCF failed to converge");
	}
}

