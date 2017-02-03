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
#include <libint2.hpp>

// Constructor
ALMOSCF::ALMOSCF(Molecule& m, Fock& f) : molecule(m), focker(f) 
{	
	// Zero out energies
	dimer_energy = e_frz = e_pol = e_ct = e_int = 0.0;
	MAX = 6;
	diis.init(MAX, m.getLog().diis());
	P = Matrix::Zero(focker.getHCore().rows(), focker.getHCore().rows());
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
	Matrix& S = focker.getS();
	Matrix P_old = P; 
	
	//Build T matrix
	int nbfs = S.rows();
	int nocc = focker.getMolecule().getNel() / 2;
	
	Matrix Tocc = Matrix::Zero(nbfs, nocc); 
	int row_offset = 0; int col_offset = 0; 
	for (auto f : fragments) {
		int f_nocc = f.getMolecule().getNel() / 2;
		int f_nbfs = f.getHCore().rows(); 
		
		Tocc.block(row_offset, col_offset, f_nbfs, f_nocc) = f.getCP().block(0, 0, f_nbfs, f_nocc);
		
		row_offset += f_nbfs; 
		col_offset += f_nocc; 
	}

	P = Tocc.transpose() * S * Tocc;
	P = Tocc * P.inverse() * Tocc.transpose(); 
	
	delta_d = (P - P_old).norm(); 
	
	Matrix P_2 = 2.0 * P; 
	// Build Fock matrix
	if (molecule.getLog().direct())
		focker.formJKdirect(focker.getIntegrals().getPrescreen(), P_2);
	else 
		focker.formJK(P_2); 
	
	focker.makeFock(); 
	Matrix& F = focker.getFockAO(); 
	
	// Calculate errors
	Matrix Q = - S * P; 
	for (int i = 0; i < nbfs; i++) Q(i, i) += 1.0; 
	Matrix QFP = Q * F; 
	if (molecule.getLog().diis()) {
		QFP = QFP * P; 
		
		std::vector<Vector> errs; 
		int offset = 0;
		for (int i = 0; i < fragments.size(); i++) {
			int f_nbfs = fragments[i].getHCore().rows();
			Matrix e = QFP.block(offset, offset, f_nbfs, f_nbfs) * fragments[i].Sxx;
			e = (e - e.transpose()).eval(); 

			Vector new_err(Eigen::Map<Vector>(e.data(), e.cols()*e.rows()));
			errs.push_back(new_err); 
			
			offset += f_nbfs; 
		}
		
		Vector weights = diis.compute(errs);
		focker.average(weights);
		F = focker.getFockAO();
		 
		QFP = Q * F; 
	} 
	
	Matrix& H = focker.getHCore();
	double new_dimer_energy = (P * (H + F)).trace() + focker.getMolecule().getEnuc();
	delta_e = new_dimer_energy - dimer_energy;
	dimer_energy = new_dimer_energy;
	
	Matrix QFQ = QFP * Q.transpose();
	QFP = QFP * P; 
	Matrix PFP = P * F * P; 
	
	for (int i = 0; i < fragments.size(); i++)
		fragments[i].buildFock(QFQ, QFP, PFP);
	
}

// Calculate the perturbative correction
void ALMOSCF::perturb(bool order4)
{
	int nbfs = focker.getHCore().rows(); 
	int nocc = focker.getMolecule().getNel() / 2;
	int nvirt = nbfs - nocc;
	
	Matrix T = Eigen::MatrixXd::Zero(nbfs, nocc); 
	Matrix V = Eigen::MatrixXd::Zero(nbfs, nvirt);
	int row_offset = 0; int occ_col_offset = 0; int virt_col_offset = 0;
	for (auto f : fragments) {
		Matrix& f_cp = f.getCP(); 
		int f_nocc = f.getMolecule().getNel() / 2;
		int f_nbfs = f.getHCore().rows(); 
		int f_nvirt = f_nbfs - f_nocc; 
		
		T.block(row_offset, occ_col_offset, f_nbfs, f_nocc) = f_cp.block(0, 0, f_nbfs, f_nocc); 
		V.block(row_offset, virt_col_offset, f_nbfs, f_nvirt) = f_cp.block(0, f_nocc, f_nbfs, f_nvirt);
		
		row_offset += f_nbfs; 
		occ_col_offset += f_nocc; 
		virt_col_offset += f_nvirt; 
	}
	
	// Form MO overlap matrix
	Matrix& S = focker.getS();
	Matrix sigma(nbfs, nbfs);
	sigma.block(0, 0, nocc, nocc) = T.transpose() * S * T; 
	sigma.block(nocc, 0, nvirt, nocc) = V.transpose() * S * T;
	sigma.block(nocc, nocc, nvirt, nvirt) = V.transpose() * S * V;
	 
	// Cholesky decompose
	LLT llt(sigma); 
	Matrix C0 = llt.matrixL();
	C0 = C0.inverse();
	
	Matrix CMI(nbfs, nbfs);
	CMI.block(0, 0, nbfs, nocc) = T;
	CMI.block(0, nocc, nbfs, nvirt) = V;
	C0 = CMI * C0.transpose(); 
	
	Matrix& F = focker.getFockAO();
	sigma = C0.transpose() * F * C0; 
	Matrix F0occ = sigma.block(0, 0, nocc, nocc);
	Matrix F0virt = sigma.block(nocc, nocc, nvirt, nvirt); 
	Matrix Fov = sigma.block(0, nocc, nocc, nvirt); 
	EigenSolver es_occ(F0occ);
	EigenSolver es_virt(F0virt);
	T = es_occ.eigenvectors();
	V = es_virt.eigenvectors(); 
	Fov = T.transpose() * Fov * V; 
	Vector eps_i = es_occ.eigenvalues(); 
	Vector eps_a = es_virt.eigenvalues();
	
	if (order4) {
		e_pert_2 = 0.0; 
		e_pert_4 = 0.0;
		
		Matrix Fvo = V.transpose() * sigma.block(nocc, 0, nvirt, nocc) * T;
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

