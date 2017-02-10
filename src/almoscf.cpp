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
ALMOSCF::ALMOSCF(Command& c, SharedMolecule m, Fock& f) : molecule(m), cmd(c), focker(f) 
{	
	// Zero out energies
	dimer_energy = e_frz = e_pol = e_ct = e_int = 0.0;
	MAX = cmd.get_option<int>("maxdiis");
	diis.init(MAX, cmd.get_option<bool>("diis"));
}

void ALMOSCF::setFragments(bool unrestricted)
{
	// Add the fragment Fock matrices and perform monomer calculations
	std::vector<Fragment>& frags = molecule->getFragments(); 
	nfrags = frags.size();
	int start = 0;
	int nf = 0;
	for (int i = 0; i < frags.size(); i++) {
		
		frags[i].buildShellBasis();
		frags[i].calcEnuc();
		std::vector<libint2::Shell>& shells = frags[i].getBasis().getIntShells();

		int f_nbfs = ints[nf].nbasis(shells);
		IntegralEngine new_engine(frags[i].shared_from_this(), focker.getIntegrals(), start, start+f_nbfs);
		ints.push_back(new_engine);
		if (unrestricted) {
			ufragments.push_back(UnrestrictedFockFragment(cmd, ints[nf], frags[i].shared_from_this(), start, start+f_nbfs)); 
			SCF hf(cmd, frags[i].shared_from_this(), ufragments[ufragments.size()-1]);
			hf.uhf_internal(false, ufragments[ufragments.size()-1]);
			molecule->control->log.print("Monomer " + std::to_string(nf) + " energy = " + std::to_string(hf.getEnergy()) + " Hartree");
			molecule->control->log.flush();   
			monomer_energies.push_back(hf.getEnergy()); 
			ufragments[ufragments.size()-1].clearDiis(); 
		} else {
			fragments.push_back(FockFragment(cmd, ints[nf], frags[i].shared_from_this(), start, start+f_nbfs));
			SCF hf(cmd, frags[i].shared_from_this(), fragments[fragments.size()-1]);
			hf.rhf(false);
			molecule->control->log.print("Monomer " + std::to_string(nf) + " energy = " + std::to_string(hf.getEnergy()) + " Hartree");
			molecule->control->log.flush();  
			monomer_energies.push_back(hf.getEnergy()); 
			fragments[fragments.size()-1].clearDiis(); 
		}
		start += f_nbfs;

		nf++;
	}
}

double ALMOSCF::makeDens(bool alpha) {
// Make inverse overlap metric
	Matrix& S = focker.getS();
	Matrix P_old; int nocc;
	
	if (alpha) {
		P_old = P_alpha; 
		nocc = focker.getMolecule()->nalpha();
	} else {
		P_old = P_beta;
		nocc = focker.getMolecule()->nbeta(); 
	}
	
	//Build T matrix
	int nbfs = S.rows();
	int i_offset = 0; int mu_offset = 0; 
	Matrix sigma(nocc, nocc);
	int f1_nocc, f2_nocc; 
	for (int i = 0; i < ufragments.size(); i++) {
		auto& f1 = ufragments[i];
		if (alpha) f1_nocc = f1.getMolecule()->nalpha();
		else f1_nocc = f1.getMolecule()->nbeta();
		int f1_nbfs = f1.getHCore().rows(); 
		
		int j_offset = 0; int nu_offset = 0;
		for (int j = 0; j <= i; j++) {
			auto& f2 = ufragments[j];
			int f2_nbfs = f2.getHCore().rows(); 
			
			if (alpha) {
				f2_nocc = f2.getMolecule()->nalpha();
				sigma.block(i_offset, j_offset, f1_nocc, f2_nocc) = f1.getCPAlpha().block(0, 0, f1_nbfs, f1_nocc).transpose() 
					* S.block(mu_offset, nu_offset, f1_nbfs, f2_nbfs) * f2.getCPAlpha().block(0, 0, f2_nbfs, f2_nocc);
			} else {
				f2_nocc = f2.getMolecule()->nbeta();
				sigma.block(i_offset, j_offset, f1_nocc, f2_nocc) = f1.getCPBeta().block(0, 0, f1_nbfs, f1_nocc).transpose() 
					* S.block(mu_offset, nu_offset, f1_nbfs, f2_nbfs) * f2.getCPBeta().block(0, 0, f2_nbfs, f2_nocc);
			}
			
			j_offset += f2_nocc;
			nu_offset += f2_nbfs;
		}
		
		mu_offset += f1_nbfs; 
		i_offset += f1_nocc; 
	}
	sigma = sigma.selfadjointView<Eigen::Lower>();
	sigma = sigma.inverse(); 
	
	i_offset = 0; mu_offset = 0;
	for (int i = 0; i < ufragments.size(); i++) {
		auto& f1 = ufragments[i];
		if (alpha) f1_nocc = f1.getMolecule()->nalpha();
		else f1_nocc = f1.getMolecule()->nbeta();		
		int f1_nbfs = f1.getHCore().rows(); 
		
		int j_offset = 0; int nu_offset = 0;
		for (int j = 0; j <= i; j++) {
			auto& f2 = ufragments[j];
			int f2_nbfs = f2.getHCore().rows(); 
			
			if (alpha) {
				f2_nocc = f2.getMolecule()->nalpha();
				P_alpha.block(mu_offset, nu_offset, f1_nbfs, f2_nbfs) = f1.getCPAlpha().block(0, 0, f1_nbfs, f1_nocc) 
					* sigma.block(i_offset, j_offset, f1_nocc, f2_nocc) * f2.getCPAlpha().block(0, 0, f2_nbfs, f2_nocc).transpose();
			} else {
				f2_nocc = f2.getMolecule()->nbeta();
				P_beta.block(mu_offset, nu_offset, f1_nbfs, f2_nbfs) = f1.getCPBeta().block(0, 0, f1_nbfs, f1_nocc) 
					* sigma.block(i_offset, j_offset, f1_nocc, f2_nocc) * f2.getCPBeta().block(0, 0, f2_nbfs, f2_nocc).transpose();
			}
			
			j_offset += f2_nocc;
			nu_offset += f2_nbfs;
		}
		
		mu_offset += f1_nbfs; 
		i_offset += f1_nocc; 
	}
	
	double dd; 
	if (alpha) {
		P_alpha = P_alpha.selfadjointView<Eigen::Lower>();
		dd = (P_alpha - P_old).norm(); 
	} else {
		P_beta = P_beta.selfadjointView<Eigen::Lower>();
		dd = (P_beta - P_old).norm(); 
	}
	
	return dd; 
}

// Routines
void ALMOSCF::rcompute() {
	
	// Make inverse overlap metric
	Matrix& S = focker.getS();
	Matrix& H = focker.getHCore(); 
	Matrix P_old = P;
	int nocc = focker.getMolecule()->getNel() / 2;
	
	//Build T matrix
	int nbfs = S.rows();
	int i_offset = 0; int mu_offset = 0; 
	Matrix sigma(nocc, nocc);
	int f1_nocc, f2_nocc; 
	for (int i = 0; i < fragments.size(); i++) {
		auto& f1 = fragments[i];
		f1_nocc = f1.getMolecule()->getNel()/2;
		int f1_nbfs = f1.getHCore().rows(); 
		
		int j_offset = 0; int nu_offset = 0;
		for (int j = 0; j <= i; j++) {
			auto& f2 = fragments[j];
			int f2_nbfs = f2.getHCore().rows(); 
			f2_nocc = f2.getMolecule()->getNel()/2;
			sigma.block(i_offset, j_offset, f1_nocc, f2_nocc) = f1.getCP().block(0, 0, f1_nbfs, f1_nocc).transpose() 
					* S.block(mu_offset, nu_offset, f1_nbfs, f2_nbfs) * f2.getCP().block(0, 0, f2_nbfs, f2_nocc);

			j_offset += f2_nocc;
			nu_offset += f2_nbfs;
		}
		
		mu_offset += f1_nbfs; 
		i_offset += f1_nocc; 
	}
	sigma = sigma.selfadjointView<Eigen::Lower>();
	sigma = sigma.inverse(); 
	
	i_offset = 0; mu_offset = 0;
	for (int i = 0; i < fragments.size(); i++) {
		auto& f1 = fragments[i];
		f1_nocc = f1.getMolecule()->getNel()/2;
		int f1_nbfs = f1.getHCore().rows(); 
		
		int j_offset = 0; int nu_offset = 0;
		for (int j = 0; j <= i; j++) {
			auto& f2 = fragments[j];
			int f2_nbfs = f2.getHCore().rows(); 
			f2_nocc = f2.getMolecule()->getNel()/2;
			P.block(mu_offset, nu_offset, f1_nbfs, f2_nbfs) = f1.getCP().block(0, 0, f1_nbfs, f1_nocc) 
					* sigma.block(i_offset, j_offset, f1_nocc, f2_nocc) * f2.getCP().block(0, 0, f2_nbfs, f2_nocc).transpose();
			
			j_offset += f2_nocc;
			nu_offset += f2_nbfs;
		}
		
		mu_offset += f1_nbfs; 
		i_offset += f1_nocc; 
	}
	P = P.selfadjointView<Eigen::Lower>();
	delta_d = (P - P_old).norm(); 
	
	// Build Fock matrix
	if (focker.getMolecule()->control->get_option<bool>("direct"))
		focker.formJKdirect(focker.getIntegrals().getPrescreen(), P, 2.0);
	else 
		focker.formJK(P, 2.0);
	
	focker.makeFock();
	Matrix& F = focker.getFockAO(); 
	
	// Calculate errors
	Matrix Q = - S * P; 
	for (int i = 0; i < nbfs; i++) Q(i, i) += 1.0; 
	Matrix QFP = Q * F; 
	
	double new_dimer_energy = (P * (H + F)).trace() + focker.getMolecule()->getEnuc();
	delta_e = new_dimer_energy - dimer_energy;
	dimer_energy = new_dimer_energy;
	
	Matrix QFQ = QFP * Q.transpose();
	QFP = QFP * P; 
	Matrix PFP = P * F * P; 
	
	std::vector<Vector> errs;
	for (int i = 0; i < fragments.size(); i++)
		errs.push_back(fragments[i].buildFock(QFQ, QFP, PFP));
	
	if (cmd.get_option<bool>("diis")) {
		Vector weights = diis.compute(errs); 
		for(int i = 0; i < fragments.size(); i++)
			fragments[i].average(weights);
	}
	
	for (int i = 0; i < fragments.size(); i++)
		fragments[i].gensolve(); 
}

void ALMOSCF::ucompute() {
	UnrestrictedFock& ufocker = dynamic_cast<UnrestrictedFock&>(focker);
	
	Matrix& S = ufocker.getS();
	Matrix& H = ufocker.getHCore();
	int nbfs = S.rows();
	
	delta_d = makeDens(true);
	delta_d += makeDens(false); 
	
	// Build Fock matrix
	if (ufocker.getMolecule()->control->get_option<bool>("direct"))
		ufocker.formJKdirect(focker.getIntegrals().getPrescreen(), P_alpha, P_beta, 2.0);
	else 
		ufocker.formJK(P_alpha, P_beta, 2.0);
	
	ufocker.makeFock();
	Matrix& Fa = ufocker.getFockAlphaAO(); 
	Matrix& Fb = ufocker.getFockBetaAO();
	
	// Calculate errors
	Matrix Q_alpha = - S * P_alpha;
	Matrix Q_beta = -S * P_beta;  
	for (int i = 0; i < nbfs; i++) {
		Q_alpha(i, i) += 1.0; 
		Q_beta(i, i) += 1.0; 
	}
	Matrix QFP_alpha = Q_alpha * Fa;
	Matrix QFP_beta = Q_beta * Fb; 
	
	double new_dimer_energy = 0.5 * (P_alpha * (H + Fa) + P_beta * (H + Fb)).trace() + ufocker.getMolecule()->getEnuc();
	delta_e = new_dimer_energy - dimer_energy;
	dimer_energy = new_dimer_energy;
	
	Matrix QFQ_alpha = QFP_alpha * Q_alpha.transpose();
	Matrix QFQ_beta = QFP_beta * Q_beta.transpose();
	QFP_alpha = QFP_alpha * P_alpha; 
	QFP_beta = QFP_beta * P_beta;
	Matrix PFP_alpha = P_alpha * Fa * P_alpha; 
	Matrix PFP_beta = P_beta * Fb * P_beta;
	
	std::vector<Vector> errs;
	for (int i = 0; i < ufragments.size(); i++) {
		errs.push_back(ufragments[i].buildFock(QFQ_alpha, QFP_alpha, PFP_alpha, true));
		errs.push_back(ufragments[i].buildFock(QFQ_beta, QFP_beta, PFP_beta, false)); 
	}
	
	if (cmd.get_option<bool>("diis")) {
		Vector weights = diis.compute(errs); 
		for(int i = 0; i < ufragments.size(); i++)
			ufragments[i].average(weights);
	}
	
	for (int i = 0; i < ufragments.size(); i++)
		ufragments[i].gensolve(); 
}

// Calculate the perturbative correction
void ALMOSCF::rperturb(bool order4)
{
	int nbfs = focker.getHCore().rows(); 
	int nocc = focker.getMolecule()->getNel() / 2;
	int nvirt = nbfs - nocc;
	
	Matrix T = Eigen::MatrixXd::Zero(nbfs, nocc); 
	Matrix V = Eigen::MatrixXd::Zero(nbfs, nvirt);
	int row_offset = 0; int occ_col_offset = 0; int virt_col_offset = 0;
	for (auto f : fragments) {
		Matrix& f_cp = f.getCP(); 
		int f_nocc = f.getMolecule()->getNel() / 2;
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
	P = Matrix::Zero(focker.getHCore().rows(), focker.getHCore().rows());
	molecule->control->log.title("Closed-shell ALMO Calculation");
	
	setFragments();
	
	molecule->control->log.initIteration(); 
	delta_d = 1.0; delta_e = 1.0; 
	bool converged = false;
	int iter = 0;
	double CONVERGE = cmd.get_option<double>("converge");
	int MAXITER = cmd.get_option<int>("maxiter");
	
	while (!converged && iter < MAXITER) {
		rcompute(); 
		molecule->control->log.iteration(iter++, dimer_energy, delta_e, delta_d);
		converged = (fabs(delta_e) < CONVERGE) && (fabs(delta_d) < CONVERGE);
	}
	
	if (converged) {
		
		double energy = dimer_energy;
		for (auto en : monomer_energies) energy -= en; 
		molecule->control->log.result("ALMO Interaction Energy", energy * Logger::TOKCAL, "kcal / mol"); 
		
		int perturbation = cmd.get_option<int>("perturb"); 
		if (perturbation > 0) {
			bool fourth = perturbation == 2; 
			rperturb(fourth); 
			e_pert_2 *= 2.0;
			molecule->control->log.result("E(2)", e_pert_2 * Logger::TOKCAL, "kcal /mol"); 
			
			if (fourth) {
				e_pert_4 *= 2.0;
				molecule->control->log.result("E(4)", e_pert_4 * Logger::TOKCAL, "kcal / mol"); 
			} else {
				e_pert_4 = 0.0;
			}
			molecule->control->log.result("Total ALMO+CT interaction energy", (energy + e_pert_2 + e_pert_4) * Logger::TOKCAL, "kcal / mol");
		}
		
	} else {
		molecule->control->log.result("ALMO SCF failed to converge");
	}
}

void ALMOSCF::uscf()
{
	int nbfs = focker.getHCore().rows();
	P_alpha = Matrix::Zero(nbfs, nbfs);
	P_beta = Matrix::Zero(nbfs, nbfs);
	molecule->control->log.title("Unrestricted Open-shell ALMO Calculation"); 
	
	setFragments(true); 
	
	molecule->control->log.initIteration(); 
	delta_d = 1.0; delta_e = 1.0; 
	bool converged = false;
	int iter =0 ;
	double CONVERGE = cmd.get_option<double>("converge");
	int MAXITER = cmd.get_option<int>("maxiter");
	
	while (!converged && iter < MAXITER) {
		ucompute(); 
		molecule->control->log.iteration(iter++, dimer_energy, delta_e, delta_d);
		converged = (fabs(delta_e) < CONVERGE) && (fabs(delta_d) < CONVERGE);
	}
	
	if (converged) {
		double energy = dimer_energy;
		for (auto en : monomer_energies) energy -= en; 
		molecule->control->log.result("ALMO Interaction Energy", energy * Logger::TOKCAL, "kcal / mol"); 
	} else {
		molecule->control->log.result("ALMO SCF failed to converge");
	}
	
}

