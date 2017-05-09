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
#include "rpa.hpp"
#include "cc.hpp"

// Constructor
ALMOSCF::ALMOSCF(Command& c, SharedMolecule m, Fock& f) : molecule(m), cmd(c), focker(f) 
{	
	// Zero out energies
	dimer_energy = e_frz = e_pol = e_ct = e_int = e_disp = e_mon_rpa = 0.0;
	MAX = cmd.get_option<int>("maxdiis");
	diis.init(MAX, cmd.get_option<bool>("diis"));
	focker.incremental_threshold = 0.0; 
}

void ALMOSCF::setFragments(bool unrestricted)
{
	// Add the fragment Fock matrices and perform monomer calculations
	std::vector<SharedFragment>& frags = molecule->getFragments(); 
	nfrags = frags.size();
	int start = 0;
	int nf = 0;
	for (int i = 0; i < frags.size(); i++) {
		
		frags[i]->buildShellBasis();
		frags[i]->calcEnuc();
		std::vector<libint2::Shell>& shells = frags[i]->getBasis().getIntShells();

		int f_nbfs = ints[nf].nbasis(shells);
		IntegralEngine new_engine(frags[i]->shared_from_this(), focker.getIntegrals(), start, start+f_nbfs);
		ints.push_back(new_engine);
		if (unrestricted) {
			ufragments.push_back(UnrestrictedFockFragment(cmd, ints[nf], frags[i], start, start+f_nbfs)); 
			SCF hf(cmd, frags[i], ufragments[ufragments.size()-1]);
			hf.uhf_internal(false, ufragments[ufragments.size()-1]);
			molecule->control->log.print("Monomer " + std::to_string(nf) + " energy = " + std::to_string(hf.getEnergy()) + " Hartree");
			molecule->control->log.flush();   
			monomer_energies.push_back(hf.getEnergy()); 
			ufragments[ufragments.size()-1].clearDiis(); 
		} else {
			fragments.push_back(FockFragment(cmd, ints[nf], frags[i], start, start+f_nbfs));
			SCF hf(cmd, frags[i], fragments[fragments.size()-1]);
			hf.rhf(false);
			molecule->control->log.print("Monomer " + std::to_string(nf+1) + " energy = " + std::to_string(hf.getEnergy()) + " Hartree");
			
			molecule->control->log.flush();  
			monomer_energies.push_back(hf.getEnergy()); 
			fragments[fragments.size()-1].clearDiis(); 
		}
		start += f_nbfs;

		nf++;
	}
	molecule->control->log.print("Monomer calculations completed.");
	molecule->control->log.localTime();
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
	sigma = Matrix::Zero(nocc, nocc);
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
	
	if (alpha) sigma_alpha = sigma;
	else sigma_beta = sigma; 
	
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
	Matrix T; 
	bool density_fitted = cmd.get_option<bool>("df"); 
	if (density_fitted) T = Matrix::Zero(H.rows(), focker.getMolecule()->getNel()/2); 
	int nocc = focker.getMolecule()->getNel() / 2;
	
	//Build T matrix
	int nbfs = S.rows();
	int i_offset = 0; int mu_offset = 0; 
	sigma = Matrix::Zero(nocc, nocc);
	int f1_nocc, f2_nocc; 
	for (int i = 0; i < fragments.size(); i++) {
		auto& f1 = fragments[i];
		f1_nocc = f1.getMolecule()->getNel()/2;
		int f1_nbfs = f1.getHCore().rows(); 
		
		if (density_fitted) T.block(mu_offset, i_offset, f1_nbfs, f1_nocc) = f1.getCP().block(0, 0, f1_nbfs, f1_nocc); 
		
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
	if (density_fitted) {
		EigenSolver es(sigma); 
		T = T * es.operatorSqrt(); 
		focker.makeFock(T, 1.0); 
	} else
		focker.makeFock(P, 2.0); 
	
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
	if (cmd.get_option<bool>("df")) {
		Matrix ca_occ = Matrix::Zero(nbfs, focker.getMolecule()->nalpha()); 
		Matrix cb_occ = Matrix::Zero(nbfs, focker.getMolecule()->nbeta()); 
		
		int ia_offset = 0; int ib_offset = 0; int mu_offset = 0; 
		int f1_nalpha, f1_nbeta; 
		for (int i = 0; i < ufragments.size(); i++) {
			auto& f1 = ufragments[i];
			f1_nalpha = f1.getMolecule()->nalpha();
			f1_nbeta = f1.getMolecule()->nbeta(); 
			int f1_nbfs = f1.getHCore().rows(); 
			
			ca_occ.block(mu_offset, ia_offset, f1_nbfs, f1_nalpha) = f1.getCPAlpha().block(0, 0, f1_nbfs, f1_nalpha); 
			cb_occ.block(mu_offset, ib_offset, f1_nbfs, f1_nbeta) = f1.getCPBeta().block(0, 0, f1_nbfs, f1_nbeta); 
		
			mu_offset += f1_nbfs; 
			ia_offset += f1_nalpha;
			ib_offset += f1_nbeta; 
		}
		
		EigenSolver esa(sigma_alpha); 
		ca_occ = ca_occ * esa.operatorSqrt(); 
		EigenSolver esb(sigma_beta);
		cb_occ = cb_occ * esb.operatorSqrt(); 
		
		ufocker.formJKdf(ca_occ, cb_occ, 1.0); 
	} else if (ufocker.getMolecule()->control->get_option<bool>("direct"))
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
	for (auto& f : fragments) {
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
	sigma = Matrix::Zero(nbfs, nbfs);
	sigma.block(0, 0, nocc, nocc) = T.transpose() * S * T; 
	sigma.block(nocc, 0, nvirt, nocc) = V.transpose() * S * T;
	sigma.block(nocc, nocc, nvirt, nvirt) = V.transpose() * S * V;
	 
	EigenSolver es(S); 
	Matrix SC = es.operatorSqrt() * T; 
	int i_offset = 0;
	for (auto& f1 : fragments) {
		int f1_nocc = f1.getMolecule()->getNel() / 2;
		
		int mu_offset = 0; int center = 1; 
		for (auto& f2 : fragments) {
			int f2_nbfs = f2.getHCore().rows(); 
			
			for (int i = i_offset; i < i_offset + f1_nocc; i++) {
				double INi = 0.0; 
				for (int mu = mu_offset; mu < mu_offset + f2_nbfs; mu++) {
					INi += SC(mu, i) * SC(mu, i);	
				}
				std::cout << "Orbital " << i << " and center " << center << ": " << INi << std::endl; 
				
			}
			center++; 
			
			mu_offset += f2_nbfs; 
		}
		i_offset += f1_nocc; 
	}
	 
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

// Calculate the perturbative correction
void ALMOSCF::uperturb(bool order4)
{
	UnrestrictedFock& ufocker = dynamic_cast<UnrestrictedFock&>(focker);
	
	int nbfs = ufocker.getHCore().rows(); 
	int nalpha = ufocker.getMolecule()->nalpha();
	int nbeta = ufocker.getMolecule()->nbeta();
	int nvirt_alpha = nbfs - nalpha;
	int nvirt_beta = nbfs - nbeta; 
	
	Matrix T_alpha = Eigen::MatrixXd::Zero(nbfs, nalpha); 
	Matrix V_alpha = Eigen::MatrixXd::Zero(nbfs, nvirt_alpha);
	Matrix T_beta = Eigen::MatrixXd::Zero(nbfs, nbeta); 
	Matrix V_beta = Eigen::MatrixXd::Zero(nbfs, nvirt_beta);
	int row_offset = 0; 
	int occ_col_offset_alpha = 0; int virt_col_offset_alpha = 0;
	int occ_col_offset_beta = 0; int virt_col_offset_beta = 0;
	for (auto& f : ufragments) {
		Matrix& f_cp_alpha = f.getCPAlpha(); 
		Matrix& f_cp_beta = f.getCPBeta(); 
		int f_nalpha = f.getMolecule()->nalpha();
		int f_nbeta = f.getMolecule()->nbeta();
		int f_nbfs = f.getHCore().rows(); 
		int f_nvirt_alpha = f_nbfs - f_nalpha; 
		int f_nvirt_beta = f_nbfs - f_nbeta; 
		
		T_alpha.block(row_offset, occ_col_offset_alpha, f_nbfs, f_nalpha) = f_cp_alpha.block(0, 0, f_nbfs, f_nalpha); 
		V_alpha.block(row_offset, virt_col_offset_alpha, f_nbfs, f_nvirt_alpha) = f_cp_alpha.block(0, f_nalpha, f_nbfs, f_nvirt_alpha);
		T_beta.block(row_offset, occ_col_offset_beta, f_nbfs, f_nbeta) = f_cp_beta.block(0, 0, f_nbfs, f_nbeta); 
		V_beta.block(row_offset, virt_col_offset_beta, f_nbfs, f_nvirt_beta) = f_cp_beta.block(0, f_nbeta, f_nbfs, f_nvirt_beta);
		
		row_offset += f_nbfs; 
		occ_col_offset_alpha += f_nalpha; 
		virt_col_offset_alpha += f_nvirt_alpha; 
		occ_col_offset_beta += f_nbeta; 
		virt_col_offset_beta += f_nvirt_beta; 
	}
	
	// Form MO overlap matrix
	Matrix& S = ufocker.getS();
	sigma_alpha = Matrix::Zero(nbfs, nbfs);
	sigma_beta = Matrix::Zero(nbfs, nbfs);
	sigma_alpha.block(0, 0, nalpha, nalpha) = T_alpha.transpose() * S * T_alpha; 
	sigma_alpha.block(nalpha, 0, nvirt_alpha, nalpha) = V_alpha.transpose() * S * T_alpha;
	sigma_alpha.block(nalpha, nalpha, nvirt_alpha, nvirt_alpha) = V_alpha.transpose() * S * V_alpha;
	sigma_beta.block(0, 0, nbeta, nbeta) = T_beta.transpose() * S * T_beta; 
	sigma_beta.block(nbeta, 0, nvirt_beta, nbeta) = V_beta.transpose() * S * T_beta;
	sigma_beta.block(nbeta, nbeta, nvirt_beta, nvirt_beta) = V_beta.transpose() * S * V_beta;
	 
	// Cholesky decompose
	LLT llt_alpha(sigma_alpha);
	LLT llt_beta(sigma_beta); 
	Matrix C0_alpha = llt_alpha.matrixL();
	C0_alpha = C0_alpha.inverse();
	Matrix C0_beta = llt_beta.matrixL();
	C0_beta = C0_beta.inverse();
	
	Matrix CMI_alpha(nbfs, nbfs);
	CMI_alpha.block(0, 0, nbfs, nalpha) = T_alpha;
	CMI_alpha.block(0, nalpha, nbfs, nvirt_alpha) = V_alpha;
	C0_alpha = CMI_alpha * C0_alpha.transpose(); 
	Matrix CMI_beta(nbfs, nbfs);
	CMI_beta.block(0, 0, nbfs, nbeta) = T_beta;
	CMI_beta.block(0, nbeta, nbfs, nvirt_beta) = V_beta;
	C0_beta = CMI_beta * C0_beta.transpose(); 
	
	Matrix& F_alpha = ufocker.getFockAlphaAO();
	sigma_alpha = C0_alpha.transpose() * F_alpha * C0_alpha; 
	Matrix F0occ_alpha = sigma_alpha.block(0, 0, nalpha, nalpha);
	Matrix F0virt_alpha = sigma_alpha.block(nalpha, nalpha, nvirt_alpha, nvirt_alpha); 
	Matrix Fov_alpha = sigma_alpha.block(0, nalpha, nalpha, nvirt_alpha); 
	EigenSolver es_occ_alpha(F0occ_alpha);
	EigenSolver es_virt_alpha(F0virt_alpha);
	T_alpha = es_occ_alpha.eigenvectors();
	V_alpha = es_virt_alpha.eigenvectors(); 
	Fov_alpha = T_alpha.transpose() * Fov_alpha * V_alpha; 
	Vector eps_i_alpha = es_occ_alpha.eigenvalues(); 
	Vector eps_a_alpha = es_virt_alpha.eigenvalues();
	
	Matrix& F_beta = ufocker.getFockBetaAO();
	sigma_beta = C0_beta.transpose() * F_beta * C0_beta; 
	Matrix F0occ_beta = sigma_beta.block(0, 0, nbeta, nbeta);
	Matrix F0virt_beta = sigma_beta.block(nbeta, nbeta, nvirt_beta, nvirt_beta); 
	Matrix Fov_beta = sigma_beta.block(0, nbeta, nbeta, nvirt_beta); 
	EigenSolver es_occ_beta(F0occ_beta);
	EigenSolver es_virt_beta(F0virt_beta);
	T_beta = es_occ_beta.eigenvectors();
	V_beta = es_virt_beta.eigenvectors(); 
	Fov_beta = T_beta.transpose() * Fov_beta * V_beta; 
	Vector eps_i_beta = es_occ_beta.eigenvalues(); 
	Vector eps_a_beta = es_virt_beta.eigenvalues();
	
	if (order4) {
		e_pert_2 = 0.0; 
		e_pert_4 = 0.0;
		
		Matrix Fvo_alpha = V_alpha.transpose() * sigma_alpha.block(nalpha, 0, nvirt_alpha, nalpha) * T_alpha;
		for (int i = 0; i < nalpha; i++) {
			double sum1 = 0.0, sum2 = 0.0;
			for (int a = 0; a < nvirt_alpha; a++) {
				auto delta = eps_i_alpha(i) - eps_a_alpha(a);
				auto val = Fov_alpha(i, a) * Fov_alpha(i, a) / delta; 
				e_pert_2 += val;
				sum1 += val; 
				sum2 += val / delta;
				
				for (int b = 0; b < nvirt_alpha; b++) {
					for (int j = 0; j < nalpha; j++) {
						if (i == j) continue;
						e_pert_4 += Fov_alpha(i, a) * Fvo_alpha(a, j) * Fov_alpha(j, b) * Fvo_alpha(b, i) / ( delta * (eps_i_alpha(i) - eps_i_alpha(j)) * (eps_i_alpha(i) - eps_a_alpha(b)) ); 
					}
				}
			}
			e_pert_4 -= sum1 * sum2; 
		}
		
		Matrix Fvo_beta = V_beta.transpose() * sigma_beta.block(nbeta, 0, nvirt_beta, nbeta) * T_beta;
		for (int i = 0; i < nbeta; i++) {
			double sum1 = 0.0, sum2 = 0.0;
			for (int a = 0; a < nvirt_beta; a++) {
				auto delta = eps_i_beta(i) - eps_a_beta(a);
				auto val = Fov_beta(i, a) * Fov_beta(i, a) / delta; 
				e_pert_2 += val;
				sum1 += val; 
				sum2 += val / delta;
				
				for (int b = 0; b < nvirt_beta; b++) {
					for (int j = 0; j < nbeta; j++) {
						if (i == j) continue;
						e_pert_4 += Fov_beta(i, a) * Fvo_beta(a, j) * Fov_beta(j, b) * Fvo_beta(b, i) / ( delta * (eps_i_beta(i) - eps_i_beta(j)) * (eps_i_beta(i) - eps_a_beta(b)) ); 
					}
				}
			}
			e_pert_4 -= sum1 * sum2; 
		}
	} else {
		e_pert_2 = 0.0;
		for (int i = 0; i < nalpha; i++)
			for (int a = 0; a < nvirt_alpha; a++)
				e_pert_2 += Fov_alpha(i, a) * Fov_alpha(i, a) / (eps_i_alpha(i) - eps_a_alpha(a)); 
		for (int i = 0; i < nbeta; i++)
			for (int a = 0; a < nvirt_beta; a++)
				e_pert_2 += Fov_beta(i, a) * Fov_beta(i, a) / (eps_i_beta(i) - eps_a_beta(a)); 
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
	double E_CONVERGE = cmd.get_option<double>("enconverge");
	double D_CONVERGE = cmd.get_option<double>("densconverge");
	int MAXITER = cmd.get_option<int>("maxiter");
	
	while (!converged && iter < MAXITER) {
		rcompute(); 
		molecule->control->log.iteration(iter++, dimer_energy, delta_e, delta_d);
		converged = (fabs(delta_e) < E_CONVERGE) && (fabs(delta_d) < D_CONVERGE);
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
		if (cmd.get_option<bool>("rpa")) {
			/*focker.transform();
			focker.diagonalise(); 
			focker.makeDens();
			
			Matrix& Pinf = focker.getDens(); 
			Matrix& H = focker.getHCore();
			Matrix& F = focker.getFockAO(); 
			
			double e_inf =(F * (Pinf - 2.0*P)).trace();
			molecule->control->log.result("E(Inf.)", e_inf * Logger::TOKCAL, "kcal /mol");*/
			
			double e_inf = e_pert_2;
	
			if (cmd.get_option<bool>("pairwise")) {
				
				double rcut = cmd.get_option<double>("rcutoff"); 
				molecule->control->log.title("Pairwise RPAxd"); 
				molecule->control->log.print("Using COM distance cutoff of " + std::to_string(rcut) + " Bohr"); 
				
				molecule->control->log.initALMOTable(); 
				
				int mu_offset = 0; 
				int f1_nocc, f2_nocc, f1_nbfs, f2_nbfs, f1_nvirt, f2_nvirt; 
				double sep, edisp, edispexch;
				Vector com1, com2;  
				for (int n1 = 0; n1 < nfrags; n1++) {
					auto& f1 = fragments[n1]; 
					f1_nocc = f1.getMolecule()->getNel()/2; 
					f1_nbfs = f1.getHCore().rows(); 
					f1_nvirt = f1_nbfs - f1_nocc; 
					 
					com1 = f1.getMolecule()->com(); 
					std::vector<libint2::Shell>& f1obs = f1.getMolecule()->getBasis().getIntShells();
					std::vector<libint2::Shell>& f1auxbs = f1.getMolecule()->getBasis().getDFShells();
					
					int nu_offset = mu_offset + f1_nbfs;
					for (int n2 = n1+1; n2 < nfrags; n2++) {
						auto& f2 = fragments[n2]; 
						f2_nocc = f2.getMolecule()->getNel()/2; 
						f2_nbfs = f2.getHCore().rows(); 
						f2_nvirt = f2_nbfs - f2_nocc; 
						
						std::vector<libint2::Shell>& f2obs = f2.getMolecule()->getBasis().getIntShells();
						std::vector<libint2::Shell>& f2auxbs = f2.getMolecule()->getBasis().getDFShells();
						com2 = f2.getMolecule()->com(); 
						
						sep = (com1 - com2).norm(); 
						edisp = edispexch = 0.0; 
						
						if (sep < rcut) {
							
							std::vector<libint2::Shell> obs = f1obs; 
							obs.insert(obs.end(), f2obs.begin(), f2obs.end()); 
							std::vector<libint2::Shell> auxbs = f1auxbs; 
							auxbs.insert(auxbs.end(), f2auxbs.begin(), f2auxbs.end()); 
							
							fInfo info(obs, auxbs); 
							
							int nbfs = f1_nbfs + f2_nbfs; 
							int nocc = f1_nocc + f2_nocc; 
							int nvirt = nbfs - nocc; 
							
							info.nocc.push_back(f1_nocc);
							info.nocc.push_back(f2_nocc);
							info.nvirt.push_back(f1_nvirt);
							info.nvirt.push_back(f2_nvirt);
							
							info.T = Matrix::Zero(nbfs, nocc);
							info.V = Matrix::Zero(nbfs, nvirt);
							info.S = Matrix::Zero(nbfs, nbfs);
							info.F = Matrix::Zero(nbfs, nbfs);  
							
							info.T.block(0, 0, f1_nbfs, f1_nocc) = f1.getCP().block(0, 0, f1_nbfs, f1_nocc); 
							info.V.block(0, 0, f1_nbfs, f1_nvirt) = f1.getCP().block(0, f1_nocc, f1_nbfs, f1_nvirt); 
							info.T.block(f1_nbfs, f1_nocc, f2_nbfs, f2_nocc) = f2.getCP().block(0, 0, f2_nbfs, f2_nocc);
							info.V.block(f1_nbfs, f1_nvirt, f2_nbfs, f2_nvirt) = f2.getCP().block(0, f2_nocc, f2_nbfs, f2_nvirt); 
							
							info.S.block(0, 0, f1_nbfs, f1_nbfs) = focker.getS().block(mu_offset, mu_offset, f1_nbfs, f1_nbfs); 
							info.S.block(f1_nbfs, f1_nbfs, f2_nbfs, f2_nbfs) = focker.getS().block(nu_offset, nu_offset, f2_nbfs, f2_nbfs); 
							info.S.block(0, f1_nbfs, f1_nbfs, f2_nbfs) = focker.getS().block(mu_offset, nu_offset, f1_nbfs, f2_nbfs); 
							info.S.block(f1_nbfs, 0, f2_nbfs, f1_nbfs) = info.S.block(0, f1_nbfs, f1_nbfs, f2_nbfs).transpose(); 
							
							info.F.block(0, 0, f1_nbfs, f1_nbfs) = focker.getFockAO().block(mu_offset, mu_offset, f1_nbfs, f1_nbfs); 
							info.F.block(f1_nbfs, f1_nbfs, f2_nbfs, f2_nbfs) = focker.getFockAO().block(nu_offset, nu_offset, f2_nbfs, f2_nbfs); 
							info.F.block(0, f1_nbfs, f1_nbfs, f2_nbfs) = focker.getFockAO().block(mu_offset, nu_offset, f1_nbfs, f2_nbfs); 
							info.F.block(f1_nbfs, 0, f2_nbfs, f1_nbfs) = info.F.block(0, f1_nbfs, f1_nbfs, f2_nbfs).transpose(); 
							
							RPA rpa(cmd, focker, nbfs, nocc); 
							rpa.fcompute(info, false); 
							
							edisp = info.edisp;
							edispexch = info.edispexch; 
							e_disp += edisp + edispexch; 
							
						}
						
						molecule->control->log.ALMORow(n1, n2, sep, edisp, edispexch); 
						
						nu_offset += f2_nbfs; 
					}

					mu_offset += f1_nbfs; 
				}
				
				
			} else { 
				int nbfs = focker.getHCore().rows(); 
				int nocc = focker.getMolecule()->getNel() / 2;
				int nvirt = nbfs - nocc;
			
				fInfo info(focker.getMolecule()->getBasis().getIntShells(), focker.getMolecule()->getBasis().getDFShells()); 
				info.T = Eigen::MatrixXd::Zero(nbfs, nocc); 
				info.V = Eigen::MatrixXd::Zero(nbfs, nvirt);
				int row_offset = 0; int occ_col_offset = 0; int virt_col_offset = 0;
				for (auto& f : fragments) {
					Matrix& f_cp = f.getCP(); 
					int f_nocc = f.getMolecule()->getNel() / 2;
					int f_nbfs = f.getHCore().rows(); 
					int f_nvirt = f_nbfs - f_nocc; 
				
					info.nocc.push_back(f_nocc);
					info.nvirt.push_back(f_nvirt); 
		
					info.T.block(row_offset, occ_col_offset, f_nbfs, f_nocc) = f_cp.block(0, 0, f_nbfs, f_nocc); 
					info.V.block(row_offset, virt_col_offset, f_nbfs, f_nvirt) = f_cp.block(0, f_nocc, f_nbfs, f_nvirt);
		
					row_offset += f_nbfs; 
					occ_col_offset += f_nocc; 
					virt_col_offset += f_nvirt; 
				}
				info.S = focker.getS();
				info.F = focker.getFockAO(); 
			
				RPA rpa(cmd, focker, nbfs, nocc); 
				rpa.fcompute(info, true); 
				e_disp = rpa.getEnergy();// - e_mon_rpa; 
			
			}
		
			molecule->control->log.result("E(Disp.)",  e_disp * Logger::TOKCAL, "kcal /mol"); 
			
			molecule->control->log.result("Total ALMO+RPA interaction energy", (energy + e_disp + e_inf) * Logger::TOKCAL, "kcal / mol"); 
	
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
	double E_CONVERGE = cmd.get_option<double>("enconverge");
	double D_CONVERGE = cmd.get_option<double>("densconverge");
	int MAXITER = cmd.get_option<int>("maxiter");
	
	while (!converged && iter < MAXITER) {
		ucompute(); 
		molecule->control->log.iteration(iter++, dimer_energy, delta_e, delta_d);
		converged = (fabs(delta_e) < E_CONVERGE) && (fabs(delta_d) < D_CONVERGE);
	}
	
	if (converged) {
		double energy = dimer_energy;
		for (auto en : monomer_energies) energy -= en; 
		molecule->control->log.result("ALMO Interaction Energy", energy * Logger::TOKCAL, "kcal / mol"); 
		
		int perturbation = cmd.get_option<int>("perturb"); 
		if (perturbation > 0) {
			bool fourth = perturbation == 2; 
			uperturb(fourth); 
			molecule->control->log.result("E(2)", e_pert_2 * Logger::TOKCAL, "kcal /mol"); 
			
			if (fourth)
				molecule->control->log.result("E(4)", e_pert_4 * Logger::TOKCAL, "kcal / mol"); 
			else
				e_pert_4 = 0.0;
			molecule->control->log.result("Total ALMO+CT interaction energy", (energy + e_pert_2 + e_pert_4) * Logger::TOKCAL, "kcal / mol");
		}
	} else {
		molecule->control->log.result("ALMO SCF failed to converge");
	}
	
}

