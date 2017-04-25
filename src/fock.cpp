/*
*
*   PURPOSE: To implement class Fock, a class containing the data and routines
*            needed for Fock space based methods.
*
*   DATE          AUTHOR              CHANGES
*  ===========================================================================
*  14/09/15       Robert Shaw         Original code.
* 
*/

//#ifndef EIGEN_USE_MKL_ALL
//#define EIGEN_USE_MKL_ALL
//#endif

#include "fock.hpp"
#include "error.hpp"
#include "atom.hpp"
#include "tensor4.hpp"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <libint2.hpp>
#include "logger.hpp"
#include "ProgramController.hpp"

// Constructor
Fock::Fock(Command& cmd, IntegralEngine& ints, SharedMolecule m) : integrals(ints), molecule(m)
{
	
	// Make the core hamiltonian matrix
	formHCore();

	// Form the orthogonalising matrix and get initial density guess
	try {
		formOrthog();
	} catch (Error e) {
		molecule->control->log.error(e);
	}

	// Retrieve whether this calculation should be done direct, or whether
	// the twoints matrix has been formed, or if the 2e integrals need to be
	// read from file.
	direct = molecule->control->get_option<bool>("direct");
	diis = cmd.get_option<bool>("diis");
	iter = 0;
	nocc = molecule->getNel() / 2; 
	MAX = cmd.get_option<int>("maxdiis");
	precision = cmd.get_option<double>("precision");
	
	twoints = false;
	if (!direct){
		Vector ests = integrals.getEstimates();
		if (ests[3] < molecule->control->get_option<double>("memory")) 
			twoints = true;
	}

	fromfile = false;
	if (!twoints && !direct)
		fromfile = true;

}

Fock::Fock(const Fock& other) : integrals(other.integrals), molecule(other.molecule) {
	hcore = other.hcore;
	jkints = other.jkints;
	jints = other.jints;
	kints = other.kints;
	orthog = other.orthog;
	fockm = other.fockm;
	focka = other.focka;
	CP = other.CP;
	forces = other.forces;
	hessian = other.hessian;
	eps = other.eps;
	focks = other.focks;
	dens = other.dens;
	direct = other.direct;
	twoints = other.twoints;
	fromfile = other.fromfile;
	diis = other.diis;
	nbfs = other.nbfs;
	nocc = other.nocc;
	iter = other.iter;
	MAX = other.MAX; 
	precision = other.precision; 
}

// Form the core hamiltonian matrix
void Fock::formHCore()
{
	// The core Hamiltonian matrix is defined to be
	// the sum of the kinetic and nuclear attraction
	// matrices
	hcore = integrals.getKinetic() + integrals.getNucAttract();
	nbfs = hcore.rows();

}

void Fock::formOrthog()
{

	Matrix& S = integrals.getOverlap();
	
	// Diagonalise the overlap matrix into lambda and U,
	// so that U(T)SU = lambda
	// We can now form S^(-1/2) - the orthogonalising matrix
	EigenSolver es(S);
	Matrix U = es.eigenvectors();
	Vector lambda = es.eigenvalues(); 
  
	orthog = Matrix::Zero(nbfs, nbfs);
	for (int i = 0; i < nbfs; i++)
		orthog(i, i) = 1.0/(std::sqrt(lambda(i)));
  
	// S^-1/2  = U(lambda^-1/2)U(T)
	orthog = U * orthog * U.transpose();
}
  
void Fock::average(Vector &w) {
	if (diis && iter > 2) {
		// Average the fock matrices according to the weights
		focka = Matrix::Zero(nbfs, nbfs);
		int offset = focks.size() - w.size();
		for (int i = offset; i < focks.size(); i++) {
			focka = focka + w[i-offset]*focks[i]; 
		} 
	}
}

// Transform the AO fock matrix to the MO basis 
void Fock::transform(bool first)
{
	if (first) { 
		// Form the core Fock matrix as (S^-1/2)(T)H(S^-1/2)
		fockm = (orthog.transpose()) * ( hcore * orthog);
	} else {
		// Form the orthogonalised fock matrix
		fockm = orthog.transpose() * (focka * orthog);
	}
}

// Diagonalise the MO fock matrix to get CP and eps
void Fock::diagonalise() 
{
	EigenSolver es(fockm);
	CP = orthog * es.eigenvectors();
	eps = es.eigenvalues();
}

// Construct the density matrix from CP, 
// for nocc number of occupied orbitals
void Fock::makeDens()
{
	// Form the density matrix
	dens = 2.0 * CP.block(0, 0, nbfs, nocc) * CP.transpose().block(0, 0, nocc, nbfs); 
}

// Make the JK matrix, depending on how two electron integrals are stored/needed
void Fock::makeJK()
{
	if (twoints){
		formJK(dens); 
	} else if (direct) {
		formJKdirect(integrals.getPrescreen(), dens);
	} else {
		try {
			formJKfile();
		} catch (Error e) {
			molecule->control->log.error(e);
		}
	}
}

// Form the 2J-K matrix, given that twoints is stored in memory
void Fock::formJK(Matrix& P, double multiplier)
{
	jints = Matrix::Zero(nbfs, nbfs);
	kints = Matrix::Zero(nbfs, nbfs);
	for (int u = 0; u < nbfs; u++){
		for (int v = 0; v < nbfs ; v++){
			for (int s = 0; s < nbfs; s++){
				for (int l = 0; l < nbfs; l++){
					jints(u, v) += multiplier * P(s, l)*integrals.getERI(u, v, s, l);
					kints(u, v) += multiplier * P(s, l)*integrals.getERI(u, l, s, v);
				}
			}
		}
	}
	jkints = jints - 0.5*kints;
}

// Form JK using integral direct methods
void Fock::formJKdirect(const Matrix& Schwarz, Matrix& D, double multiplier)
{
	using libint2::Shell;
	using libint2::Engine;
	using libint2::Operator;
	
	std::vector<Shell>& shells = molecule->getBasis().getIntShells();
	const auto n = integrals.nbasis(shells);
	jints = Matrix::Zero(n, n);
	kints = Matrix::Zero(n, n);
	
	const auto do_schwarz_screen = Schwarz.cols() != 0 && Schwarz.rows() != 0;
	EMatrix D_shblk_norm = integrals.compute_shellblock_norm(shells, D);
	auto fock_precision = precision;
	auto maxnprim =  integrals.max_nprim(shells);
	auto maxnprim4 = maxnprim * maxnprim * maxnprim * maxnprim;
	auto engine_precision = std::min(fock_precision / D_shblk_norm.maxCoeff(), std::numeric_limits<double>::epsilon()) / maxnprim4;
	
	Engine engine(Operator::coulomb, integrals.max_nprim(shells), integrals.max_l(shells), 0);
	engine.set_precision(engine_precision);
	auto shell2bf = integrals.map_shell_to_basis_function(shells); 
	
	const auto& buf = engine.results();
	
	for(auto s1=0; s1!=shells.size(); ++s1) {

		auto bf1_first = shell2bf[s1]; // first basis function in this shell
		auto n1 = shells[s1].size();   // number of basis functions in this shell

		for(auto s2=0; s2<=s1; ++s2) {

			auto bf2_first = shell2bf[s2];
			auto n2 = shells[s2].size();
		
			const auto Dnorm12 = do_schwarz_screen ? D_shblk_norm(s1, s2) : 0.0;

			for(auto s3=0; s3<=s1; ++s3) {

				auto bf3_first = shell2bf[s3];
				auto n3 = shells[s3].size();
		  
				const auto Dnorm123 = do_schwarz_screen ? std::max(D_shblk_norm(s1, s3),
				std::max(D_shblk_norm(s2, s3), Dnorm12)) : 0.0;

				const auto s4_max = (s1 == s3) ? s2 : s3;
				for(auto s4=0; s4<=s4_max; ++s4) {

					const auto Dnorm1234 = do_schwarz_screen ? std::max(D_shblk_norm(s1, s4),
					std::max(D_shblk_norm(s2, s4), std::max(D_shblk_norm(s3, s4), Dnorm123))) : 0.0;
						
					if(do_schwarz_screen && Dnorm1234 * Schwarz(s1, s2) * Schwarz(s3, s4) < fock_precision) 
						continue;
					
					auto bf4_first = shell2bf[s4];
					auto n4 = shells[s4].size();

					// compute the permutational degeneracy (i.e. # of equivalents) of the given shell set
					auto s12_deg = (s1 == s2) ? 1.0 : 2.0;
					auto s34_deg = (s3 == s4) ? 1.0 : 2.0;
					auto s12_34_deg = (s1 == s3) ? (s2 == s4 ? 1.0 : 2.0) : 2.0;
					auto s1234_deg = s12_deg * s34_deg * s12_34_deg;

					engine.compute(shells[s1], shells[s2], shells[s3], shells[s4]);
					const auto* buf_1234 = buf[0];
					if (buf_1234 == nullptr)
						continue; // if all integrals screened out, skip to next quartet

					for(auto f1=0, f1234=0; f1!=n1; ++f1) {
						const auto bf1 = f1 + bf1_first;
						for(auto f2=0; f2!=n2; ++f2) {
							const auto bf2 = f2 + bf2_first;
							for(auto f3=0; f3!=n3; ++f3) {
								const auto bf3 = f3 + bf3_first;
								for(auto f4=0; f4!=n4; ++f4, ++f1234) {
									const auto bf4 = f4 + bf4_first;

									const auto value = buf_1234[f1234];

									const auto value_scal_by_deg = value * s1234_deg;

									jints(bf1,bf2) += multiplier * D(bf3,bf4) * value_scal_by_deg;
									jints(bf3,bf4) += multiplier * D(bf1,bf2) * value_scal_by_deg;
									kints(bf1,bf3) += 0.5 * multiplier * D(bf2,bf4) * value_scal_by_deg;
									kints(bf2,bf4) += 0.5 * multiplier * D(bf1,bf3) * value_scal_by_deg;
									kints(bf1,bf4) += 0.5 * multiplier * D(bf2,bf3) * value_scal_by_deg;
									kints(bf2,bf3) += 0.5 * multiplier * D(bf1,bf4) * value_scal_by_deg;
								}
							}
						}
					}

				}
			}
		}
	}
	
	// symmetrize the result
	jints = 0.25 * (jints + jints.transpose()).eval();
	kints = 0.25 * (kints + kints.transpose()).eval();
	jkints = jints - 0.5*kints;
}

// Form the JK matrix from two electron integrals stored on file
void Fock::formJKfile()
{
}
		

void Fock::makeFock()
{
	focka = hcore + jkints;
	if (diis) { // Archive for averaging
		if (iter >= MAX) {
			focks.erase(focks.begin());
		}
		focks.push_back(focka);
		iter++;
	}		
}

void Fock::simpleAverage(Matrix& D0, double weight)
{
	dens = weight*dens + (1.0-weight)*D0;
}


void Fock::clearDiis() {
	iter = 0;
	focks.clear(); 
}

FockFragment::FockFragment(Command& cmd, IntegralEngine& ints, SharedMolecule m, int _start, int _end) : Fock(cmd, ints, m), start(_start), end(_end) { 
	Sxx = ints.getOverlap();
}

FockFragment::FockFragment(const FockFragment& other) : Fock(other) {
	Sxx = other.Sxx;
	start = other.start;
	end = other.end;
}

Vector FockFragment::buildFock(Matrix& qfq, Matrix& qfp, Matrix& pfp, bool alpha) 
{
	int nbfs = end - start; 
	
	focka = qfp.block(start, start, nbfs, nbfs) * Sxx; 
	Matrix e = focka - focka.transpose(); 
	Vector err(Eigen::Map<Vector>(e.data(), e.cols()*e.rows()));
	
	focka = focka + focka.transpose() + qfq.block(start, start, nbfs, nbfs) + Sxx * pfp.block(start, start, nbfs, nbfs) * Sxx; 
	
	if (diis) { // Archive for averaging
		if (iter >= MAX) {
			focks.erase(focks.begin());
		}
		focks.push_back(focka);
	}		
	
	return err; 
}

void FockFragment::gensolve()
{
	GeneralizedEigenSolver es(focka, Sxx); 
	CP = es.eigenvectors();
	eps = es.eigenvalues();
	iter++;
}

UnrestrictedFock::UnrestrictedFock(Command& cmd, IntegralEngine& ints, SharedMolecule m) : Fock(cmd, ints, m)
{
	nalpha = molecule->nalpha();
	nbeta = molecule->nbeta(); 
}

UnrestrictedFock::UnrestrictedFock(const UnrestrictedFock& other) : Fock(other)
{
	fock_alpha_ao = other.fock_alpha_ao;
	fock_beta_ao = other.fock_alpha_ao;
	fock_alpha_mo = other.fock_alpha_ao;
	fock_beta_mo = other.fock_alpha_ao;
	dens_alpha = other.dens_alpha;
	dens_beta = other.dens_beta;
	CP_alpha = other.CP_alpha;
	eps_alpha = other.eps_alpha;
	CP_beta = other.CP_beta;
	eps_beta = other.eps_beta;
	jkints_alpha = other.jkints_alpha;
	jkints_beta = other.jkints_beta;
	nalpha = other.nalpha;
	nbeta = other.nbeta; 
}

// Form the 2J-K matrix, given that twoints is stored in memory
void UnrestrictedFock::formJK(Matrix& Pa, Matrix& Pb, double multiplier)
{
	jkints_alpha = Matrix::Zero(nbfs, nbfs);
	jkints_beta = Matrix::Zero(nbfs, nbfs);
	for (int u = 0; u < nbfs; u++){
		for (int v = 0; v < nbfs ; v++){
			for (int s = 0; s < nbfs; s++){
				for (int l = 0; l < nbfs; l++){
					auto ival = multiplier * integrals.getERI(u, v, s, l); 
					jkints_alpha(u, v) += Pa(s, l) * ival;
					jkints_beta(u, v) += Pb(s, l) * ival;
					ival = 0.5 * multiplier * integrals.getERI(u, l, s, v);
					jkints_alpha(u, v) -= Pa(s, l) * ival;
					jkints_beta(u, v) -= Pb(s, l) * ival;
				}
			}
		}
	}
}

// Form JK using integral direct methods
void UnrestrictedFock::formJKdirect(const Matrix& Schwarz, Matrix& Da, Matrix& Db, double multiplier)
{
	using libint2::Shell;
	using libint2::Engine;
	using libint2::Operator;
	
	std::vector<Shell>& shells = molecule->getBasis().getIntShells();
	const auto n = integrals.nbasis(shells);
	jkints_alpha = Matrix::Zero(n, n);
	jkints_beta = Matrix::Zero(n, n);
	
	const auto do_schwarz_screen = Schwarz.cols() != 0 && Schwarz.rows() != 0;
	Matrix D_tot = Da + Db; 
	EMatrix D_shblk_norm = integrals.compute_shellblock_norm(shells, D_tot);
	auto fock_precision = precision;
	auto maxnprim =  integrals.max_nprim(shells);
	auto maxnprim4 = maxnprim * maxnprim * maxnprim * maxnprim;
	auto engine_precision = std::min(fock_precision / D_shblk_norm.maxCoeff(), std::numeric_limits<double>::epsilon()) / maxnprim4;
	
	Engine engine(Operator::coulomb, integrals.max_nprim(shells), integrals.max_l(shells), 0);
	engine.set_precision(engine_precision);
	auto shell2bf = integrals.map_shell_to_basis_function(shells); 
	
	const auto& buf = engine.results();
	
	for(auto s1=0; s1!=shells.size(); ++s1) {

		auto bf1_first = shell2bf[s1]; // first basis function in this shell
		auto n1 = shells[s1].size();   // number of basis functions in this shell

		for(auto s2=0; s2<=s1; ++s2) {

			auto bf2_first = shell2bf[s2];
			auto n2 = shells[s2].size();
		
			const auto Dnorm12 = do_schwarz_screen ? D_shblk_norm(s1, s2) : 0.0;

			for(auto s3=0; s3<=s1; ++s3) {

				auto bf3_first = shell2bf[s3];
				auto n3 = shells[s3].size();
		  
				const auto Dnorm123 = do_schwarz_screen ? std::max(D_shblk_norm(s1, s3),
				std::max(D_shblk_norm(s2, s3), Dnorm12)) : 0.0;

				const auto s4_max = (s1 == s3) ? s2 : s3;
				for(auto s4=0; s4<=s4_max; ++s4) {

					const auto Dnorm1234 = do_schwarz_screen ? std::max(D_shblk_norm(s1, s4),
					std::max(D_shblk_norm(s2, s4), std::max(D_shblk_norm(s3, s4), Dnorm123))) : 0.0;
						
					if(do_schwarz_screen && Dnorm1234 * Schwarz(s1, s2) * Schwarz(s3, s4) < fock_precision) 
						continue;
					
					auto bf4_first = shell2bf[s4];
					auto n4 = shells[s4].size();

					// compute the permutational degeneracy (i.e. # of equivalents) of the given shell set
					auto s12_deg = (s1 == s2) ? 1.0 : 2.0;
					auto s34_deg = (s3 == s4) ? 1.0 : 2.0;
					auto s12_34_deg = (s1 == s3) ? (s2 == s4 ? 1.0 : 2.0) : 2.0;
					auto s1234_deg = s12_deg * s34_deg * s12_34_deg;

					engine.compute(shells[s1], shells[s2], shells[s3], shells[s4]);
					const auto* buf_1234 = buf[0];
					if (buf_1234 == nullptr)
						continue; // if all integrals screened out, skip to next quartet

					for(auto f1=0, f1234=0; f1!=n1; ++f1) {
						const auto bf1 = f1 + bf1_first;
						for(auto f2=0; f2!=n2; ++f2) {
							const auto bf2 = f2 + bf2_first;
							for(auto f3=0; f3!=n3; ++f3) {
								const auto bf3 = f3 + bf3_first;
								for(auto f4=0; f4!=n4; ++f4, ++f1234) {
									const auto bf4 = f4 + bf4_first;

									const auto value = buf_1234[f1234];

									const auto value_scal_by_deg = multiplier * value * s1234_deg;

									jkints_alpha(bf1,bf2) += Da(bf3,bf4) * value_scal_by_deg;
									jkints_alpha(bf3,bf4) += Da(bf1,bf2) * value_scal_by_deg;
									jkints_alpha(bf1,bf3) -= 0.25 * Da(bf2,bf4) * value_scal_by_deg;
									jkints_alpha(bf2,bf4) -= 0.25 * Da(bf1,bf3) * value_scal_by_deg;
									jkints_alpha(bf1,bf4) -= 0.25 * Da(bf2,bf3) * value_scal_by_deg;
									jkints_alpha(bf2,bf3) -= 0.25 * Da(bf1,bf4) * value_scal_by_deg;
									
									jkints_beta(bf1,bf2) += Db(bf3,bf4) * value_scal_by_deg;
									jkints_beta(bf3,bf4) += Db(bf1,bf2) * value_scal_by_deg;
									jkints_beta(bf1,bf3) -= 0.25 * Db(bf2,bf4) * value_scal_by_deg;
									jkints_beta(bf2,bf4) -= 0.25 * Db(bf1,bf3) * value_scal_by_deg;
									jkints_beta(bf1,bf4) -= 0.25 * Db(bf2,bf3) * value_scal_by_deg;
									jkints_beta(bf2,bf3) -= 0.25 * Db(bf1,bf4) * value_scal_by_deg;
								}
							}
						}
					}

				}
			}
		}
	}
	
	// symmetrize the result
	jkints_alpha = 0.25 * (jkints_alpha + jkints_alpha.transpose()).eval();
	jkints_beta = 0.25 * (jkints_beta + jkints_beta.transpose()).eval();
}

// Transform the AO fock matrix to the MO basis 
void UnrestrictedFock::transform(bool first)
{
	if (first) { 
		// Form the core Fock matrix as (S^-1/2)(T)H(S^-1/2)
		fock_alpha_mo = (orthog.transpose()) * ( hcore * orthog);
		fock_beta_mo = fock_alpha_mo; 
	} else {
		// Form the orthogonalised fock matrix
		fock_alpha_mo = orthog.transpose() * (fock_alpha_ao * orthog);
		fock_beta_mo = orthog.transpose() * (fock_beta_ao * orthog);
	}
}

// Diagonalise the MO fock matrix to get CP and eps
void UnrestrictedFock::diagonalise() 
{
	EigenSolver es_alpha(fock_alpha_mo);
	EigenSolver es_beta(fock_beta_mo);
	CP_alpha = orthog * es_alpha.eigenvectors();
	CP_beta = orthog * es_beta.eigenvectors(); 
	eps_alpha = es_alpha.eigenvalues();
	eps_beta = es_beta.eigenvalues();
}

// Construct the density matrix from CP, 
// for nocc number of occupied orbitals
void UnrestrictedFock::makeDens()
{
	// Form the density matrix
	dens_alpha = 2.0 * CP_alpha.block(0, 0, nbfs, nalpha) * CP_alpha.transpose().block(0, 0, nalpha, nbfs);
	dens_beta = 2.0 * CP_beta.block(0, 0, nbfs, nbeta) * CP_beta.transpose().block(0, 0, nbeta, nbfs); 
}

// Make the JK matrix, depending on how two electron integrals are stored/needed
void UnrestrictedFock::makeJK()
{
	if (twoints){
		formJK(dens_alpha, dens_beta); 
	} else if (direct) {
		formJKdirect(integrals.getPrescreen(), dens_alpha, dens_beta);
	} else {
		try {
			formJKfile();
		} catch (Error e) {
			molecule->control->log.error(e);
		}
	}
}

void UnrestrictedFock::makeFock()
{
	fock_alpha_ao = hcore + jkints_alpha;
	fock_beta_ao = hcore + jkints_beta; 
	if (diis) { // Archive for averaging
		if (iter >= MAX) {
			alpha_focks.erase(alpha_focks.begin());
			beta_focks.erase(beta_focks.begin());
		}
		alpha_focks.push_back(fock_alpha_ao);
		beta_focks.push_back(fock_beta_ao);
		iter++;
	}		
}

void UnrestrictedFock::average(Vector &w) {
	if (diis && iter > 2) {
		// Average the fock matrices according to the weights
		fock_alpha_ao = Matrix::Zero(nbfs, nbfs);
		fock_beta_ao = Matrix::Zero(nbfs, nbfs);
		int offset = alpha_focks.size() - w.size();
		for (int i = offset; i < alpha_focks.size(); i++) {
			fock_alpha_ao = fock_alpha_ao + w[i-offset]*alpha_focks[i];
			fock_beta_ao = fock_beta_ao + w[i-offset]*beta_focks[i]; 
		} 
	}
}

void UnrestrictedFock::formJKfile() {
	std::cerr << "File-read not implemented yet!" << std::endl;  
}

void UnrestrictedFock::clearDiis() {
	iter = 0;
	alpha_focks.clear();
	beta_focks.clear();
}

UnrestrictedFockFragment::UnrestrictedFockFragment(Command& cmd, IntegralEngine& ints, SharedMolecule m, int _start, int _end) : UnrestrictedFock(cmd, ints, m), start(_start), end(_end) { 
	Sxx = ints.getOverlap();
}

UnrestrictedFockFragment::UnrestrictedFockFragment(const UnrestrictedFockFragment& other) : UnrestrictedFock(other) {
	Sxx = other.Sxx;
	start = other.start;
	end = other.end;
}

Vector UnrestrictedFockFragment::buildFock(Matrix& qfq, Matrix& qfp, Matrix& pfp, bool alpha) 
{
	int nbfs = end - start; 
	
	focka = qfp.block(start, start, nbfs, nbfs) * Sxx; 
	
	Matrix e = focka - focka.transpose(); 
	Vector err(Eigen::Map<Vector>(e.data(), e.cols()*e.rows()));
	
	focka = focka + focka.transpose() + qfq.block(start, start, nbfs, nbfs) + Sxx * pfp.block(start, start, nbfs, nbfs) * Sxx; 
	
	if (alpha) {
		fock_alpha_ao = focka;
		if (diis) { // Archive for averaging
			if (iter >= MAX)
				alpha_focks.erase(alpha_focks.begin());
			alpha_focks.push_back(fock_alpha_ao);
		}	
	} else { 
		fock_beta_ao = focka;
		if (diis) { // Archive for averaging
			if (iter >= MAX)
				beta_focks.erase(beta_focks.begin());
			beta_focks.push_back(fock_beta_ao);
		}	
	}
	return err; 
}

void UnrestrictedFockFragment::gensolve()
{
	GeneralizedEigenSolver es_alpha(fock_alpha_ao, Sxx); 
	CP_alpha = es_alpha.eigenvectors();
	eps_alpha = es_alpha.eigenvalues(); 
	
	GeneralizedEigenSolver es_beta(fock_beta_ao, Sxx);
	CP_beta = es_beta.eigenvectors();
	eps_beta = es_beta.eigenvalues(); 
	
	iter++;
}

