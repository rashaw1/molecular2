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

typedef std::vector<libint2::Shell> BasisSet; 
 
void Fock::makeJK(Matrix& P, double multiplier)
{
	
	if (twoints){
		formJK(P, multiplier); 
	} else if (density_fitted) {
		formJKdf(P, multiplier); 
	} else if (direct) {
		formJKdirect(integrals.getPrescreen(), P, multiplier);
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
	
void Fock::formJKdf(Matrix& Cocc, double multiplier) {
	jkints = multiplier * compute_2body_fock_df(Cocc); 
}	

void Fock::makeFock(Matrix& P, double multiplier)
{
	if (not density_fitted) {
		focka = fockinc; 
		if (not started_incremental &&
		rms_error < incremental_threshold) {
			started_incremental = true;
			reset_incremental = false;
			next_reset = rms_error / 1e1;
			last_reset = iter - 1; 
		}
		if (reset_incremental || not started_incremental) {
			focka = hcore;
			dens_diff = dens;
		}
		if (reset_incremental && started_incremental) {
			reset_incremental = false;
			last_reset = iter; 
			next_reset = rms_error / 1e1;
		}
	} else 
		focka = hcore; 
	
	makeJK(P, multiplier);  
	focka += jkints; 
	if (not density_fitted) fockinc = focka; 
	addDiis(); 
}

void Fock::addDiis() {
	if (diis) { // Archive for averaging
		if (iter >= MAX) {
			focks.erase(focks.begin());
		}
		focks.push_back(focka);
		iter++;
	}		
}

void Fock::makeFock() {
	if (not density_fitted)
		makeFock(dens_diff); 
	else {
		Matrix cocc = CP.block(0, 0, nbfs, nocc); 
		makeFock(cocc); 
	}
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
// Form the 2J-K matrix, given that twoints is stored in memory
void UnrestrictedFock::formJK(Matrix& Pa, Matrix& Pb, double multiplier)
{
	kints_alpha = Matrix::Zero(nbfs, nbfs);
	kints_beta = Matrix::Zero(nbfs, nbfs); 
	jints_alpha = Matrix::Zero(nbfs, nbfs);
	jints_beta = Matrix::Zero(nbfs, nbfs);
	for (int u = 0; u < nbfs; u++){
		for (int v = 0; v < nbfs ; v++){
			for (int s = 0; s < nbfs; s++){
				for (int l = 0; l < nbfs; l++){
					auto ival = multiplier * integrals.getERI(u, v, s, l); 
					jints_alpha(u, v) += Pa(s, l) * ival;
					jints_beta(u, v) += Pb(s, l) * ival;
					ival =  multiplier * integrals.getERI(u, l, s, v);
					kints_alpha(u, v) += Pa(s, l) * ival;
					kints_beta(u, v) += Pb(s, l) * ival;
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
	kints_alpha = Matrix::Zero(n, n);
	kints_beta = Matrix::Zero(n, n); 
	jints_alpha = Matrix::Zero(n, n);
	jints_beta = Matrix::Zero(n, n);
		
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
									jints_alpha(bf1,bf2) += Da(bf3, bf4) * value_scal_by_deg;
									jints_alpha(bf3,bf4) += Da(bf1, bf2) * value_scal_by_deg;
									kints_alpha(bf1,bf3) += 0.5 * Da(bf2,bf4) * value_scal_by_deg;
									kints_alpha(bf2,bf4) += 0.5 * Da(bf1,bf3) * value_scal_by_deg;
									kints_alpha(bf1,bf4) += 0.5 * Da(bf2,bf3) * value_scal_by_deg;
									kints_alpha(bf2,bf3) += 0.5 * Da(bf1,bf4) * value_scal_by_deg;
										
									jints_beta(bf1,bf2) += Db(bf3, bf4) * value_scal_by_deg;
									jints_beta(bf3,bf4) += Db(bf1, bf2) * value_scal_by_deg;
									kints_beta(bf1,bf3) += 0.5 * Db(bf2,bf4) * value_scal_by_deg;
									kints_beta(bf2,bf4) += 0.5 * Db(bf1,bf3) * value_scal_by_deg;
									kints_beta(bf1,bf4) += 0.5 * Db(bf2,bf3) * value_scal_by_deg;
									kints_beta(bf2,bf3) += 0.5 * Db(bf1,bf4) * value_scal_by_deg;
								}
							}
						}
					}
				}
			}
		}
	}
	
	// symmetrize the result
	kints_alpha = 0.25 * (kints_alpha + kints_alpha.transpose()).eval();
	kints_beta = 0.25 * (kints_beta + kints_beta.transpose()).eval();
	jints_alpha = 0.25 * (jints_alpha + jints_alpha.transpose()).eval();
	jints_beta = 0.25 * (jints_beta + jints_beta.transpose()).eval();
}

void UnrestrictedFock::formJKdf(Matrix& ca_occ, Matrix& cb_occ, double multiplier) {
	/*jkints_alpha = multiplier * compute_2body_fock_df(ca_occ); 
	jkints_beta = multiplier * compute_2body_fock_df(cb_occ); */
}

// Make the JK matrix, depending on how two electron integrals are stored/needed
void UnrestrictedFock::makeJK()
{
	if (twoints){
		formJK(dens_alpha, dens_beta); 
	} else if (density_fitted) {
		Matrix ca_occ = CP_alpha.block(0, 0, nbfs, nalpha); 
		Matrix cb_occ = CP_beta.block(0, 0, nbfs, nbeta); 
		formJKdf(ca_occ, cb_occ); 
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
	fock_alpha_ao = hcore + jints_alpha + jints_beta;
	fock_beta_ao = fock_alpha_ao - kints_beta; 
	fock_alpha_ao -= kints_alpha; 
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

// an Fock builder that can accept densities expressed a separate basis
Matrix Fock::compute_2body_fock_general(const BasisSet& obs, const Matrix& D, const BasisSet& D_bs,
bool D_is_shelldiagonal, double precision)
{
		
	const auto n = integrals.nbasis(obs);
	const auto nshells = obs.size();
	const auto n_D = integrals.nbasis(D_bs);
	assert(D.cols() == D.rows() && D.cols() == n_D);

	Matrix G = Matrix::Zero(n, n);

	// construct the 2-electron repulsion integrals engine
	using libint2::Engine;
	Engine engine(libint2::Operator::coulomb,
	std::max(integrals.max_nprim(obs), integrals.max_nprim(D_bs)),
	std::max(integrals.max_l(obs), integrals.max_l(D_bs)), 0);
	engine.set_precision(precision);  // shellset-dependent precision control
	// will likely break positive
	// definiteness
	// stick with this simple recipe

	auto shell2bf = integrals.map_shell_to_basis_function(obs);
	auto shell2bf_D = integrals.map_shell_to_basis_function(D_bs);

	const auto& buf = engine.results();

	// loop over permutationally-unique set of shells
	for (auto s1 = 0l, s1234 = 0l; s1 != nshells; ++s1) {
		auto bf1_first = shell2bf[s1];  // first basis function in this shell
		auto n1 = obs[s1].size();       // number of basis functions in this shell

		for (auto s2 = 0; s2 <= s1; ++s2) {
			auto bf2_first = shell2bf[s2];
			auto n2 = obs[s2].size();

			for (auto s3 = 0; s3 < D_bs.size(); ++s3) {
				auto bf3_first = shell2bf_D[s3];
				auto n3 = D_bs[s3].size();

				auto s4_begin = D_is_shelldiagonal ? s3 : 0;
				auto s4_fence = D_is_shelldiagonal ? s3 + 1 : D_bs.size();

				for (auto s4 = s4_begin; s4 != s4_fence; ++s4, ++s1234) {
					auto bf4_first = shell2bf_D[s4];
					auto n4 = D_bs[s4].size();

					// compute the permutational degeneracy (i.e. # of equivalents) of
					// the given shell set
					auto s12_deg = (s1 == s2) ? 1.0 : 2.0;

					if (s3 >= s4) {
						auto s34_deg = (s3 == s4) ? 1.0 : 2.0;
						auto s1234_deg = s12_deg * s34_deg;
						// auto s1234_deg = s12_deg;
						engine.compute2<libint2::Operator::coulomb, libint2::BraKet::xx_xx, 0>(
							obs[s1], obs[s2], D_bs[s3], D_bs[s4]);
						const auto* buf_1234 = buf[0];
						if (buf_1234 != nullptr) {
							for (auto f1 = 0, f1234 = 0; f1 != n1; ++f1) {
								const auto bf1 = f1 + bf1_first;
								for (auto f2 = 0; f2 != n2; ++f2) {
									const auto bf2 = f2 + bf2_first;
									for (auto f3 = 0; f3 != n3; ++f3) {
										const auto bf3 = f3 + bf3_first;
										for (auto f4 = 0; f4 != n4; ++f4, ++f1234) {
											const auto bf4 = f4 + bf4_first;

											const auto value = buf_1234[f1234];
											const auto value_scal_by_deg = value * s1234_deg;
											G(bf1, bf2) += 2.0 * D(bf3, bf4) * value_scal_by_deg;
										}
									}
								}
							}
						}
					}

					engine.compute2<libint2::Operator::coulomb, libint2::BraKet::xx_xx, 0>(
						obs[s1], D_bs[s3], obs[s2], D_bs[s4]);
					const auto* buf_1324 = buf[0];
					if (buf_1324 == nullptr)
						continue; // if all integrals screened out, skip to next quartet

					for (auto f1 = 0, f1324 = 0; f1 != n1; ++f1) {
						const auto bf1 = f1 + bf1_first;
						for (auto f3 = 0; f3 != n3; ++f3) {
							const auto bf3 = f3 + bf3_first;
							for (auto f2 = 0; f2 != n2; ++f2) {
								const auto bf2 = f2 + bf2_first;
								for (auto f4 = 0; f4 != n4; ++f4, ++f1324) {
									const auto bf4 = f4 + bf4_first;

									const auto value = buf_1324[f1324];
									const auto value_scal_by_deg = value * s12_deg;
									G(bf1, bf2) -= D(bf3, bf4) * value_scal_by_deg;
								}
							}
						}
					}
				}
			}
		}
	}

	// symmetrize the result and return
	return 0.5 * (G + G.transpose());		
}

Matrix Fock::compute_2body_fock_df(const Matrix& Cocc) {

	BasisSet& obs = molecule->getBasis().getIntShells(); 
	BasisSet& dfbs = molecule->getBasis().getDFShells(); 
  
	const auto n = integrals.nbasis(obs); 
	const auto ndf = integrals.nbasis(dfbs); 
	
	Logger& log = molecule->control->log; 

	// using first time? compute 3-center ints and transform to inv sqrt
	// representation
	if (xyK.rows() == 0) {
		
		integrals.compute_eris_3index(obs, dfbs, xyK); 
		Matrix V = integrals.compute_eris_2index(dfbs);
	
		Eigen::LLT<Matrix> V_LLt(V);
		Matrix I = Matrix::Identity(ndf, ndf);
		auto L = V_LLt.matrixL();
		Linv = L.solve(I).transpose();

		Vector row; 
		for (int x = 0; x < n; x++) {
			for (int y = 0; y<= x; y++) {
				int xy = x*n+y; 
				int yx = y*n+x; 
				row = xyK.row(xy); 
			
				for (int K = 0; K < ndf; K++)
					xyK(xy, K) = xyK(yx, K) = row.dot(Linv.col(K)); 
			}
		} 
	}  // if (xyK.size() == 0)

	// compute exchange
	double start = log.getGlobalTime(); 
	std::cout << "Exchange: "; 
	Matrix xiK = Matrix::Zero(n*nocc, ndf); 
	for (int x = 0; x < n; x++) {
		for (int i = 0; i < nocc; i++) {
			for (int K = 0; K < ndf; K++) {
				for (int y = 0; y < n; y++) 
					xiK(x*nocc+i, K) += xyK(x*n+y, K) * Cocc(y, i); 
			}
		}
	}

	Matrix G = Matrix::Zero(n, n); 
	for (int x = 0; x < n; x++) {
		for (int y = 0; y <= x; y++) {
			for (int i = 0; i < nocc; i++)
				for (int K = 0; K < ndf; K++)
					G(x, y) -= xiK(x*nocc+i, K) * xiK(y*nocc+i, K); 
		}
	}
	double end = log.getGlobalTime(); 
	std::cout << end - start << " seconds" << std::endl; 

	// compute Coulomb
	start = end;
	std::cout << "Coulomb: "; 
	Vector Jtmp = Vector::Zero(ndf); 
	for (int K = 0; K < ndf; K++) 
		for (int x = 0; x < n; x++)
			for (int i = 0; i < nocc; i++)
				Jtmp[K] += xiK(x*nocc+i, K) * Cocc(x, i); 
	xiK.resize(0, 0); 
  
	for (int x = 0; x < n; x++) {
		for (int y = 0; y <= x; y++) { 
			for (int K = 0; K < ndf; K++)  
				G(x, y) += 2.0 * xyK(x*n+y, K) * Jtmp[K]; 
			G(y, x) = G(x, y); 
		}
	}
	end = log.getGlobalTime(); 
	std::cout << end - start << " seconds" << std::endl << std::endl; 

	return G; 
}

Matrix Fock::compute_2body_fock_df_local(Matrix& Cocc, const Matrix& sigmainv, Matrix& Pt, std::vector<FragmentInfo>& finfo) {

	BasisSet& obs = molecule->getBasis().getIntShells(); 
	BasisSet& dfbs = molecule->getBasis().getDFShells(); 
  
	const auto n = integrals.nbasis(obs); 
	const auto ndf = integrals.nbasis(dfbs); 
	
	Logger& log = molecule->control->log; 

	EigenSolver es(sigmainv); 
	Cocc = Cocc * es.operatorSqrt(); 
	
	// using first time? compute 3-center ints and transform to inv sqrt
	// representation
	if (xyK.rows() == 0) {
		
		integrals.compute_eris_3index(obs, dfbs, xyK); 
		Matrix V = integrals.compute_eris_2index(dfbs);
	
		Eigen::LLT<Matrix> V_LLt(V);
		Matrix I = Matrix::Identity(ndf, ndf);
		auto L = V_LLt.matrixL();
		Linv = L.solve(I).transpose();

		/*Vector row; 
		for (int x = 0; x < n; x++) {
			for (int y = 0; y<= x; y++) {
				int xy = x*n+y; 
				int yx = y*n+x; 
				row = xyK.row(xy); 
			
				for (int K = 0; K < ndf; K++)
					xyK(xy, K) = xyK(yx, K) = row.dot(Linv.col(K)); 
			}
		} */
	
		int i_offset = 0; 
		for (int B = 0; B < finfo.size(); B++) {
			for (int i = i_offset; i < i_offset+finfo[B].occ; i++) {
				Domain d; 
				
				std::vector<int> notcentres; 
				int mu_offset = 0;
				for (int A = 0; A < finfo.size(); A++) {
					double sum = 0.0; 
					for (int mu = mu_offset; mu < mu_offset + finfo[A].occ; mu++) 
						sum += Cocc(mu, i) * Cocc(mu, i); 
				
					if (sum > 1e-6) { 
						d.starts.push_back(finfo[A].start); 
						d.sizes.push_back(finfo[A].nbfs); 
						d.centres.push_back(A); 
					} else 
						notcentres.push_back(A); 
				
					mu_offset += finfo[A].nbfs; 
				}
				
				Domain dao; 
				dao.starts = d.starts;
				dao.sizes = d.sizes; 
				dao.centres = d.centres; 
				mu_offset = 0; 
				for (auto A : d.centres) {
					std::vector<int> newnc; 
					for (int b = 0; b < notcentres.size(); b++){
						int C = notcentres[b]; 
						double sep = (finfo[C].com - finfo[A].com).norm(); 
						if (finfo[C].radius > sep) {
							dao.centres.push_back(C); 
							dao.starts.push_back(finfo[C].start);
							dao.sizes.push_back(finfo[C].nbfs); 
						} else 
							newnc.push_back(C); 
					} 
					notcentres = newnc; 
				}
			
				d.sumsizes();
				dao.sumsizes(); 
				lmo_domains.push_back(d); 
				ao_domains.push_back(dao); 
			}
			i_offset += finfo[B].occ; 
		}
		
		EigenSolver es2(integrals.getOverlap()); 
		Matrix SC = es2.operatorSqrt() * Cocc; 
		i_offset = 0;
		for (auto& f1 : finfo) { 
			
			for (int i = i_offset; i < i_offset + f1.occ; i++) {
				Domain d;
				
				int mu_offset = 0; int center = 0; 
				for (auto& f2 : finfo) { 
					
					double sep = (f2.com - f1.com).norm(); 
					
					double INi = 0.0; 
					for (int mu = mu_offset; mu < mu_offset + f2.nbfs; mu++) {
						INi += SC(mu, i) * SC(mu, i);	
					}
					
					if (INi > 5e-2 || sep < 7.0) {
						d.centres.push_back(center);
						d.sizes.push_back(f2.naux); 
						d.starts.push_back(f2.auxstart); 
					}
				
					center++; 
					mu_offset += f2.nbfs; 
				}
				
				d.sumsizes();
				fit_domains.push_back(d); 
			}
			
			i_offset += f1.occ; 
		}
		
		for (int i = 0; i < nocc; i++) {
			
			auto& lmo_d = lmo_domains[i]; 
			auto& fit_d = fit_domains[i];
			int fitsize = fit_d.centres.size();  
			
			Matrix GAB = Matrix::Zero(fit_d.totalsize, fit_d.totalsize); 
			int A = 0; 
			for (int a = 0; a < fitsize; a++) {
				for (int Ax = fit_d.starts[a]; Ax < fit_d.starts[a] + fit_d.sizes[a]; Ax++) {
					
					int B = 0; 
					for (int b = 0; b < fitsize; b++) {
						for (int Bx = fit_d.starts[b]; Bx < fit_d.starts[b] + fit_d.sizes[b]; Bx++) {
							GAB(A, B) = V(Ax, Bx); 
							B++; 
						}
					}
					
					A++; 
				}
			}
			
			Eigen::LLT<Matrix> G_LLt(GAB);
			Matrix I = Matrix::Identity(fit_d.totalsize, fit_d.totalsize);
			auto GL = G_LLt.matrixL();
			lmo_d.G = GL.solve(I).transpose();
			
		}
		
	}  // if (xyK.size() == 0) 
	
	Matrix G = Matrix::Zero(n, n); 
	
	// compute exchange
	int nfrags = finfo.size(); 
	double start = log.getGlobalTime(); 

	for (int i = 0; i < nocc; i++) {
		auto& lmo_d = lmo_domains[i]; 
		auto& ao_d = ao_domains[i];
		auto& fit_d = fit_domains[i]; 
		
		int nmo = lmo_d.centres.size();
		int nao = ao_d.centres.size();
		int nfit = fit_d.centres.size(); 
		
		Matrix uA = Matrix::Zero(ao_d.totalsize, fit_d.totalsize); 
		
		int Al = 0; 
		int fitstart, aostart, mostart, fitsize, aosize, mosize; 
		for (int fit = 0; fit < nfit; fit++) {
			fitstart = fit_d.starts[fit]; 
			fitsize = fit_d.sizes[fit]; 
			
			for (int A = fitstart; A < fitstart + fitsize; A++) {
				
				int ul = 0;
				for (int ao = 0; ao < nao; ao++) {
					aostart = ao_d.starts[ao]; 
					aosize = ao_d.sizes[ao]; 
					
					for (int u = aostart; u < aostart + aosize; u++) {
						for (int mo = 0; mo < nmo; mo++) {
							mostart = lmo_d.starts[mo];
							mosize = lmo_d.sizes[mo]; 
							
							for (int nu = mostart; nu < mostart + mosize; nu++)  
								uA(ul, Al) += xyK(u*n+nu, A) * Cocc(nu, i); 
							
						}
						ul++; 
						
					}
				}
				Al++; 
			}
		}
		
		uA *= lmo_d.G; 
		uA = uA * uA.transpose(); 
		
		int ao1start, ao2start, ao1size, ao2size; 
		int u = 0; 
		for (int ao1 = 0; ao1 < nao; ao1++) {
			ao1start = ao_d.starts[ao1];
			ao1size = ao_d.sizes[ao1]; 
			
			for (int mu = ao1start; mu < ao1start + ao1size; mu++) {
				
				int v = 0; 
				for (int ao2 = 0; ao2 < nao; ao2++) {
					ao2start = ao_d.starts[ao2]; 
					ao2size = ao_d.sizes[ao2]; 
					 
					for (int nu = ao2start; nu < ao2start + ao2size; nu++) {

						G(mu, nu) -= uA(u, v); 
						
						v++; 
					}
				}
				u++; 
			}
		}
	}

	double end = log.getGlobalTime(); 
	std::cout << end - start << " seconds" << std::endl; 
	
	// compute Coulomb
	start = log.getGlobalTime();
	std::cout << "Coulomb: "; 
	
	Vector Pv(Eigen::Map<Vector>(Pt.data(), Pt.cols()*Pt.rows()));
	Pv = xyK.transpose() * Pv;
	Pv = Linv * Linv.transpose() * Pv;  
	Pv = 2.0 * xyK *  Pv; 
	for (int x = 0; x < n; x++)
		for (int y = 0; y <=x; y++) {
			G(x, y) += Pv[x*n+y];
			G(y, x) = G(x, y); 
		}
		
	end = log.getGlobalTime(); 
	std::cout << end - start << " seconds" << std::endl << std::endl; 

	return G; 
}