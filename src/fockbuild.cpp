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

void UnrestrictedFock::formJKdf(Matrix& ca_occ, Matrix& cb_occ, double multiplier) {
	jkints_alpha = multiplier * compute_2body_fock_df(ca_occ); 
	jkints_beta = multiplier * compute_2body_fock_df(cb_occ); 
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

  // using first time? compute 3-center ints and transform to inv sqrt
  // representation
  if (xyK.rows() == 0) {

    const auto nshells = obs.size();
    const auto nshells_df = dfbs.size();


    Matrix Zxy = integrals.compute_eris_3index(obs, dfbs); 
	Matrix V = integrals.compute_eris_2index(dfbs);
	
    Eigen::LLT<Matrix> V_LLt(V);
    Matrix I = Matrix::Identity(ndf, ndf);
    auto L = V_LLt.matrixL();
    Matrix V_L = L;
    Matrix Linv = L.solve(I).transpose();

    xyK = Zxy * Linv; 
	Zxy.resize(0, 0); 
	
  }  // if (xyK.size() == 0)

  // compute exchange

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

  // compute Coulomb

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

  return G; 
}
