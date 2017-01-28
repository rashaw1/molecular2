
#include "cc.hpp"
#include "logger.hpp"
#include "mvector.hpp"
#include "fock.hpp"
#include "molecule.hpp"
#include "integrals.hpp"
#include "error.hpp"
#include <iostream>
#include <thread>

// Constructor
CCSD::CCSD(MP2& _mp2, bool _triples) : mp2(_mp2), withTriples(_triples) {
	N = 2*mp2.getN();
	nocc = 2*mp2.getNocc();
	energy = 0.0;
	delta_e = 0.0;
	delta_singles = 0.0;
	delta_doubles = 0.0;
}

void CCSD::build_fock() {
	// Transform hcore to MO basis
	Matrix hcore_spatial = mp2.getFock().getCP().transpose() * mp2.getFock().getHCore() * mp2.getFock().getCP();
	Matrix hcore(N, N, 0.0);
	for (int p = 0; p < N; p++)
		for (int q = p; q < N; q++) {
			hcore(p, q) = hcore_spatial(p/2, q/2) * (p%2 == q%2);
			hcore(q, p) = hcore(p, q);
		}
	
	Tensor4& spinInts = mp2.getInts();
	
	// Assign the spin-Fock matrix
	spinFock.assign(N, N, 0.0);
	// Build 
	for (int p = 0; p < N; p++) {
		for (int q = p; q < N; q++) {
			spinFock(p, q) = hcore(p, q);
			for (int m = 0; m < nocc; m++) spinFock(p, q) += spinInts(p, m, q, m);
			spinFock(q, p) = spinFock(p, q);
		}
	}

	// Now build denominator arrays
	Dia.resize(nocc, N-nocc);
	Dijab.resize(nocc, nocc, N-nocc, N-nocc);
	for (int i = 0; i < nocc; i++)
		for (int a = 0; a < N-nocc; a++) {
		
			Dia(i, a) = spinFock(i, i) - spinFock(a+nocc, a+nocc);
			
			for (int j = i; j < nocc; j++)
				for (int b = a; b < N-nocc; b++)  {
					Dijab(i, j, a, b) = Dia(i, a) + spinFock(j, j) - spinFock(b+nocc, b+nocc);  
					Dijab(j, i, a, b) = Dijab(i, j, b, a) = Dijab(j, i, b, a) = Dijab(i, j, a, b);
				}
		}
}

void CCSD::build_guess() {
	// Assign singles amplitudes as zero
	singles.assign(nocc, N-nocc, 0.0);
	
	// Assign and build doubles amplitudes
	Tensor4& spinInts = mp2.getInts();
	Vector& eps = mp2.getFock().getEps();
	doubles.assign(nocc, nocc, N-nocc, N-nocc, 0.0);
	for (int i = 0; i < nocc; i++) {
		for (int j = i; j < nocc; j++) {
			auto eocc = eps[i/2] + eps[j/2];
			
			for (int a = nocc; a < N; a++)
				for (int b = a; b < N; b++) {
					auto resolvent = eocc - eps[a/2] - eps[b/2];
					doubles(i, j, a-nocc, b-nocc) = doubles(j, i, b-nocc, a-nocc) = spinInts(i, j, a, b) / resolvent;
					doubles(j, i, a-nocc, b-nocc) = doubles(i, j, b-nocc, a-nocc) = -doubles(i, j, a-nocc, b-nocc);
				}
			
		}
	}
	
	mp2.calculateEnergy(doubles);	
}

void CCSD::build_intermediates(Matrix& F, Tensor4& W, Tensor4& tau, Tensor4& tautilde) {
	// Construct tau tensors from amplitudes	
	for (int i = 0; i < nocc; i++)
		for (int j = i; j < nocc; j++)
			for (int a = 0; a < N-nocc; a++)
				for (int b = a; b < N-nocc; b++) {
					auto delta = singles(i, a) * singles(j, b) - singles(i, b) * singles(j, a);
					tau(i, j, a, b) = tau(j, i, b, a) = doubles(i, j, a, b) + delta;
					tau(j, i, a, b) = tau(i, j, b, a) = -tau(i, j, a, b);
					tautilde(i, j, a, b) = tautilde(j, i, b, a) = doubles(i, j, a, b) + 0.5*delta;
					tautilde(j, i, a, b) = tautilde(i, j, b, a) = -tautilde(i, j, a, b);
				}
					
	Tensor4& spinInts = mp2.getInts();		
	// Build F
	for (int a = nocc; a < N; a++) {
		for (int e = nocc; e < N; e++) {
			F(a, e) = a != e ? spinFock(a, e) : 0.0;
			double sum1 = 0.0, sum2 = 0.0, sum3 = 0.0;
			for (int m = 0; m < nocc; m++) {
				sum1 += spinFock(m, e) * singles(m, a-nocc);
				
				for (int f = nocc; f < N; f++) {
					sum2 += singles(m, f-nocc) * spinInts(m, a, f, e);
					
					for (int n = m; n < nocc; n++)
						sum3 += (2 - (m == n)) * tautilde(m, n, a-nocc, f-nocc) * spinInts(m, n, e, f);
				}
			}
			F(a, e) += sum2 - 0.5 * (sum1 + sum3);
		}
	}
	
	for (int m = 0; m < nocc; m++) {
		for (int i = 0; i < nocc; i++) {
			F(m, i) = m != i ? spinFock(m, i) : 0.0;
			double sum1 = 0.0, sum2 = 0.0, sum3 = 0.0;
			for (int e = nocc; e < N; e++) {
				sum1 += singles(i, e-nocc) * spinFock(m, e);
				
				for (int n = 0; n < nocc; n++) {
					sum2 += singles(n, e-nocc) * spinInts(m, n, i, e);
					
					for (int f = e; f < N; f++)
						sum3 += (2 - (f==e)) * tautilde(i, n, e-nocc, f-nocc) * spinInts(m, n, e, f);
				}
			}
			F(m, i) += 0.5*(sum1 + sum3) + sum2;
		}
	}	
	
	for (int m = 0; m < nocc; m++) {
		for (int e = nocc; e < N; e++) {
			F(m, e) = spinFock(m, e);
			
			double sum = 0.0;
			for (int n = 0; n < nocc; n++) {
				for (int f = nocc; f < N; f++)
					sum += singles(n, f-nocc) * spinInts(m, n, e, f);
			}
			
			F(m, e) += sum;
		}
	}
	
	// Build W
	
	W = spinInts;
	for (int m = 0; m < nocc; m++) {
		for (int n = 0; n < m; n++) {
			for (int i = 0; i < nocc; i++) {
				for (int j = 0; j < nocc; j++) {
					double sum = 0.0;
					for (int e = nocc; e < N; e++) {
						sum += singles(j, e-nocc) * spinInts(m, n, i, e);
						sum -= singles(i, e-nocc) * spinInts(m, n, j, e);
						
						double sum2 = 0.0;
						for (int f = nocc; f < N; f++) 
							sum2 += tau(i, j, e-nocc, f-nocc) * spinInts(m, n, e, f);
						sum += 0.25 * sum2; 
					}
					W(m, n, i, j) += sum;
					W(n, m, i, j) = - W(m, n, i, j);
				}
			}
		}
	} 
	
	for (int a = nocc; a < N; a++) {
		for (int b = nocc; b < N; b++) {
			for (int e = nocc; e < N; e++) {
				for (int f = nocc; f < e; f++) {
					double sum = 0.0;
					for (int m = 0; m < nocc; m++) {
						sum -= singles(m, b-nocc) * spinInts(a, m, e, f);
						sum += singles(m, a-nocc) * spinInts(b, m, e, f);
						
						double sum2 = 0.0;
						for (int n = 0; n < nocc; n++)
							sum2 += tau(m, n, a-nocc, b-nocc) * spinInts(m, n, e, f);
						sum += 0.25 * sum2;
					}
					W(a, b, e, f) += sum;
					W(a, b, f, e) = -W(a, b, e, f);
				}
				
			}
		}
	}
	
	for (int m = 0; m < nocc; m++) {
		for (int b = nocc; b < N; b++) {
			for (int e = nocc; e < N; e++) {
				for (int j = 0; j < nocc; j++) {
					double sum  = 0.0;
					for (int f = nocc; f < N; f++)
						sum += singles(j, f-nocc) * spinInts(m, b, e, f);
					
					for (int n = 0; n < nocc; n++) {
						sum -= singles(n, b-nocc) * spinInts(m, n, e, j);
						
						for (int f = nocc; f < N; f++)
							sum -= (0.5*doubles(j, n, f-nocc, b-nocc) + singles(j, f-nocc) * singles(n, b-nocc)) * spinInts(m, n, e, f);
					}
					
					W(m, b, e, j) += sum;
				}
			}
		}
	}

}	

void CCSD::build_amplitudes(Matrix& F, Tensor4& W, Tensor4& tau, Tensor4& tautilde) {
	Tensor4& spinInts = mp2.getInts();
	
	// Build new singles amplitudes
	Matrix newSingles(nocc, N-nocc, 0.0);
	for (int i = 0; i < nocc; i++) {
		for (int a = nocc; a < N; a++) {
			newSingles(i, a-nocc) = spinFock(i, a);
			
			double sum =0.0;
			for (int e = nocc; e < N; e++) sum += singles(i, e-nocc) * F(a, e);
			for (int m = 0; m < nocc; m++){
				 sum -= singles(m, a-nocc) * F(m, i);
				 for (int e = nocc; e < N; e++) {
					 sum += doubles(i, m, a-nocc, e-nocc) * F(m, e);
					 sum -= singles(m, e-nocc) * spinInts(m, a, i, e);
					 
					 double sum2 = 0.0;
					 for (int f = nocc; f < N; f++) sum2 += doubles(i, m, e-nocc, f-nocc) * spinInts(m, a, e, f);
					 for (int n = 0; n < nocc; n++) sum2 += doubles(m, n, a-nocc, e-nocc) * spinInts(n, m, e, i);
					 sum -= 0.5*sum2;
				 }
			 }
			 
			 newSingles(i, a-nocc) += sum;
			 newSingles(i, a-nocc) /= Dia(i, a-nocc);
		}
	}
	
	// Build new doubles amplitudes
	Tensor4 newDoubles(nocc, nocc, N-nocc, N-nocc, 0.0);
	for (int i = 0; i < nocc; i++) {
		for (int j = i; j < nocc; j++) {
			for (int a = nocc; a < N; a++) {
				for (int b = a; b < N; b++) {
					
					auto value = spinInts(i, j, a, b);
					
					double sum = 0.0;
					for (int e = nocc; e < N; e++) {
						double sum2 = 0.0, sum3 = 0.0;
						for (int m = 0; m < nocc; m++) {
							sum2 += singles(m, b-nocc) * F(m, e);
							sum3 += singles(m, a-nocc) * F(m, e);
						}
						sum += doubles(i, j, a-nocc, e-nocc) * (F(b, e) - 0.5*sum2);
						sum -= doubles(i, j, b-nocc, e-nocc) * (F(a, e) - 0.5*sum3);
					}
					
					for (int m = 0; m < nocc; m++) {
						double sum2 = 0.0, sum3 = 0.0;
						for (int e = nocc; e < N; e++) {
							sum2 += singles(j, e-nocc) * F(m, e);
							sum3 += singles(i, e-nocc) * F(m, e);
						}
						sum -= doubles(i, m, a-nocc, b-nocc) * (F(m, j) + 0.5*sum2);
						sum += doubles(j, m, a-nocc, b-nocc) * (F(m, i) + 0.5*sum3);
					}
					
					double sum2 = 0.0;
					for (int m = 0; m < nocc; m++) 
						for (int n = 0; n < nocc; n++)
							sum2 += tau(m, n, a-nocc, b-nocc) * W(m, n, i, j);
					
					for (int e = nocc; e < N; e++)
						for (int f = nocc; f < N; f++)
							sum2 += tau(i, j, e-nocc, f-nocc) * W(a, b, e, f);
					sum += 0.5*sum2;
					
					for (int m = 0; m < nocc; m++) {
						for (int e = nocc; e < N; e++) {
							sum +=  doubles(i, m, a-nocc, e-nocc) * W(m, b, e, j) - singles(i, e-nocc) * singles(m, a-nocc) * spinInts(m, b, e, j);
							sum += -doubles(j, m, a-nocc, e-nocc) * W(m, b, e, i) + singles(j, e-nocc) * singles(m, a-nocc) * spinInts(m, b, e, i);
							sum += -doubles(i, m, b-nocc, e-nocc) * W(m, a, e, j) + singles(i, e-nocc) * singles(m, b-nocc) * spinInts(m, a, e, j);
							sum +=  doubles(j, m, b-nocc, e-nocc) * W(m, a, e, i) - singles(j, e-nocc) * singles(m, b-nocc) * spinInts(m, a, e, i);
						}
					}
					
					for (int e = nocc; e < N; e++) {
						sum += singles(i, e-nocc) * spinInts(a, b, e, j);
						sum -= singles(j, e-nocc) * spinInts(a, b, e, i);
					}
					
					for (int m = 0; m < nocc; m++) {
						sum -= singles(m, a-nocc) * spinInts(m, b, i, j);
						sum += singles(m, b-nocc) * spinInts(m, a, i, j);
					}
					
					
					value += sum;
					value /= Dijab(i, j, a-nocc, b-nocc);
					newDoubles(i, j, a-nocc, b-nocc) = value;
					newDoubles(j, i, a-nocc, b-nocc) = -value;
					newDoubles(i, j, b-nocc, a-nocc) = -value;
					newDoubles(j, i, b-nocc, a-nocc) = value;
				
				}
			}
		}
	}
	
	// Compute RMS differences
	delta_singles = fnorm(newSingles - singles);
	singles = newSingles;
	delta_doubles = fnorm(newDoubles - doubles);
	doubles = newDoubles;
}

void CCSD::calculateEnergy() { 
	double newEnergy = 0.0;
	
	Tensor4& spinInts = mp2.getInts();
	
	double sum1 = 0.0, sum2 = 0.0;
	for (int i = 0; i < nocc; i++) {
		for (int a = nocc; a < N; a++) {
			newEnergy += spinFock(i, a) * singles(i, a-nocc);
			
			for (int  j = 0; j < nocc; j++) {
				for (int b = nocc; b < N; b++) {
					sum1 += spinInts(i, j, a, b) * doubles(i, j, a-nocc, b-nocc);
					sum2 += spinInts(i, j, a, b) * singles(i, a-nocc) * singles(j, b-nocc);
				}
			}
		}
	}
	newEnergy += 0.25*sum1 + 0.5*sum2; 
	
	delta_e = newEnergy - energy;
	energy = newEnergy; 
}

void CCSD::calculateTriples() {
	triples_energy = 0.0; 
	
	int nvirt = N-nocc;
	int nvirt2 = nvirt * nvirt;
	double Dijk[nvirt * nvirt2];
	double tijkd[nvirt * nvirt2];
	double tijkc[nvirt * nvirt2];
	
	Tensor4& spinInts = mp2.getInts();
	int A, B, C, a2x, a1x, b2x, b1x, c2x, c1x;
	for (int i = 0; i < nocc; i++)
		for (int j = 0; j < nocc; j++)
			for (int k = 0; k < nocc; k++) {
				auto fijk = spinFock(i, i) + spinFock(j, j) + spinFock(k, k);
				
				// Build Dijk and connected/disconnected tijk
				for (int a = 0; a < nvirt; a++) {
					a2x = a*nvirt2; a1x = a*nvirt; A = a+nocc;
					for (int b = 0; b <= a; b++) {
						b2x = b*nvirt2; b1x = b*nvirt; B = b+nocc;
						for (int c = 0; c <= b; c++) {
							c2x = c*nvirt2; c1x = c*nvirt; C = c+nocc;
							auto val = fijk - spinFock(A, A) - spinFock(B, B) - spinFock(C, C);
							Dijk[a2x + b1x + c] = Dijk[a2x + c1x + b] = Dijk[b2x + a1x + c] = val;
							Dijk[b2x + c1x + a] = Dijk[c2x + a1x + b] = Dijk[c2x + b1x + a] = val;
							
							auto val2 = singles(i, a) * spinInts(j, k, B, C) - singles(j, a) * spinInts(i, k, B, C)
									- singles(k, a) * spinInts(j, i, B, C) - singles(i, b) * spinInts(j, k, A, C)
										+ singles(j, b) * spinInts(i, k, A, C) + singles(k, b) * spinInts(j, i, A, C)
											- singles(i, c) * spinInts(j, k, B, A) + singles(j, c) * spinInts(i, k, B, A)
												+ singles(k, c) * spinInts(j, i, B, A);
							val2 /= val;
							tijkd[a2x + b1x + c] = tijkd[b2x + c1x + a] = tijkd[c2x + a1x + b] = val2;
							tijkd[a2x + c1x + b] = tijkd[b2x + a1x + c] = tijkd[c2x + b1x + a] = -val2;
							
							auto val3 = 0.0;
							int E;
							for (int e = 0; e < nvirt; e++) {
								E = e + nocc;
								val3 += doubles(j, k, a, e)*spinInts(E, i, B, C) - doubles(i, k, a, e)*spinInts(E, j, B, C)
										- doubles(j, i, a, e)*spinInts(E, k, B, C) - doubles(j, k, b, e)*spinInts(E, i, A, C)
											+ doubles(i, k, b, e)*spinInts(E, j, A, C) + doubles(j, i, b, e)*spinInts(E, k, A, C)
												- doubles(j, k, c, e)*spinInts(E, i, B, A) + doubles(i, k, c, e)*spinInts(E, j, B, A)
													+ doubles(j, i, c, e)*spinInts(E, k, B, A);
							}
							for (int m = 0; m < nocc; m++) {
								val3 += -doubles(i, m, b, c)*spinInts(m, A, j, k) + doubles(j, m, b, c)*spinInts(m, A, i, k)
										+ doubles(k, m, b, c)*spinInts(m, A, j, i) + doubles(i, m, a, c)*spinInts(m, B, j, k)
											- doubles(j, m, a, c)*spinInts(m, B, i, k) - doubles(k, m, a, c)*spinInts(m, B, j, i)
												+ doubles(i, m, b, a)*spinInts(m, C, j, k) - doubles(j, m, b, a)*spinInts(m, C, i, k)
													- doubles(k, m, b, a)*spinInts(m, C, j, i);
							}
							
							val3 /= val;
							tijkc[a2x + b1x + c] = tijkc[b2x + c1x + a] = tijkc[c2x + a1x + b] = val3;
							tijkc[a2x + c1x + b] = tijkc[b2x + a1x + c] = tijkc[c2x + b1x + a] = -val3;
						}
					}
				}
				
				// Determine contribution to triples energy
				for (int a = 0; a < nvirt; a++) {
					a2x = a * nvirt2;
					for (int b = 0; b < nvirt; b++) {
						b1x = b *nvirt;
						for (int c = 0; c < nvirt; c++) {
							c1x = a2x + b1x + c;
							triples_energy += tijkc[c1x] * Dijk[c1x] * ( tijkc[c1x] + tijkd[c1x]);
						}
					}
				}
							
			}
			
	triples_energy /= 36.0;
				
}

void CCSD::compute()
{
	Logger& log = mp2.getFock().getMolecule().getLog();
	
	log.title("CCSD CALCULATION");
	
	log.print("Transforming MO integrals to spin basis\n");
	mp2.spatialToSpin(); 
	log.localTime();
	
	log.print("Building guess amplitudes and spin-basis Fock matrix\n");
	build_guess();
	build_fock();
	log.localTime();

	mp2.calculateEnergy(doubles);
	log.print("MP2 Energy = " + std::to_string(mp2.getEnergy()) + "\n");
	log.localTime();
	
	log.print("\nBeginning CC iterations:\n\n");
	bool converged = false; 
	int MAXITER = log.maxiter();
	int iter = 1;
	Matrix F(N, N, 0.0);
	Tensor4 W(N, N, N, N, 0.0);
	Tensor4 tau(nocc, nocc, N-nocc, N-nocc, 0.0);
	Tensor4 tautilde(nocc, nocc, N-nocc, N-nocc, 0.0);
	
	double time_interms, time_amps, time_en;
	while (!converged && iter < MAXITER) {
		
		build_intermediates(F, W, tau, tautilde);
		time_interms = log.getLocalTime();
		build_amplitudes(F, W, tau, tautilde);
		time_amps = log.getLocalTime();
		calculateEnergy();
		time_en = log.getLocalTime();
		
		log.iteration(iter, energy, delta_e, delta_doubles);
		std::cout << "Iteration " << iter << ": " << delta_singles << " " << time_interms << " " << time_amps << " " << time_en << std::endl;
		converged = (delta_e < log.converge()) && (delta_doubles < log.converge());
		iter++;
	}
	
	if (converged){
		log.result("CCSD correlation energy", energy, "Hartree");
		if (withTriples) {
			log.print("Calculating triples... \n");
			calculateTriples();
			log.localTime();
			log.result("(T) correction", triples_energy, "Hartree");
			log.result("CCSD(T) correlation energy", energy + triples_energy, "Hartree");
		}
	} else log.result("CCSD failed to converge");
		
}

					
