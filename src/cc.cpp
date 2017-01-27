
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
CCSD::CCSD(MP2& _mp2) : mp2(_mp2) {
	N = mp2.getN();
	nocc = mp2.getNocc();
	energy = 0.0;
	delta_e = 0.0;
	delta_singles = 0.0;
	delta_doubles = 0.0;
}

void CCSD::build_fock() {
	Matrix& hcore = mp2.getFock().getHCore();
	Tensor4& spinInts = mp2.getInts();
	
	// Assign the spin-Fock matrix
	spinFock = hcore;
	// Build 
	for (int p = 0; p < N; p++)
		for (int q = 0; q < N; q++) 
			for (int m = 0; m < nocc; m++) spinFock(p, q) += spinInts(p, m, q, m);
	
	// Now build denominator arrays
	Dia.resize(nocc, N-nocc);
	Dijab.resize(nocc, nocc, N-nocc, N-nocc);
	for (int i = 0; i < nocc; i++)
		for (int a = 0; a < N-nocc; a++) {
		
			Dia(i, a) = spinFock(i, i) - spinFock(a+nocc, a+nocc);
			
			for (int j = 0; j < nocc; j++)
				for (int b = 0; b < N-nocc; b++) 
					Dijab(i, j, a, b) = Dia(i, a) + spinFock(j, j) - spinFock(b+nocc, b+nocc);  
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
		for (int j = 0; j < nocc; j++) {
			auto eocc = eps[i] + eps[j];
			
			for (int a = nocc; a < N; a++)
				for (int b = nocc; b < N; b++) {
					auto resolvent = eocc - eps[a] - eps[b];
					doubles(i, j, a - nocc, b-nocc) = spinInts(i, j, a, b) / resolvent;
				}
			
		}
	}
	
	mp2.calculateEnergy(doubles);	
}

void CCSD::build_intermediates(Matrix& F, Tensor4& W, Tensor4& tau, Tensor4& tautilde) {
	// Construct tau tensors from amplitudes	
	for (int i = 0; i < nocc; i++)
		for (int j = 0; j < nocc; j++)
			for (int a = 0; a < N-nocc; a++)
				for (int b = 0; b < N-nocc; b++) {
					auto delta = singles(i, a) * singles(j, b) - singles(i, b) * singles(j, a);
					tau(i, j, a, b) = doubles(i, j, a, b) + delta;
					tautilde(i, j, a, b) = doubles(i, j, a, b) + 0.5*delta;
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
					
					for (int n = 0; n < nocc; n++)
						sum3 += tautilde(m, n, a-nocc, f-nocc) * spinInts(m, n, e, f);
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
					
					for (int f = nocc; f < N; f++)
						sum3 += tautilde(i, n, e-nocc, f-nocc);
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
		for (int n = 0; n < nocc; n++) {
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
				}
			}
		}
	}
	
	for (int a = nocc; a < N; a++) {
		for (int b = nocc; b < N; b++) {
			for (int e = nocc; e < N; e++) {
				for (int f = nocc; f < N; f++) {
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
		for (int j = 0; j < nocc; j++) {
			for (int a = nocc; a < N; a++) {
				for (int b = nocc; b < N; b++) {
					
					newDoubles(i, j, a-nocc, b-nocc) = spinInts(i, j, a, b);
					
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
							sum += doubles(i, m, a-nocc, e-nocc) * W(m, b, e, j) - singles(i, e-nocc) * singles(m, a-nocc) * spinInts(m, b, e, j);
							sum -= doubles(j, m, a-nocc, e-nocc) * W(m, b, e, i) - singles(j, e-nocc) * singles(m, a-nocc) * spinInts(m, b, e, i);
							sum -= doubles(i, m, b-nocc, e-nocc) * W(m, a, e, j) - singles(i, e-nocc) * singles(m, b-nocc) * spinInts(m, a, e, j);
							sum += doubles(j, m, b-nocc, e-nocc) * W(m, a, e, i) - singles(j, e-nocc) * singles(m, b-nocc) * spinInts(m, a, e, i);
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
					
					newDoubles(i, j, a-nocc, b-nocc) += sum;
					newDoubles(i, j, a-nocc, b-nocc) /= Dijab(i, j, a-nocc, b-nocc);
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
	
	for (int i = 0; i < nocc; i++) {
		for (int a = nocc; a < N; a++) {
			newEnergy += spinFock(i, a) * singles(i, a-nocc);
			
			double sum1 = 0.0, sum2 = 0.0;
			for (int  j = 0; j < nocc; j++) {
				for (int b = nocc; b < N; b++) {
					sum1 += spinInts(i, j, a, b) * doubles(i, j, a-nocc, b-nocc);
					sum2 += spinInts(i, j, a, b) * singles(i, a-nocc) * singles(j, b-nocc);
				}
			}
			newEnergy += 0.25 * sum1 + 0.5 * sum2;
		}
	}
	
	delta_e = newEnergy - energy;
	energy = newEnergy; 
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
	while (!converged && iter < MAXITER) {
		
		build_intermediates(F, W, tau, tautilde);
		build_amplitudes(F, W, tau, tautilde);
		calculateEnergy();
		
		log.iteration(iter, energy, delta_e, delta_doubles);
		std::cout << "Iteration " << iter << ": " << delta_singles << std::endl;
		converged = (delta_e < log.converge()) && (delta_doubles < log.converge());
		iter++;
	}
}

					
