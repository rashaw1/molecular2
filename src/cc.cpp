
#include "cc.hpp"
#include "logger.hpp"
#include "fock.hpp"
#include "molecule.hpp"
#include "integrals.hpp"
#include "error.hpp"
#include <iostream>
#include <thread>

// Constructor
CCSD::CCSD(Command& c, MP2& _mp2) : cmd(c), mp2(_mp2) {
	N = mp2.getN();
	nocc = mp2.getNocc(); 
	energy = 0.0;
	triples_energy = 0.0;
	delta_e = 0.0;
	delta_singles = 0.0;
	delta_doubles = 0.0;
	
	maxDiis = cmd.get_option<int>("maxdiis");
	doDiis = cmd.get_option<bool>("diis");
	withTriples = cmd.get_option<bool>("triples"); 
	
	diis.init(maxDiis, doDiis);
}

double divide(double a, double b){
	return a/b;
}

/*
void CCSD::build_intermediates(Matrix& F, Tensor4& W, S4OddTensor4& tau, S4OddTensor4& tautilde) {
	// Construct tau tensors from amplitudes	
	for (int i = 0; i < nocc; i++)
		for (int j = 0; j <= i; j++)
			for (int a = 0; a < N-nocc; a++)
	for (int b = 0; b <= a; b++) {
		auto delta = singles(i, a) * singles(j, b) - singles(i, b) * singles(j, a);
		tau.set(i, j, a, b, doubles(i, j, a, b) + delta );
		tautilde.set(i, j, a, b, doubles(i, j, a, b) + 0.5*delta );
	}
					
	S8OddTensor4& spinInts = mp2.getSpinInts();		
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
			std::cout << m << " " << i << " " << F(m, i) << std::endl;
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
	
	for (int m = 0; m < nocc; m++) {
		for (int n = 0; n < m; n++) {
			for (int i = 0; i < nocc; i++) {
				for (int j = 0; j < nocc; j++) {
					W(m, n, i, j) = spinInts(m, n, i, j);
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
					W(a, b, e, f) = spinInts(a, b, e, f);
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
					W(m, b, e, j) = spinInts(m, b, e, j);
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

void CCSD::build_amplitudes(Matrix& F, Tensor4& W, S4OddTensor4& tau, S4OddTensor4& tautilde) {
	S8OddTensor4& spinInts = mp2.getSpinInts();
	
	// Build new singles amplitudes
	Matrix newSingles = Matrix::Zero(nocc, N-nocc);
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
	S4OddTensor4 newDoubles(nocc, N-nocc, 0.0);
	for (int i = 0; i < nocc; i++) {
		for (int j = 0; j <= i; j++) {
			for (int a = nocc; a < N; a++) {
				for (int b = nocc; b <= a; b++) {
					
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
					newDoubles.set(i, j, a-nocc, b-nocc, value);
				
				}
			}
		}
	}
	
	// Compute RMS differences
	calculateError(newSingles, newDoubles);
	singles = newSingles;
	doubles = newDoubles;
}*/

void CCSD::calculateError(Matrix& newSingles, S4OddTensor4& newDoubles) {
/*	Matrix ds = newSingles - singles;
	S4OddTensor4 dd = newDoubles -  doubles;
	
	delta_singles = ds.norm();
	delta_doubles = fnorm(dd);
	if (doDiis) {
		// Save amplitudes
		if (singles_cache.size() == maxDiis) {
			singles_cache.erase(singles_cache.begin());
			doubles_cache.erase(doubles_cache.begin());
		}
		singles_cache.push_back(newSingles);
		doubles_cache.push_back(newDoubles);
		
		// Compute error vector
		std::vector<Vector> errs;
		Vector new_err(Eigen::Map<Vector>(ds.data(), ds.cols()*ds.rows()));
		errs.push_back(new_err);
		errs.push_back(tensorToVector(dd, Tensor4::ODD_4));
		Vector weights = diis.compute(errs);
		if (iter > 2) {
			newSingles = Matrix::Zero(nocc, N-nocc);
			newDoubles.assign(nocc, N-nocc, 0.0);
			int offset = singles_cache.size() - weights.size();
			for (int i = offset; i < singles_cache.size(); i++) {
				newSingles = newSingles + weights[i-offset]*singles_cache[i];
				newDoubles = newDoubles + weights[i-offset]*doubles_cache[i];
			} 
		}
		
	}
	*/
}

void slow_iteration(Vector& eps, S8EvenTensor4& moInts, Matrix& Tai, Tensor4& Tabij, int nocc, int nvirt) {
	
	Tensor4 Tau(nvirt, nvirt, nocc, nocc, 0.0); 
	for (int a = 0; a < nvirt; a++)
		for (int b = 0; b < nvirt; b++)
			for (int i = 0; i < nocc; i++)
				for (int j = 0; j < nocc; j++)
					Tau(a, b, i, j) = Tabij(a, b, i, j) + Tai(a, i)*Tai(b, j); 
	
	Matrix Fia = Matrix::Zero(nocc, nvirt);
	Matrix Fab = Matrix::Zero(nvirt, nvirt);
	Matrix Fij = Matrix::Zero(nocc, nocc);
	Matrix Fpij = Matrix::Zero(nocc, nocc); 
	
	for (int i = 0; i < nocc; i++) {
		for (int a = 0; a < nvirt; a++) {
			for (int m = 0; m < nocc; m++) {
				for (int e = 0; e < nvirt; e++) {
					Fia(i, a) += (2.0*moInts(i, a+nocc, m, e+nocc) - moInts(i, e+nocc, m, a+nocc)) * Tai(e, m); 
				}
			}
		}
	}
	
	
	for (int a = 0; a < nvirt; a++) {
		for (int b = 0; b < nvirt; b++) {
			for (int m = 0; m < nocc; m++) {
				for (int e = 0; e < nvirt; e++) {
					Fab(a, b) += (2.0 * moInts(a+nocc, b+nocc, m, e+nocc) - moInts(m, b+nocc, a+nocc, e+nocc))*Tai(e, m);
					
					for (int n = 0; n < nocc; n++)
						Fab(a, b) -= (2.0*moInts(m, e+nocc, n, b+nocc) - moInts(m, b+nocc, n, e+nocc)) * Tau(e, a, m, n); 
				}
			}
		}
	}
	
	for (int i = 0; i < nocc; i++) {
		for (int j = 0; j < nocc; j++) {
			for (int m = 0; m < nocc; m++) {
				for (int e = 0; e < nvirt; e++) {
					Fpij(i, j) += (2.0*moInts(i, j, m, e+nocc) - moInts(i, e+nocc, m, j))*Tai(e, m); 
					
					for (int f = 0; f < nvirt; f++)
						Fpij(i, j) += (2.0*moInts(m, e+nocc, i, f+nocc) - moInts(i, e+nocc, m, f+nocc))*Tabij(e, f, m, j); 
				}
			}
			
			Fij(i, j) = Fpij(i, j);
			for (int e = 0; e < nvirt; e++)
				Fij(i, j) += Fia(i, e) * Tai(e, j); 
		
		}
	}
	
	Tensor4 Wijkl(nocc, nocc, nocc, nocc, 0.0);
	Tensor4 Wiajb(nocc, nvirt, nocc, nvirt, 0.0);
	Tensor4 Wiabj(nocc, nvirt, nvirt, nocc, 0.0); 
	Tensor4 Wabci(nvirt, nvirt, nvirt, nocc, 0.0);
	Tensor4 Wiajk(nocc, nvirt, nocc, nocc, 0.0); 
	
	for (int i = 0; i < nocc; i++) {
		for (int j = 0; j < nocc; j++) {
			for (int k = 0; k < nocc; k++) {
				for (int l = 0; l < nocc; l++) {
					
					Wijkl(i, j, k, l) = moInts(i, k, j, l); 
					
					for (int e = 0; e < nvirt; e++) {
						for (int f = 0; f < nvirt; f++)
							Wijkl(i, j, k, l) += moInts(i, e+nocc, j, f+nocc) * Tau(e, f, k, l); 
						
						Wijkl(i, j, k, l) += Tai(e, k) * moInts(i, e+nocc, j, l); 
						Wijkl(i, j, k, l) += Tai(e, l) * moInts(j, e+nocc, i, k); 
					}
					
				}
			}
		}
	}
	
	for (int i = 0; i < nocc; i++) {
		for (int a = 0; a < nvirt; a++) {
			for (int j = 0; j < nocc; j++) {
				for (int b = 0; b < nvirt; b++) {
					Wiajb(i, a, j, b) = moInts(i, j, a+nocc, b+nocc); 
					Wiabj(i, a, b, j) = moInts(i, b+nocc, a+nocc, j); 
					
					
					for (int m = 0; m < nocc; m++) {
						for (int e = 0; e < nvirt; e++) {
							Wiajb(i, a, j, b) -= 0.5*moInts(i, e+nocc, m, b+nocc)*(Tau(e, a, j, m) + Tai(e, j)*Tai(a, m)); 
							Wiabj(i, a, b, j) -= 0.5*moInts(i, b+nocc, m, e+nocc)*(Tau(a, e, m, j) + Tai(a, m)*Tai(e, j)); 
							Wiabj(i, a, b, j) += 0.5*(2.0*moInts(i, b+nocc, m, e+nocc) - moInts(i, e+nocc, m, b+nocc))*Tabij(e, a, m, j); 
						}
						
						Wiajb(i, a, j, b) -= moInts(i, j, m, b+nocc) * Tai(a, m); 
						Wiabj(i, a, b, j) -= moInts(i, b+nocc, m, j) * Tai(a, m); 
					}
					
					for (int e = 0; e < nvirt; e++){
						Wiajb(i, a, j, b) += moInts(i, e+nocc, a+nocc, b+nocc)*Tai(e, j); 
						Wiabj(i, a, b, j) += moInts(i, b+nocc, a+nocc, e+nocc)*Tai(e, j); 
					}
					
				}
			}
		}
	}
	
	for (int a = 0; a < nvirt; a++) {
		for (int b = 0; b < nvirt; b++) {
			for (int c = 0; c < nvirt; c++) {
				for (int i = 0; i < nocc; i++) {
					Wabci(a, b, c, i) = moInts(a+nocc, c+nocc, b+nocc, i); 
					
					for (int m = 0; m < nocc; m++) {
						Wabci(a, b, c, i) -= moInts(a+nocc, c+nocc, m, i) * Tai(b, m); 
						Wabci(a, b, c, i) -= Tai(a, m) * moInts(m, c+nocc, b+nocc, i); 
					}
					
				}
			}
		}
	}
	
	for (int i = 0; i < nocc; i++) {
		for (int a = 0; a < nvirt; a++) {
			for (int j = 0; j < nocc; j++) {
				for (int k = 0; k < nocc; k++) {
					Wiajk(i, a, j, k) = moInts(i, j, a+nocc, k); 
					
					for (int e = 0; e < nvirt; e++) {
						for (int f = 0; f < nvirt; f++) {
							Wiajk(i, a, j, k) += moInts(i, e+nocc, a+nocc, f+nocc)*Tau(e, f, j, k); 
						}
					}
				}
			}
		}
	}
	
	Matrix Zai = Matrix::Zero(nvirt, nocc); 
	Tensor4 Zabij(nvirt, nvirt, nocc, nocc, 0.0); 
	
	for (int a = 0; a < nvirt; a++) {
		for (int i = 0; i < nocc; i++) {
			
			for (int e = 0; e < nvirt; e++) {
				Zai(a, i) += Fab(a, e)*Tai(e, i); 
				
				for (int m = 0; m < nocc; m++) {
					Zai(a, i) += Fia(m, e) * (2.0*Tabij(e, a, m, i) - Tabij(e, a, i, m)); 
					Zai(a, i) += Tai(e, m) * (2.0*moInts(m, e+nocc, a+nocc, i) - moInts(a+nocc, e+nocc, m, i)); 
					
					for (int n = 0; n < nocc; n++) 
						Zai(a, i) -= moInts(m, e+nocc, n, i) * (2.0 * Tabij(e, a, m, n) - Tabij(a, e, m, n)); 
					
					for (int f = 0; f < nvirt; f++) 
						Zai(a, i) += moInts(m, e+nocc, a+nocc, f+nocc) * (2.0 * Tabij(e, f, m, i) - Tabij(e, f, i, m));
					
				}
			}
			
			for (int m = 0; m < nocc; m++)
				Zai(a, i) -= Fpij(m, i) * Tai(a, m); 
			
		}
	}
	
	for (int a = 0; a < nvirt; a++) {
		for (int b = 0; b < nvirt; b++) {
			for (int i = 0; i < nocc; i++) {
				for (int j = 0; j < nocc; j++) {
					
					Zabij(a, b, i, j) = moInts(a+nocc, i, b+nocc, j); 
					
					for (int e = 0; e < nvirt; e++) {
						Zabij(a, b, i, j) += Tabij(a, e, i, j) * Fab(b, e); 
						Zabij(b, a, j, i) += Tabij(b, e, j, i) * Fab(a, e); 
						
						for (int f = 0; f < nvirt; f++) {
							Zabij(a, b, i, j) += 0.5*moInts(a+nocc, e+nocc, b+nocc, f+nocc) * Tau(e, f, i, j); 
							Zabij(a, b, i, j) += 0.5*moInts(b+nocc, e+nocc, a+nocc, f+nocc) * Tau(e, f, j, i);
						}
						
						for (int m = 0; m < nocc; m++) {
							Zabij(a, b, i, j) -= Tabij(a, e, m, j) * Wiajb(m, b, i, e); 
							Zabij(a, b, i, j) -= Tabij(b, e, m, i) * Wiajb(m, a, j, e); 
							
							Zabij(a, b, i, j) -= Wiajb(m, a, i, e) * Tabij(e, b, m, j); 
							Zabij(a, b, i, j) -= Wiajb(m, b, j, e) * Tabij(e, a, m, i); 
							
							Zabij(a, b, i, j) += (2.0*Tabij(e, a, m, i) - Tabij(e, a, i, m))*Wiabj(m, b, e, j); 
							Zabij(a, b, i, j) += (2.0*Tabij(e, b, m, j) - Tabij(e, b, j, m))*Wiabj(m, a, e, i); 
						}
						
						Zabij(a, b, i, j) += Tai(e, i)*Wabci(a, b, e, j); 
						Zabij(a, b, i, j) += Tai(e, j)*Wabci(b, a, e, i); 
					} 
					
					for (int m = 0; m < nocc; m++) {
						Zabij(a, b, i, j) -= Tabij(a, b, i, m) * Fij(m, j); 
						Zabij(a, b, i, j) -= Tabij(b, a, j, m) * Fij(m, i); 
						
						for (int n = 0; n < nocc; n++) {
							Zabij(a, b, i, j) += 0.5*Tau(a, b, m, n) * Wijkl(m, n, i, j); 
							Zabij(a, b, i, j) += 0.5*Tau(b, a, m, n) * Wijkl(m, n, j, i); 
						}
						
						Zabij(a, b, i, j) -= Tai(a, m) * Wiajk(m, b, i, j); 
						Zabij(a, b, i, j) -= Tai(b, m) * Wiajk(m, a, j, i); 			
					}

				}
			}
		}
	}
	
	double sres, dres; 
	for (int a = 0; a < nvirt; a++) {
		for (int i = 0; i < nocc; i++) {
			sres = eps[i] - eps[a+nocc]; 
			Tai(a, i) = Zai(a, i) / sres; 
			
			for (int b = 0; b < nvirt; b++) {
				for (int j = 0; j < nocc; j++) {
					dres = sres + eps[j] - eps[b+nocc]; 
					
					Tabij(a, b, i, j) = Zabij(a, b, i, j) / dres;  
				}
			}
		}
	}
	
	for (int a = 0; a < nvirt; a++)
		for (int b = 0; b < nvirt; b++)
			for (int i = 0; i < nocc; i++)
				for (int j = 0; j < nocc; j++)
					Tau(a, b, i, j) = Tabij(a, b, i, j) + Tai(a, i)*Tai(b, j); 
	
	double en = 0.0; 
	for (int a = 0; a < nvirt; a++)
		for (int b = 0; b < nvirt; b++)
			for (int i = 0; i < nocc; i++)
	for (int j = 0; j < nocc; j++) {
					en += Tau(a, b, i, j) * (2.0*moInts(i, a+nocc, j, b+nocc) - moInts(i, b+nocc, j, a+nocc)); 
					if (fabs(Tabij(a, b, i, j)) > 0.05) std::cout << a << " " << b << " " << i << " " << j << " " << Tabij(a, b, i, j) << std::endl;
				}
	std::cout << en << std::endl; 

}

void CCSD::compute()
{
	Logger& log = mp2.getFock().getMolecule()->control->log;
	
	log.title("CCSD CALCULATION");
	
	log.print("MP2 Energy = " + std::to_string(mp2.getEnergy()) + "\n");
	log.localTime();
			
	log.print("\nBeginning CC iterations:\n\n");
	log.flush();
	bool converged = false; 
	int MAXITER = cmd.get_option<int>("maxiter");
	double CONVERGE = cmd.get_option<double>("converge");
	iter = 1;
		
	Amplitudes& T = *(mp2.getAmplitudes());
	Integrals& V = *(mp2.getSpinInts()); 
		
	log.initIterationCC();
	double time_tot; 
	double new_en; 
	double nt1 = T.ai->norm2();
	T["abij"] = T["abij"]; 
	double nt2 = T.abij->norm2();
	double newnt1, newnt2; 
	
	int oovv[4] = {nocc, nocc, N-nocc, N-nocc};
	int vvoo[4] = {N-nocc, N-nocc, nocc, nocc};
	int nsns[4] = {NS, NS, NS, NS}; 
	CTF::Tensor<> Aijab(4, oovv, nsns, mp2.dw);
	Aijab["ijab"] = 2.0 * V["iajb"] - V["ibja"]; 
	CTF::Tensor<> Tau(4, vvoo, nsns, mp2.dw);
	
	/*mp2.transformIntegrals(); 
	S8EvenTensor4& moInts = mp2.getMOInts(); 
	Vector& eps = mp2.getFock().getEps(); 
	Matrix Tai = Matrix::Zero(N-nocc, nocc); 
	Tensor4 Tabij(N-nocc, N-nocc, nocc, nocc, 0.0); 
	for (int a = 0; a < N-nocc; a++) {
		for (int b = 0; b < N-nocc; b++) {
			for (int i = 0; i < nocc; i++) {
				for (int j = 0; j < nocc; j++) {
					Tabij(a, b, i, j) = moInts(i, a+nocc, j, b+nocc) / (eps[i] + eps[j] - eps[a+nocc] - eps[b+nocc]); 
				}
			}
		}
	}
	
	double mp2en = 0.0;
	for (int a = 0; a < N-nocc; a++) {
		for (int b = 0; b < N-nocc; b++) {
			for (int i = 0; i < nocc; i++) {
				for (int j = 0; j < nocc; j++) {
					mp2en += Tabij(a, b, i, j) * (2.0*moInts(i, a+nocc, j, b+nocc) - moInts(i, b+nocc, j, a+nocc)); 
				}
			}
		}
	}
	std::cout << mp2en << std::endl;  */
	
	while (!converged && iter < MAXITER) {
		ccsd_iteration(V, T, 0); 
		//slow_iteration(eps, moInts, Tai, Tabij, nocc, N-nocc);  
		Tau["abij"] = T["abij"] + T["ai"]*T["bj"];
		new_en = Aijab["ijab"]*Tau["abij"]; 

		delta_e = new_en - energy;
		energy = new_en; 
		newnt1 = T.ai->norm2();
		newnt2 = T.abij->norm2();
		time_tot = log.getLocalTime();
		delta_singles = nt1 - newnt1;
		delta_doubles = nt2 - newnt2;
		nt1 = newnt1;
		nt2 = newnt2; 
		log.iterationCC(iter, energy, delta_e, delta_singles, delta_doubles, -1, -1, time_tot); 
		converged = (fabs(delta_e) < CONVERGE) && (fabs(delta_doubles) < CONVERGE); 
		iter++;
	}
	
	if (converged){
		log.result("CCSD correlation energy", energy, "Hartree");
		if (withTriples) {
			log.print("Calculating triples... \n");
			calculateTriples(V, T);
			log.localTime();
			log.result("(T) correction", triples_energy, "Hartree");
			log.result("CCSD(T) correlation energy", energy + triples_energy, "Hartree");
		}
	} else log.result("CCSD failed to converge");
		
}

void CCSD::ccsd_iteration(Integrals   &V,
Amplitudes  &T,
int sched_nparts){
	int rank;   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	int nvirt = N - nocc; 
	
	// Intermediates
	
	CTF::Matrix<> FIA(nocc, nvirt, NS, mp2.dw); 
	CTF::Matrix<> FAB(nvirt, nvirt, NS, mp2.dw);
	CTF::Matrix<> FIJ(nocc, nocc, NS, mp2.dw); 
	CTF::Matrix<> FPIJ(nocc, nocc, NS, mp2.dw);
	
	int oooo[4] = {nocc, nocc, nocc, nocc}; 
	int ovov[4] = {nocc, nvirt, nocc, nvirt};
	int vvvo[4] = {nvirt, nvirt, nvirt, nocc};
	int ovoo[4] = {nocc, nvirt, nocc, nocc};
	int vvoo[4] = {nvirt, nvirt, nocc, nocc}; 
	int ovvo[4] = {nocc, nvirt, nvirt, nocc}; 
	
	int sysy[4] = {SY, NS, SY, NS};
	int nssy[4] = {NS, NS, SY, NS};
	int syns[4] = {SY, NS, NS, NS};
	int nsns[4] = {NS, NS, NS, NS};
	
	CTF::Tensor<> Tau(4, vvoo, nsns, mp2.dw);
	Tau["abij"] = T["abij"];
	Tau["abij"] += T["ai"]*T["bj"]; 
	
	CTF::Tensor<> WIJKL(4, oooo, nsns, mp2.dw);
	CTF::Tensor<> WIAJB(4, ovov, nsns, mp2.dw);
	CTF::Tensor<> WIABJ(4, ovvo, nsns, mp2.dw);
	CTF::Tensor<> WABCI(4, vvvo, nsns, mp2.dw);
	CTF::Tensor<> WIAJK(4, ovoo, nsns, mp2.dw); 
	
	FIA["ia"] = 2.0*V["iame"]*T["em"];
	FIA["ia"] -= V["iema"]*T["em"];
	
	FAB["ab"] = 2.0*V["abme"]*T["em"]; 
	FAB["ab"] -= V["mbae"]*T["em"]; 
	FAB["ab"] -= 2.0*V["menb"]*Tau["eamn"]; 
	FAB["ab"] += V["mbne"]*Tau["eamn"]; 
	
	FPIJ["ij"] = 2.0*V["ijme"]*T["em"];
	FPIJ["ij"] -= V["iemj"]*T["em"];
	FPIJ["ij"] += 2.0*V["meif"]*T["efmj"]; 
	FPIJ["ij"] -= V["iemf"]*T["efmj"];
	
	FIJ["ij"] = FPIJ["ij"];
	FIJ["ij"] += FIA["ie"]*T["ej"]; 
	
	WIJKL["ijkl"] = V["ikjl"];
	WIJKL["ijkl"] += V["iejf"]*Tau["efkl"]; 
	WIJKL["ijkl"] += T["ek"]*V["iejl"];
	WIJKL["ijkl"] += T["el"]*V["jeik"];
	
	WIAJB["iajb"] = V["ijab"];
	WIAJB["iajb"] -= 0.5*V["iemb"]*Tau["eajm"];
	WIAJB["iajb"] -= 0.5*V["iemb"]*T["ej"]*T["am"]; 
	WIAJB["iajb"] -= V["ijmb"]*T["am"]; 
	WIAJB["iajb"] += V["ieab"]*T["ej"]; 
	
	WIABJ["iabj"] = V["ibaj"]; 
	WIABJ["iabj"] -= 0.5*V["ibme"]*T["aemj"];
	WIABJ["iabj"] -= V["ibme"]*T["am"]*T["ej"]; 
	WIABJ["iabj"] += V["ibae"]*T["ej"]; 
	WIABJ["iabj"] -= V["ibmj"]*T["am"]; 
	WIABJ["iabj"] += V["ibme"]*T["eamj"]; 
	WIABJ["iabj"] -= 0.5*V["iemb"]*T["eamj"]; 
	
	WABCI["abci"] = V["acbi"];
	WABCI["abci"] -= V["acmi"]*T["bm"]; 
	WABCI["abci"] -= T["am"]*V["mcbi"]; 
	
	WIAJK["iajk"] = V["ijak"]; 
	WIAJK["iajk"] += V["ieaf"]*Tau["efjk"];
	
	// Amplitudes
	
	CTF::Matrix<> Zai(nvirt, nocc, NS, mp2.dw); 
	CTF::Tensor<> Zabij(4, vvoo, nsns, mp2.dw); 
	CTF::Tensor<> temp(4, vvoo, nsns, mp2.dw);
	
	temp["abij"] = 2.0*T["abij"]; 
	temp["abij"] -= T["abji"]; 
	
	Zai["ai"] = FAB["ae"]*T["ei"]; 
	Zai["ai"] -= FPIJ["mi"]*T["am"]; 
	Zai["ai"] += 2.0*T["em"]*V["meai"];
	Zai["ai"] -= T["em"]*V["aemi"];
	Zai["ai"] += FIA["me"]*temp["eami"];
	Zai["ai"] += V["meaf"]*temp["efmi"];  
	
	Zabij["abij"] = V["aibj"];
	Zabij["abij"] += T["aeij"]*FAB["be"];
	Zabij["abij"] -= T["abim"]*FIJ["mj"];
	Zabij["abij"] += 0.5*V["aebf"]*Tau["efij"]; 
	Zabij["abij"] += 0.5*Tau["abmn"]*WIJKL["mnij"];  
	Zabij["abij"] -= T["aemj"]*WIAJB["mbie"]; 
	Zabij["abij"] -= WIAJB["maie"]*T["ebmj"];
	Zabij["abij"] += temp["eami"]*WIABJ["mbej"]; 
	Zabij["abij"] += T["ei"]*WABCI["abej"];
	Zabij["abij"] -= T["am"]*WIAJK["mbij"]; 
	
	Zabij["abij"] += T["beji"]*FAB["ae"];
	Zabij["abij"] -= T["bajm"]*FIJ["mi"];
	Zabij["abij"] += 0.5*V["beaf"]*Tau["efji"]; 
	Zabij["abij"] += 0.5*Tau["bamn"]*WIJKL["mnji"];  
	Zabij["abij"] -= T["bemi"]*WIAJB["maje"]; 
	Zabij["abij"] -= WIAJB["mbje"]*T["eami"];
	Zabij["abij"] += temp["ebmj"]*WIABJ["maei"]; 
	Zabij["abij"] += T["ej"]*WABCI["baei"];
	Zabij["abij"] -= T["bm"]*WIAJK["maji"]; 
	
	temp["abij"] = 2.0*T["abij"];
	temp["abij"] -= T["baij"];
	
	Zai["ai"] -= V["meni"]*temp["eamn"];
	
	CTF::Matrix<> Dai(nvirt, nocc, NS, mp2.dw); 
	CTF::Tensor<> Dabij(4, vvoo, nsns, mp2.dw);
	
	Dai["ai"]  = V["ii"]; 
	Dai["ai"] -= V["aa"];
	Dabij["abij"] = V["ii"];
	Dabij["abij"] += V["jj"];
	Dabij["abij"] -= V["aa"];
	Dabij["abij"] -= V["bb"]; 
	
	CTF::Function<> fctr(&divide);
	T.ai->contract(1.0, Zai, "ai", Dai, "ai", 0.0, "ai", fctr); 
	T.abij->contract(1.0, Zabij, "abij", Dabij, "abij", 0.0, "abij", fctr);
	  
} 

void CCSD::calculateTriples(Integrals &V, Amplitudes &T) {
	
	int nvirt = N-nocc;
	int shapeNSNS[6] = {NS,NS,NS,NS,NS,NS}; 
	int vvvooo[6] = {nvirt, nvirt, nvirt, nocc, nocc, nocc}; 

	CTF::Tensor<> Dabcijk(6, vvvooo, shapeNSNS, *V.dw);
 
	Dabcijk["abcijk"] += V["ii"];
	Dabcijk["abcijk"] += V["jj"];
	Dabcijk["abcijk"] += V["kk"];
	Dabcijk["abcijk"] -= V["aa"];
	Dabcijk["abcijk"] -= V["bb"];
	Dabcijk["abcijk"] -= V["cc"];

	CTF::Function<> fctr(&divide);
	
	CTF::Tensor<> T3(6, vvvooo, shapeNSNS, mp2.dw, "Tabcijk", 1); 
	
	T3["abcijk"]  = V["beck"]*T["aeij"];
	T3["abcijk"] -= T["abim"]*V["mjck"]; 
	T3["abcijk"] += V["cebj"]*T["aeik"];
	T3["abcijk"] -= T["acim"]*V["mkbj"]; 
	T3["abcijk"] += V["aeck"]*T["beji"];
	T3["abcijk"] -= T["bajm"]*V["mick"]; 
	T3["abcijk"] += V["ceai"]*T["bejk"];
	T3["abcijk"] -= T["bcjm"]*V["mkai"]; 
	T3["abcijk"] += V["beai"]*T["cekj"];
	T3["abcijk"] -= T["cbkm"]*V["mjai"]; 
	T3["abcijk"] += V["aebj"]*T["ceki"];
	T3["abcijk"] -= T["cakm"]*V["mibj"]; 
 
  	{
		CTF::Tensor<> T3bar(6, vvvooo, shapeNSNS, mp2.dw, "Tbarabcijk", 1); 
		T3bar["abcijk"] = (4.0/3.0) * T3["abcijk"];
		T3bar["abcijk"] -= 2.0 * T3["acbijk"]; 
		T3bar["abcijk"] += (2.0/3.0) * T3["bcaijk"];  
		
		T3bar.contract(1.0, T3bar, "abcijk", Dabcijk, "abcijk", 0.0, "abcijk", fctr);
		
		triples_energy = T3bar["abcijk"]*T3["abcijk"];
	}
	
	CTF::Tensor<> Z3(6, vvvooo, shapeNSNS, mp2.dw, "Zabcijk", 1); 
	
	Z3["abcijk"] = 4.0*T["ai"]*V["jbkc"]; 
	Z3["abcijk"] += 4.0*T["bj"]*V["iakc"];
	Z3["abcijk"] += 4.0*T["ck"]*V["iajb"]; 
	Z3["abcijk"] -= 6.0*T["ai"]*V["jckb"];
	Z3["abcijk"] -= 6.0*T["cj"]*V["iakb"];
	Z3["abcijk"] -= 6.0*T["bk"]*V["iajc"]; 
	Z3["abcijk"] += 2.0*T["bi"]*V["jcka"];
	Z3["abcijk"] += 2.0*T["cj"]*V["ibka"];
	Z3["abcijk"] += 2.0*T["ak"]*V["ibjc"]; 
	
	Z3.contract(1.0, Z3, "abcijk", Dabcijk, "abcijk", 0.0, "abcijk", fctr); 

	triples_energy += (1.0 / 3.0) * T3["abcijk"] * Z3["abcijk"]; 
}
					
