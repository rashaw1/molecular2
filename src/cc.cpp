
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

void CCSD::calculateTriples(Integrals &V, Amplitudes &T) {
	
	int nvirt = N-nocc;
	int shapeASAS[6] = {AS,NS,NS,AS,NS,NS};
	
	int vvvooo[6] = {nvirt, nvirt, nvirt, nocc, nocc, nocc}; 
	CTF::Tensor<> T3(6, vvvooo, shapeASAS, mp2.dw, "Tabcijk", 1); 
	CTF::Tensor<> Z3(6, vvvooo, shapeASAS, mp2.dw, "Tabcijk", 1); 
	
	Z3["abcijk"]  = V["bcek"]*T["aeij"];
	Z3["abcijk"] -= V["amij"]*T["bcmk"];

	T3 = Z3;
	
	int sh_sym[6] = {SH, NS, NS, SH, NS, NS};
	CTF::Tensor<> Dabcijk(6, T3.lens, sh_sym, *V.dw);
 
	Dabcijk["abcijk"] += V["i"];
	Dabcijk["abcijk"] += V["j"];
	Dabcijk["abcijk"] += V["k"];
	Dabcijk["abcijk"] -= V["a"];
	Dabcijk["abcijk"] -= V["b"];
	Dabcijk["abcijk"] -= V["c"];

	CTF::Function<> fctr(&divide);

	T3.contract(1.0, Z3, "abcijk", Dabcijk, "abcijk", 0.0, "abcijk", fctr);
	
	Z3["abcijk"] += V["abij"]*T["ck"]; 
	
	triples_energy = (1.0/36.0) * T3["efgmno"] * Z3["efgmno"]; 
}

void CCSD::compute()
{
	Logger& log = mp2.getFock().getMolecule()->control->log;
	
	log.title("CCSD CALCULATION");
	
	mp2.calculateEnergy();
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
	while (!converged && iter < MAXITER) {
		ccsd_iteration(V, T, 0); 
		double singles_en = V["ai"]*T["ai"];
		double disconn_en = 0.5*V["abij"]*T["ai"]*T["bj"]; 
		double doubles_en = 0.25*V["abij"]*T["abij"]; 
		new_en = singles_en + disconn_en + doubles_en; 
		std::cout << singles_en << " " << disconn_en << " " << doubles_en << std::endl;
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

	int ns[2] = {NS, NS}; 
	int nsns[4] = {NS, NS, NS, NS}; 
	int asns[4] = {AS, NS, NS, NS}; 
	int nsas[4] = {NS, NS, AS, NS}; 
	 
	CTF::Tensor<> Tau = CTF::Tensor<>(T.abij); 
	CTF::Tensor<> TauTilde = CTF::Tensor<>(T.abij); 
	CTF::Tensor<> FAE(2, V.ab->lens, ns, *V.dw); 
	CTF::Tensor<> FMI(2, V.ij->lens, ns, *V.dw);  
	CTF::Tensor<> FME = CTF::Tensor<>(V.ia);
	CTF::Tensor<> WMNIJ(4, V.ijkl->lens, asns, *V.dw);
	CTF::Tensor<> WABEF(4, V.abcd->lens, nsas, *V.dw);
	CTF::Tensor<> WBMEJ = CTF::Tensor<>(V.aibj);
	
	Tau["abij"] = T["ai"]*T["bj"]; 
	Tau["abij"] -= T["bi"]*T["aj"]; 
	TauTilde["abij"] = T["abij"];
	TauTilde["abij"] += 0.5*Tau["abij"];
	Tau["abij"] += T["abij"]; 
	
	// Intermediates
	FME["me"] = V["me"]; 
	FME["me"] += V["mnef"]*T["fn"]; 
	
	FMI["mi"] = V["mi"]; 
	FMI["mi"] += 0.5*V["me"]*T["ei"]; 
	FMI["mi"] -= V["mnei"]*T["en"]; 
	FMI["mi"] += 0.5*V["mnef"]*TauTilde["efin"]; 
 	
	FAE["ae"] = V["ae"]; 
	FAE["ae"] -= 0.5*V["me"]*T["am"];
	FAE["ae"] += V["amef"]*T["fm"]; 
	FAE["ae"] -= 0.5*V["mnef"]*TauTilde["afmn"]; 
	
	WMNIJ["mnij"] = V["mnij"]; 
	WMNIJ["mnij"] -= V["mnei"]*T["ej"]; 
	WMNIJ["mnij"] += V["mnej"]*T["ei"];
	WMNIJ["mnij"] += 0.25*V["mnef"]*Tau["efij"]; 
	
	WABEF["abef"] = V["abef"];
	WABEF["abef"] -= V["amef"]*T["bm"];
	WABEF["abef"] += V["bmef"]*T["am"]; 
	WABEF["abef"] += 0.25*V["mnef"]*Tau["abmn"]; 
	
	WBMEJ["bmej"] = V["bmej"]; 
	WBMEJ["bmej"] += V["bmef"]*T["fj"];
	WBMEJ["bmej"] += V["mnej"]*T["bn"];
	WBMEJ["bmej"] -= 0.5*V["mnef"]*T["bfjn"]; 
	WBMEJ["bmej"] += V["mnef"]*T["fj"]*T["bn"]; 
		
	// Iteration
	
	CTF::Tensor<> Zai = CTF::Tensor<>(V.ai);
	CTF::Tensor<> Zabij = CTF::Tensor<>(V.abij); 
	
	Zai["ai"] = V["ai"];
	Zai["ai"] += FAE["ae"]*T["ei"]; 
	Zai["ai"] -= FMI["mi"]*T["am"];
	Zai["ai"] += FME["me"]*T["aeim"];
	Zai["ai"] -= V["anfi"]*T["fn"];
	Zai["ai"] += 0.5*V["amef"]*T["efim"]; 
	Zai["ai"] -= 0.5*V["nmei"]*T["aemn"]; 
	
	
	CTF::Tensor<> tempAE(2, V.ab->lens, ns, *V.dw);
	CTF::Tensor<> tempMJ(2, V.ij->lens, ns, *V.dw);
	Zabij["abij"] = 4.0 * V["abij"];
	tempAE["be"] = FAE["be"] - 0.5*FME["me"]*T["bm"]; 
	Zabij["abij"] += 2.0 * tempAE["be"] * T["aeij"]; 
	tempAE["ae"] = FAE["ae"] - 0.5*FME["me"]*T["am"];
	Zabij["abij"] -= 2.0* tempAE["ae"] * T["beij"];
	tempMJ["mj"] = FMI["mj"] + 0.5*FME["me"]*T["ej"]; 
	Zabij["abij"] -= 2.0 * tempMJ["mj"] * T["abim"]; 
	tempMJ["mi"] = FMI["mi"] + 0.5*FME["me"]*T["ei"];
	Zabij["abij"] += 2.0 * tempMJ["mi"] * T["abjm"]; 
	Zabij["abij"] += WMNIJ["mnij"]*Tau["abmn"]; 
	Zabij["abij"] += WABEF["abef"]*Tau["efij"]; 
	Zabij["abij"] += 2.0 * V["abej"] * T["ei"];
	Zabij["abij"] -= 2.0 * V["abei"] * T["ej"]; 
	Zabij["abij"] += 2.0 * V["bmij"] * T["am"];
	Zabij["abij"] -= 2.0 * V["amij"] * T["bm"]; 
	Zabij["abij"] -= WBMEJ["bmej"] * T["aeim"]; 
	Zabij["abij"] += V["bmej"]*T["ei"]*T["am"]; 
	Zabij["abij"] += WBMEJ["bmei"] * T["aejm"];
	Zabij["abij"] -= V["bmei"]*T["ej"]*T["am"]; 
	Zabij["abij"] += WBMEJ["amej"]*T["beim"]; 
	Zabij["abij"] -= V["amej"]*T["bm"];
	Zabij["abij"] -= WBMEJ["amei"]*T["bejm"]; 
	Zabij["abij"] += V["amei"]*T["ej"]*T["bm"];
	
	CTF::Tensor<> Dai(2, V.ai->lens, V.ai->sym, *V.dw);
	int sh_sym[4] = {SH, NS, SH, NS};
	CTF::Tensor<> Dabij(4, V.abij->lens, sh_sym, *V.dw);
	Dai["ai"] += V["i"];
	Dai["ai"] -= V["a"];
 
	Dabij["abij"] += V["i"];
	Dabij["abij"] += V["j"];
	Dabij["abij"] -= V["a"];
	Dabij["abij"] -= V["b"];


	CTF::Function<> fctr(&divide);

	T.ai->contract(1.0, Zai, "ai", Dai, "ai", 0.0, "ai", fctr);
	T.abij->contract(0.125, Zabij, "abij", Dabij, "abij", 0.0, "abij", fctr);
} 

					
