
#include "mp2.hpp"
#include "logger.hpp"
#include "eigen_wrapper.hpp"
#include "molecule.hpp"
#include "integrals.hpp"
#include "error.hpp"
#include <iostream>
#include <thread>
#include "ProgramController.hpp"

// Constructor
MP2::MP2(Fock& _focker) : spinBasis(false), focker(_focker)
{
	N = focker.getDens().rows();
	nocc = focker.getMolecule()->getNel()/2;
	energy = 0.0;
	moInts.assign(N, 0.0);
}

// Integral transformation
void MP2::transformIntegrals()
{
	if (focker.getMolecule()->control->get_option<bool>("direct")) {
		Error e("MP2TRANS", "Integral direct MP2 not implemented yet.");
		focker.getMolecule()->control->log.error(e);
		nocc = 0;
	} else {
	
		// Multithread
		int nthreads = focker.getMolecule()->control->get_option<int>("nthreads");

		std::vector<Tensor4> moTemps(nthreads);
		std::vector<std::thread> thrds(nthreads);
		std::vector<int> startPoints(nthreads);
		std::vector<int> endPoints(nthreads);

		int spacing = N/nthreads;
		spacing = (double(N)/double(nthreads) - spacing > 0.5 ? spacing+1 : spacing);
		startPoints[0] = 0;
		endPoints[0] = spacing;
	
		for (int i = 0; i < nthreads; i++){
			moTemps[i] = Tensor4(endPoints[i]-startPoints[i], N, N, N, 0.0);
			thrds[i] = std::thread(&MP2::transformThread, *this, startPoints[i], endPoints[i], std::ref(moTemps[i]));
			if (i < nthreads - 1) { 
				startPoints[i+1] = endPoints[i];
				endPoints[i+1] = (i == (nthreads - 2) ? N : endPoints[i] + spacing);
			}
		}
		for (int i = 0; i < nthreads; i++) {
			thrds[i].join();

			// Copy in temporary matrix to moInts
			for (int p = startPoints[i]; p < endPoints[i]; p++){
				for (int q = 0; q < N; q++){
					for (int r = 0; r < N; r++){
						for (int s = 0; s < N; s++){
							moInts(p, q, r, s) = moTemps[i](p-startPoints[i], q, r, s);
						}
					}
				}
			}
		}
		
		focker.getIntegrals().clearTwoInts();
	}
}

void MP2::spatialToSpin() {
	if (!spinBasis) {
		if (focker.getMolecule()->control->get_option<bool>("direct")) {
			Error e("MP2TRANS", "Integral direct CC not implemented yet.");
			focker.getMolecule()->control->log.error(e);
			nocc = 0;
		} else if (moInts.getW() == 0) {
			Error e("SPINTRANS", "Integrals have not yet been transformed to the MO basis");
			focker.getMolecule()->control->log.error(e);
			nocc = 0;
		} else {
			spinInts.assign(2*N, 0.0);
			for (int p = 0; p < N; p++){ 
				for (int q = 0; q < N; q++) {
					for (int r = 0; r < N; r++) {
						for (int s = 0; s < N; s++) {
							auto val1 = moInts(p, r, q, s);
							auto val2 = moInts(p, s, q, r);
							auto diff = val1 - val2; 
							
							int P = 2*p; int Q = 2*q; int R = 2*r; int S = 2*s;
							
							spinInts.set(P, Q, R, S, diff);
							spinInts.set(P, Q+1, R, S+1, val1);
							spinInts.set(P, Q+1, R+1, S, -val2);
							spinInts.set(P+1, Q, R+1, S, val1);
							spinInts.set(P+1, Q, R, S+1, -val2);
							spinInts.set(P+1, Q+1, R+1, S+1, diff); 
							
						}
						
					}
				}
			}
			spinBasis = true;
			moInts.resize(0);
		}
	
	}
}
void MP2::transformThread(int start, int end, Tensor4& moTemp)
{
	Matrix& C = focker.getCP();
	IntegralEngine& aoInts = focker.getIntegrals();

	int offset = end - start;
	Tensor4 temp1(offset, N, N, N, 0.0);
	Tensor4 temp2(offset, N, N, N, 0.0);
	Tensor4 temp3(offset, N, N, N, 0.0);

	// Transform as four 'quarter-transforms'
	for (int p = start; p < end; p++){
		
		for (int mu = 0; mu < N; mu++){
			for (int a = 0; a < N; a++){
				for (int b = 0; b < N; b++){
					for (int c = 0; c < N; c++)
						temp1(p-start, a, b, c) += C(mu, p)*aoInts.getERI(mu, a, b, c);
				} // b
			} // a
		} // mu
		
		for (int q = 0; q < N; q++){
			for (int nu = 0; nu < N; nu++){
				for (int b = 0; b < N; b++){
					for (int c = 0; c < N; c++)
						temp2(p-start, q, b, c) += C(nu, q)*temp1(p-start, nu, b, c);
				} // b
			} // nu

			for (int r = 0; r < N; r++){
				for (int lam = 0; lam < N; lam++){
					for (int c = 0; c < N; c++)
						temp3(p-start, q, r, c) += C(lam, r)*temp2(p-start, q, lam, c);
				} // lam
				
				for (int s = 0; s < N; s++){
					for (int sig = 0; sig < N; sig++)
						moTemp(p-start, q, r, s) += C(sig, s)*temp3(p-start, q, r, sig);
				} // s
			} // r
		} // q
	} // p
}

// Determine the MP2 energy
void MP2::calculateEnergy()
{
	Vector& eps = focker.getEps();

	double ediff, etemp;
	energy = 0.0;
	int occ_min = 0; 
	for (int i = 0; i < nocc; i++)
		if (eps[i] >= -1.0) { occ_min = i; break; } 
	
	for (int i = occ_min; i < nocc; i++){
		for (int j = occ_min; j < nocc; j++){
			
			for (int a = nocc; a < N; a++){	
				for (int b = nocc; b < N; b++){
					ediff = eps[i] + eps[j] - eps[a] - eps[b];
					etemp = moInts(i, a, j, b)*(2.0*moInts(a, i, b, j) - moInts(a, j, b, i));
					
					energy += etemp/ediff;
				} // b
			} // a
		} // j
	} // i
}

void MP2::calculateEnergy(const S4OddTensor4& t) {
	energy = 0.0;
	
	for (int i = 0; i < 2*nocc; i++)
		for (int j = 0; j < 2*nocc; j++)
			for (int a = 2*nocc; a < 2*N; a++)
				for (int b = 2*nocc; b < 2*N; b++) energy += spinInts(i, j, a, b) * t(i, j, a-2*nocc, b-2*nocc);
	
	energy *= 0.25;
}	
					
