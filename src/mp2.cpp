
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
void MP2::transformIntegrals(bool withSpin)
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
	
	if (withSpin) {
		S8OddTensor4 tempSpinInts(2*N, 0.0); 
	
		for (int p = 0; p < N; p++){ 
			for (int q = 0; q < N; q++) {
				for (int r = 0; r < N; r++) {
					for (int s = 0; s < N; s++) {
						auto val1 = moInts(p, r, q, s);
						auto val2 = moInts(p, s, q, r);
						auto diff = val1 - val2; 
							
						int P = 2*p; int Q = 2*q; int R = 2*r; int S = 2*s;
							
						tempSpinInts.set(P, Q, R, S, diff);
						tempSpinInts.set(P, Q+1, R, S+1, val1);
						tempSpinInts.set(P, Q+1, R+1, S, -val2);
						tempSpinInts.set(P+1, Q, R+1, S, val1);
						tempSpinInts.set(P+1, Q, R, S+1, -val2);
						tempSpinInts.set(P+1, Q+1, R+1, S+1, diff); 
							
					}
						
				}
			}
		}
		spinBasis = true;
		moInts.resize(0);
	
		N *= 2;
		nocc *= 2;
		int nvirt = N-nocc; 
		int offset = 0;
	
		if(!focker.getMolecule()->control->get_option<bool>("withcore")) {
			offset = nocc - focker.getMolecule()->getNValence(); 
			nocc -= offset; 
			N = nocc + nvirt; 
		}
	
		int offnocc = nocc + offset; 
		int offN = N + offset; 
	
		focker.getMolecule()->control->log.print("No. of occ. orbitals: " + std::to_string(nocc));
		focker.getMolecule()->control->log.print("No. of virt. orbitals: " + std::to_string(nvirt)); 
	
		spinInts = std::make_shared<Integrals>(nocc, nvirt,  dw);
		Integrals& V = *spinInts; 
	
		int64_t sz, *indices;
		double *values;
		int ctr;
	
		V.abcd->read_local(&sz, &indices, &values);
		ctr = 0;
		for (int a = offnocc; a < offN; a++) 
			for (int b = offnocc; b < a; b++)
				for (int c = offnocc; c < offN; c++)
		for (int d = offnocc; d < c; d++) {
			values[ctr] = tempSpinInts(d, c, b, a);
			ctr++; 
		}				
		V.abcd->write(sz, indices, values);

		V.abci->read_local(&sz, &indices, &values);
		ctr = 0;
		for (int i = offset; i < offnocc; i++) 
			for (int c = offnocc; c < offN; c++)
				for (int b = offnocc; b < offN; b++)
		for (int a = offnocc; a < b; a++) {
			values[ctr] = tempSpinInts(a, b, c, i);
			ctr++; 
		}
		V.abci->write(sz, indices, values);

		V.aibc->read_local(&sz, &indices, &values);
		ctr = 0;
		for (int c = offnocc; c < offN; c++) 
			for (int b = offnocc; b < c; b++)
				for (int i = offset; i < offnocc; i++)
		for (int a = offnocc; a < offN; a++) {
			values[ctr] = tempSpinInts(a, i, b, c);
			ctr++; 
		}
		V.aibc->write(sz, indices, values);

		V.aibj->read_local(&sz, &indices, &values);
		ctr = 0;
		for (int j = offset; j < offnocc; j++) 
			for (int b = offnocc; b < offN; b++)
				for (int i = offset; i < offnocc; i++)
		for (int a = offnocc; a < offN; a++) {
			values[ctr] = tempSpinInts(a, i, b, j);
			ctr++; 
		}
		V.aibj->write(sz, indices, values);

		V.abij->read_local(&sz, &indices, &values);
		ctr = 0;
		for (int j = offset; j < offnocc; j++) 
			for (int i = offset; i < j; i++)
				for (int b = offnocc; b < offN; b++)
		for (int a = offnocc; a < b; a++) {
			values[ctr] = tempSpinInts(a, b, i, j);
			ctr++; 
		}
		V.abij->write(sz, indices, values);
		
		V.ijab->read_local(&sz, &indices, &values);
		ctr = 0;
		for (int b = offnocc; b < offN; b++) 
			for (int a = offnocc; a < b; a++)
				for (int j = offset; j < offnocc; j++)
		for (int i = offset; i < j; i++) {
			values[ctr] = tempSpinInts(i, j, a, b);
			ctr++; 
		}
		V.ijab->write(sz, indices, values);	

		V.aijk->read_local(&sz, &indices, &values);
		ctr = 0;
		for (int k = offset; k < offnocc; k++) 
			for (int j = offset; j < k; j++)
				for (int i = offset; i < offnocc; i++)
		for (int a = offnocc; a < offN; a++) {
			values[ctr] = tempSpinInts(a, i, j, k);
			ctr++; 
		}
		V.aijk->write(sz, indices, values);

		V.ijak->read_local(&sz, &indices, &values);
		ctr = 0;
		for (int k = offset; k < offnocc; k++) 
			for (int a = offnocc; a < offN; a++)
				for (int j = offset; j < offnocc; j++)
		for (int i = offset; i < j; i++) {
			values[ctr] = tempSpinInts(i, j, a, k);
			ctr++; 
		}
		V.ijak->write(sz, indices, values);

		V.ijkl->read_local(&sz, &indices, &values);
		ctr = 0;
		for (int l = offset; l < offnocc; l++) 
			for (int k = offset; k < l; k++)
				for (int j = offset; j < offnocc; j++)
		for (int i = offset; i < j; i++) {
			values[ctr] = tempSpinInts(i, j, k, l);
			ctr++; 
		}
		V.ijkl->write(sz, indices, values);
	
		// Transform core Hamiltonian
		Matrix hcore_spatial = focker.getCP().transpose() * focker.getHCore() * focker.getCP();
		Matrix hcore = Matrix::Zero(N, N);
		for (int p = 0; p < N; p++)
		for (int q = p; q < N; q++) {
			hcore(p, q) = hcore_spatial((p+offset)/2, (q+offset)/2) * (p%2 == q%2);
			hcore(q, p) = hcore(p, q);
		}
	
		// Assign the spin-Fock matrix
		Matrix F = Matrix::Zero(N, N);
		// Build 
		for (int p = 0; p < N; p++) {
			for (int q = p; q < N; q++) {
				F(p, q) = hcore(p, q);
				for (int m = 0; m < offnocc; m++) F(p, q) += tempSpinInts(p+offset, m, q+offset, m);
				F(q, p) = F(p, q);
			}
		}
	
		V.aa->read_local(&sz, &indices, &values);
		for (int a = 0; a < nvirt; a++) values[a] = F(a+nocc, a+nocc); 
		V.aa->write(sz, indices, values); 

		V.ii->read_local(&sz, &indices, &values); 
		for (int i = 0; i < nocc; i++) values[i] = F(i, i); 
		V.ii->write(sz, indices, values); 

		V.ab->read_local(&sz, &indices, &values);
		ctr = 0;  
		for (int a = 0; a < nvirt; a++) { 
			for (int b = 0; b < a; b++) {
				values[ctr] = F(b+nocc, a+nocc); 
				ctr++; 
			}
		}
		V.ab->write(sz, indices, values); 

		V.ai->read_local(&sz, &indices, &values);
		ctr = 0;  
		for (int i = 0; i < nocc; i++) { 
			for (int a = 0; a < nvirt; a++) {
				values[ctr] = F(a+nocc, i); 
				ctr++; 
			}
		}
		V.ai->write(sz, indices, values); 

		V.ia->read_local(&sz, &indices, &values);
		ctr = 0;  
		for (int a = 0; a < nvirt; a++) { 
			for (int i = 0; i < nocc; i++) {
				values[ctr] = F(i, a+nocc);
				ctr++; 
			}
		}
		V.ia->write(sz, indices, values); 

		V.ij->read_local(&sz, &indices, &values);
		ctr = 0;  
		for (int i = 0; i < nocc; i++) { 
			for (int j = 0; j < i; j++) {
				values[ctr] = F(j, i); 
				ctr++; 
			}
		}
		V.ij->write(sz, indices, values);
	
		Vector& eps = focker.getEps();
		amplitudes = std::make_shared<Amplitudes>(nocc, nvirt, dw); 
		Amplitudes& T = *amplitudes;
	 
		T.ai->read_local(&sz, &indices, &values); 
		for (int i = 0; i < sz; i++) values[i] = 0.0;
		T.ai->write(sz, indices, values); 

		T.abij->read_local(&sz, &indices, &values); 
		ctr = 0;
		for (int j = offset; j < offnocc; j++){
			for (int i = offset; i < j; i++){
				auto eocc = eps[i/2] + eps[j/2];
			
				for (int b = offnocc; b < offN; b++){ 
					for (int a = offnocc; a < b; a++) {
						auto resolvent = eocc - eps[a/2] - eps[b/2];
						values[ctr++] = tempSpinInts(i, j, a, b) / resolvent;
					}
				}
			}
		}
		T.abij->write(sz, indices, values);
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
	Amplitudes& T = *amplitudes;
	Integrals& V = *spinInts; 
	energy = 0.25 * V["abij"] * T["abij"]; 
}

					
