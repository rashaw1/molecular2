
#include "cc.hpp"
#include "logger.hpp"
#include "mvector.hpp"
#include "matrix.hpp"
#include "molecule.hpp"
#include "integrals.hpp"
#include "error.hpp"
#include <iostream>
#include <thread>

// Constructor
CCSD::CCSD(MP2& _mp2) : mp2(_mp2) {
	N = mp2.getN();
	nocc = mp2.getNocc();
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
					t(i, j, a - nocc, b-nocc) = spinInts(i, j, a, b) / resolvent;
				}
			
		}
	}
	
	mp2.calculateEnergy(doubles);	
}


					
