#ifndef CCHEADERDEF
#define CCHEADERDEF

#include "tensor4.hpp"
#include "matrix.hpp"
#include "mp2.hpp"
#include "diis.hpp"
#include <vector>

class IntegralEngine;

class CCSD
{
private:
	int N, nocc, iter;
	MP2& mp2;
	Matrix spinFock;
	Matrix singles;
	S4OddTensor4 doubles;
	Matrix Dia;
	
	DIISEngine diis;
	bool doDiis;
	int maxDiis;
	std::vector<Matrix> singles_cache;
	std::vector<S4OddTensor4> doubles_cache;

	S4EvenTensor4 Dijab;
	double energy, delta_e, triples_energy;
	double delta_singles, delta_doubles;
	bool withTriples;
public:
	CCSD(MP2& _mp2, bool _triples = false, bool _doDIIS = true);
	void build_fock();
	void build_guess();
	void build_intermediates(Matrix& F, Tensor4& W, S4OddTensor4& tau, S4OddTensor4& tautilde);
	void build_amplitudes(Matrix& F, Tensor4& W, S4OddTensor4& tau, S4OddTensor4& tautilde);
	void calculateEnergy();
	void calculateTriples();
	void calculateError(Matrix& newSingles, S4OddTensor4& newDoubles);
	void compute();
};

#endif 
