#ifndef CCHEADERDEF
#define CCHEADERDEF

#include "tensor4.hpp"
#include "matrix.hpp"
#include "mp2.hpp"

class IntegralEngine;

class CCSD
{
private:
	int N, nocc;
	MP2& mp2;
	Matrix spinFock;
	Matrix singles;
	Tensor4 doubles;
	Matrix Dia;
	Tensor4 Dijab;
	double energy, delta_e;
	double delta_singles, delta_doubles;
public:
	CCSD(MP2& _mp2);
	void build_fock();
	void build_guess();
	void build_intermediates(Matrix& F, Tensor4& W, Tensor4& tau, Tensor4& tautilde);
	void build_amplitudes(Matrix& F, Tensor4& W, Tensor4& tau, Tensor4& tautilde);
	void calculateEnergy();
	void compute();
};

#endif 
