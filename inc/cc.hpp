#ifndef CCHEADERDEF
#define CCHEADERDEF

#include "tensor4.hpp"
#include "fock.hpp"

class IntegralEngine;

class CCSD
{
private:
	int N, nocc;
	MP2& mp2;
	Matrix spinFock;
	Matrix singles;
	Tensor4 doubles;
public:
	CCSD(MP2& _mp2);
	void build_fock();
	void build_guess();
	
};

#endif 
