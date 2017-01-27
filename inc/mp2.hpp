#ifndef MP2HEADERDEF
#define MP2HEADERDEF

#include "tensor4.hpp"
#include "fock.hpp"

class IntegralEngine;

class MP2
{
private:
	int N, nocc;
	double energy;
	Tensor4 moInts;
	bool spinBasis;
	Fock& focker;
public:
	MP2(Fock& _focker);
	void transformIntegrals();
	void spatialToSpin();
	void transformThread(int start, int end, Tensor4& moTemp);
	void calculateEnergy();
	void calculateEnergy(const Tensor4& amplitudes);
	double getEnergy() const { return energy; }
	Tensor4& getInts() { return moInts; }
	int getN() const { return N; }
	int getNocc() const { return nocc; }
	Fock& getFock() { return focker; }
};

#endif 
