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
	Fock& focker;
public:
	MP2(Fock& _focker);
	void transformIntegrals();
	void spatialToSpint();
	void transformThread(int start, int end, Tensor4& moTemp);
	void calculateEnergy();
	double getEnergy() const { return energy; }
	double getEnergy(const Tensor4& amplitudes);
	Tensor4& getInts() { return moInts; }
	int getN() const { return N; }
	int getNocc() const { return nocc; }
	Fock& getFock() { return focker; }
};

#endif 
