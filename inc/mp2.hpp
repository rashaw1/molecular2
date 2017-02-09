#ifndef MP2HEADERDEF
#define MP2HEADERDEF

#include "tensor4.hpp"
#include "fock.hpp"
#include "ProgramController.hpp"

class IntegralEngine;

class MP2
{
private:
	Command& cmd;
	int N, nocc;
	double energy;
	S8EvenTensor4 moInts;
	S8OddTensor4 spinInts;
	bool spinBasis;
	Fock& focker;
public:
	MP2(Fock& _focker);
	void transformIntegrals();
	void spatialToSpin();
	void transformThread(int start, int end, Tensor4& moTemp);
	void calculateEnergy();
	void calculateEnergy(const S4OddTensor4& amplitudes);
	double getEnergy() const { return energy; }
	S8OddTensor4& getSpinInts() { 
		return spinInts;
	}
	S8EvenTensor4& getMOInts() {
		return moInts;
	}
	int getN() const { return N; }
	int getNocc() const { return nocc; }
	Fock& getFock() { return focker; }
};

#endif 
