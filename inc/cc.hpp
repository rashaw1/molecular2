#ifndef CCHEADERDEF
#define CCHEADERDEF

#include "tensor4.hpp"
#include <ctf.hpp>
#include "eigen_wrapper.hpp"
#include "mp2.hpp"
#include "diis.hpp"
#include <vector>
#include "ProgramController.hpp"

class IntegralEngine;

class CCSD
{
private:
	int N, nocc, iter;
	Command& cmd;
	MP2& mp2;
	
	DIISEngine diis;
	bool doDiis;
	int maxDiis;
	std::vector<Matrix> singles_cache;
	std::vector<S4OddTensor4> doubles_cache;

	double energy, delta_e, triples_energy;
	double delta_singles, delta_doubles;
	bool withTriples;
public:
	CCSD(Command& c, MP2& _mp2);
	void calculateTriples(Integrals& V, Amplitudes& T);
	void calculateError(Matrix& newSingles, S4OddTensor4& newDoubles);
	void compute();
  	double getEnergy() const { return energy; }
  	double getETriples() const { return triples_energy; }
	
	void ccsd_iteration(Integrals& V, Amplitudes& aT, int sched_nparts = 0);
};

#endif 
