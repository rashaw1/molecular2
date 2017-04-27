#ifndef RPAHEADERDEF
#define RPAHEADERDEF

#include <fock.hpp>
#include "ProgramController.hpp"
#include <vector>

class RPA
{
private:
	Command& cmd; 
	Fock& focker; 
	
	int N, nocc, nvirt; 
	double energy; 

public:
	CTF::World dw;
	
	RPA(Command& c, Fock& f); 
	
	void longrange_eris(const std::vector<libint2::Shell>&, std::vector<CTF::Tensor<> >&); 
	void full_eris(std::vector<CTF::Tensor<> >&); 
	
	void compute(); 
	
	double getEnergy() const { return energy; }
};

#endif 
