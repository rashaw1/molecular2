#ifndef RPAHEADERDEF
#define RPAHEADERDEF

#include <fock.hpp>
#include "ProgramController.hpp"
#include <vector>
#include "eigen_wrapper.hpp"
#include <libint2.hpp>

struct fInfo {

	std::vector<int> nocc, nvirt;
	Matrix F, S, T, V; 
	std::vector<libint2::Shell>& shells;
	std::vector<libint2::Shell>& df_shells; 
	
	double eintra, edisp, edispexch, eionic, ebsse;  
	
	fInfo(std::vector<libint2::Shell>& s1, std::vector<libint2::Shell>& s2) : shells(s1), df_shells(s2) {} 

};


class RPA
{
private:
	Command& cmd; 
	Fock& focker; 
	
	int N, nocc, nvirt; 
	double energy; 

public:
	CTF::World dw;
	
	RPA(Command& c, Fock& f, int _N, int _nocc);
	
	void eris(const std::vector<libint2::Shell>&, std::vector<CTF::Tensor<> >&, Matrix& T, Matrix& V, bool longrange = false); 
	void df_eris(const std::vector<libint2::Shell>& obs, const std::vector<libint2::Shell>& auxbs, Matrix& T, Matrix& V, Matrix& BmnP); 
	
	void compute(bool print = true); 
	void fcompute(fInfo& info, bool print = false); 
	
	double getEnergy() const { return energy; }
};

#endif 
