#ifndef GSHELL_HEAD
#define GSHELL_HEAD

#include <vector>

struct GaussianShell {
	std::vector<double> exps;
	std::vector<double> coeffs;
	double* centerVec;
	int l;
	
	GaussianShell(double* A, int l);
	void addPrim(double exp, double c);
	
	int nprimitive() const { return exps.size(); }
	int ncartesian() const { return ((l+1)*(l+2))/2; }
	double* center() const { return centerVec; };
	double exp(int i) const { return exps[i]; }
	double coef(int i) const { return coeffs[i]; }
	int am() const { return l; }
};

#endif