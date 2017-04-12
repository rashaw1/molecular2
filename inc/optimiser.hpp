#ifndef OPTIMISER_HEADER
#define OPTIMISER_HEADER	


#include "molecule.hpp"
#include "fock.hpp"
#include "eigen_wrapper.hpp"
#include "ProgramController.hpp"
#include <vector>

double rfo(Vector& dx, Vector& grad, Matrix& hessian, double alpha, double stepsize); 
double trust_newton(Vector& dx, Vector& grad, Matrix& hessian, double stepsize); 
std::pair<double, double> rfo_prime(Vector& dy, Vector &grad, Matrix& hessian, double alpha); 

struct RHFCalculator { 
	Command& cmd; 
	SharedMolecule mol; 
	double old_e, delta_e, grad_norm, energy; 
	int iter; 
	
	RHFCalculator(Command& _cmd, SharedMolecule _m) : cmd(_cmd), mol(_m), old_e(0.0), iter(0) { }
	
	double operator()(const Vector& dx, Vector& grad, Matrix& hessian);
};

struct RHFOptimiser {

	Command& cmd; 
	SharedMolecule mol;
	RHFCalculator calc; 
	
	RHFOptimiser(Command& _cmd, SharedMolecule _m) : cmd(_cmd), mol(_m), calc(_cmd, _m) { }
	
	void optimise(); 
	void frequencies(Matrix& hessian); 
	
}; 
#endif