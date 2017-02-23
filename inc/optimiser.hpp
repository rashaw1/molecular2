#ifndef OPTIMISER_HEADER
#define OPTIMISER_HEADER	


#include "molecule.hpp"
#include "fock.hpp"
#include "eigen_wrapper.hpp"
#include "ProgramController.hpp"
#include <vector>

void quadratic_scf(Command& cmd, SharedMolecule mol);
void conjugate_scf(Command& cmd, SharedMolecule mol);
Vector quadratic(Matrix& hessian, Vector& gradient);

struct ConjugateGradient {
	std::vector<double> grad_norms;
	Vector prev_grad, prev_step; 
	double step_size; 
	bool first; 
	
	ConjugateGradient() : first(true) {}
	Vector compute(Vector& gk);
};

#endif