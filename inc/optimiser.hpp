#ifndef OPTIMISER_HEADER
#define OPTIMISER_HEADER	


#include "molecule.hpp"
#include "fock.hpp"
#include "eigen_wrapper.hpp"
#include "ProgramController.hpp"

void quadratic_scf(Command& cmd, SharedMolecule mol, Fock& f);
Vector quadratic(Matrix& hessian, Vector& gradient);

#endif