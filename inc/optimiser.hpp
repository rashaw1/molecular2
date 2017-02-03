#ifndef OPTIMISER_HEADER
#define OPTIMISER_HEADER	


#include "molecule.hpp"
#include "fock.hpp"
#include "eigen_wrapper.hpp"

void quadratic_scf(Molecule& mol, Fock& f);
Vector quadratic(Matrix& hessian, Vector& gradient);

#endif