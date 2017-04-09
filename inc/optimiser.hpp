#ifndef OPTIMISER_HEADER
#define OPTIMISER_HEADER	


#include "molecule.hpp"
#include "fock.hpp"
#include "eigen_wrapper.hpp"
#include "ProgramController.hpp"
#include <vector>

void optimise(Command& cmd, SharedMolecule m); 

#endif