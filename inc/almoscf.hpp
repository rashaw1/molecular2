

#ifndef ALMOSCFHEADERDEF
#define ALMOSCFHEADERDEF

#include "fock.hpp"
#include "matrix.hpp"
#include "molecule.hpp"
#include "diis.hpp"
#include <vector>

// Declare forward dependencies
class IntegralEngine;
class Vector;

// Begin class
class ALMOSCF
{
private:
  Molecule& molecule;
  Fock& focker;
  std::vector<FockFragment> fragments; 
  DIISEngine diis;
  double dimer_energy, e_frz, e_pol, e_ct, e_int; 
  std::vector<double> monomer_energies;
  Matrix Sinv; 
  int nfrags; 
public:
  // Constructor
  ALMOSCF(Molecule& m, Fock& f);
  // Routines
  void calcE();
  double getDimerEnergy() const { return dimer_energy; } 
  std::vector<double>& getMonomerEnergies() { return monomer_energies; }
	
  void rscf();
};

#endif
