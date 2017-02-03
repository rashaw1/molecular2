

#ifndef ALMOSCFHEADERDEF
#define ALMOSCFHEADERDEF

#include "fock.hpp"
#include "eigen_wrapper.hpp"
#include "molecule.hpp"
#include "diis.hpp"
#include <vector>
#include "integrals.hpp"

// Begin class
class ALMOSCF
{
private:
	Molecule& molecule;
	Fock& focker;
	std::vector<FockFragment> fragments;
	std::vector<IntegralEngine> ints; 
	DIISEngine diis;
	double dimer_energy, e_frz, e_pol, e_ct, e_int, e_pert_2, e_pert_4; 
	double delta_e, delta_d;
	std::vector<double> monomer_energies;
	int nfrags, MAX; 
	Matrix P; 
public:
	// Constructor
	ALMOSCF(Molecule& m, Fock& f);
	// Routines
	void perturb(bool order4 = false);
	void setFragments(bool unrestricted = false);
	double getDimerEnergy() const { return dimer_energy; } 
	std::vector<double>& getMonomerEnergies() { return monomer_energies; }
	
	void rscf();
	void rcompute();
};

#endif
