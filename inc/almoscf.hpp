

#ifndef ALMOSCFHEADERDEF
#define ALMOSCFHEADERDEF

#include "fock.hpp"
#include "matrix.hpp"
#include "molecule.hpp"
#include "diis.hpp"
#include <vector>
#include "integrals.hpp"
#include <Eigen/Dense>

// Declare forward dependencies
class Vector;

// Begin class
class ALMOSCF
{
private:
	Molecule& molecule;
	Fock& focker;
	std::vector<FockFragment> fragments; 
	std::vector<IntegralEngine> ints; 
	DIISEngine diis;
	double dimer_energy, e_frz, e_pol, e_ct, e_int; 
	double delta_e, delta_d;
	std::vector<double> monomer_energies;
	int nfrags; 
	Eigen::MatrixXd P; 
public:
	// Constructor
	ALMOSCF(Molecule& m, Fock& f);
	// Routines
	void calcE();
	double getDimerEnergy() const { return dimer_energy; } 
	std::vector<double>& getMonomerEnergies() { return monomer_energies; }
	
	void rscf();
	void compute();
};

#endif
