

#ifndef ALMOSCFHEADERDEF
#define ALMOSCFHEADERDEF

#include "fock.hpp"
#include "eigen_wrapper.hpp"
#include "molecule.hpp"
#include "diis.hpp"
#include <vector>
#include "integrals.hpp"
#include "ProgramController.hpp"

// Begin class
class ALMOSCF
{
private:
	SharedMolecule molecule;
	Command& cmd; 
	Fock& focker;
	std::vector<FockFragment> fragments;
	std::vector<UnrestrictedFockFragment> ufragments; 
	std::vector<IntegralEngine> ints; 
	DIISEngine diis;
	double dimer_energy, e_frz, e_pol, e_ct, e_int, e_disp, e_pert_2, e_pert_4, e_mon_rpa; 
	double delta_e, delta_d;
	std::vector<double> monomer_energies;
	int nfrags, MAX; 
	Matrix P, P_alpha, P_beta, sigma, sigma_alpha, sigma_beta; 
public:
	// Constructor
	ALMOSCF(Command& c, SharedMolecule m, Fock& f);
	// Routines
	void rperturb(bool order4 = false);
	void uperturb(bool order4 = false);
	void setFragments(bool unrestricted = false);
	double getDimerEnergy() const { return dimer_energy; } 
	std::vector<double>& getMonomerEnergies() { return monomer_energies; }
	
	void rscf();
	void uscf();
	void rcompute();
	void ucompute();
	double makeDens(bool alpha);
};

#endif
