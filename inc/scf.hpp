/*
*
*   PURPOSE: To declare a class SCF which is used for doing HF-SCF calculations
*            of both the rhf and uhf kind.
*
*       class SCF
*             owns: focker - a Fock class instance, for doing the bulk of the legwork
*                   molecule - a reference to the molecule in question
*             data: last_dens - the previous density matrix, for convergence checking
*                   last_CP - the previous coefficient matrix, for DIIS
*                   energy - the SCF energy
*             routines: calcE - calculates the energy
*                       rhf - does a restricted HF calculation
*                       uhf - does an unrestricted HF calculation
*
*   DATE             AUTHOR               CHANGES
*   =================================================================================
*   15/09/15         Robert Shaw          Original code.
*
*/

#ifndef SCFHEADERDEF
#define SCFHEADERDEF

#include "fock.hpp"
#include "eigen_wrapper.hpp"
#include "molecule.hpp"
#include "diis.hpp"


// Declare forward dependencies
class IntegralEngine;
class Command; 

// Begin class
class SCF
{
private:
	Command& cmd; 
	Molecule& molecule;
	Fock& focker;
	DIISEngine diis;
	double energy, last_energy, one_E, two_E, error, last_err;
public:
	// Constructor
	SCF(Command& c, Molecule& m, Fock& f);
	// Routines
	void calcE();
	double getEnergy() const { return energy; } 
	double calcE(const Matrix& hcore, const Matrix& dens, const Matrix& fock); 
	Vector calcErr(const Matrix& F, const Matrix& D, const Matrix& S, const Matrix& orthog);
	Vector calcErr();
	bool testConvergence(double val);
	void rhf(bool print = true);
	void uhf(bool print = true);
	void uhf_internal(bool print, UnrestrictedFock& ufocker);
};

#endif
