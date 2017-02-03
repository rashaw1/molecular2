/*
*
*     Implements molecule.hpp, defining class Molecule
*
*     DATE         AUTHOR           CHANGES  
*     ==============================================================
*     26/08/15     Robert Shaw      Original code.
*
*/

#include <iostream>
#include "logger.hpp"
#include "molecule.hpp"
#include "error.hpp"
#include "basisreader.hpp"
#include <cmath>
#include <vector>
#include <libint2.hpp>

// An initialisation function, for code reuse purposes
void Molecule::init()
{
	// Set the data variables
	charge = log.getCharge();
	multiplicity = log.getMultiplicity();
	parent = true;

	// Get the basis set     
	bfset = log.getBasis();

	// Now declare the array of atoms
	natoms = log.getNatoms();
	if (natoms != 0){
    
		atoms = new Atom[natoms];
		nel = 0; // Initialise to no electrons
	
		// Populate the array, and give the atoms
		// their basis functions. At the same time,
		// calculate the number of electrons
		for (int i = 0; i < natoms; i++) {
			atoms[i] = log.getAtom(i);
			atoms[i].setBasis(bfset);
			nel += atoms[i].getCharge();
		}

		// Account for overall charge
		nel -= charge;

	} else { // Nothing to do
		Error e("INIT", "There are no atoms!");
		log.error(e);
		atoms = NULL;
	}
}

void Fragment::init(Atom* as, int nat, int q, int mult)
{
	charge = q;
	parent = false;
	multiplicity = mult;
	atoms = as; 
	natoms = nat;
	if (nat <= 0) {
		Error e("FRAGINIT", "There are no atoms in this fragment!");
		log.error(e);
	} else {
		bfset = log.getBasis();
		
		nel = 0;
		for (int i = 0; i < natoms; i++) 
			nel += atoms[i].getCharge();
	}
}

// Constructor
Molecule::Molecule(Logger& logger, int q, bool doInit) : log(logger)
{
	// Set the log and then initialise
	if(doInit) init();
}

// Copy constructor
Molecule::Molecule(const Molecule& other) : log(other.log)
{
	init();
}


// Destructor
Molecule::~Molecule()
{
	if(log.getNatoms()!=0 && parent){
		delete [] atoms;
	}
}

Fragment::Fragment(Logger& logger, Atom* as, int nat, int q, int mult) : Molecule(logger, q, false)  {
	init(as, nat, q, mult); 
}

Fragment::Fragment(const Fragment& other) : Molecule(other.log, other.charge, false) {
	init(other.atoms, other.natoms, other.charge, other.multiplicity);
}
Fragment::~Fragment() { }


// Routines

// rotate(Matrix U) rotates the coordinate system according
// to the transformation specified by the unitary matrix U
// U needs to be 3x3
void Molecule::rotate(const Matrix& U)
{
	// Check if 3 x 3
	if(U.rows() == 3 && U.cols() == 3){
		for (int i = 0; i < log.getNatoms(); i++){
			atoms[i].rotate(U); // Do the rotation
		}
		updateBasisPositions();
	} else {
		Error e("ROTATE", "Unsuitable rotation matrix.");
		log.error(e);
	}
}

// translate(x, y, z) translates the coordinate system by
// adding the vector (x, y, z) to the coordinate vector
// of each atom.
void Molecule::translate(double x, double y, double z)
{
	for (int i = 0; i < log.getNatoms(); i++){
		atoms[i].translate(x, y, z);
	}
	updateBasisPositions();
}

// nalpha and nbeta return the number of alpha and beta
// (spin up/spin down) electrons
int Molecule::nalpha() const
{
	// multiplicity = 2s+1; each unpaired alpha electron contributes
	// s = 1/2, therefore s = nunpairedalpha*(1/2), so
  
	return (nel-multiplicity+1)/2 + multiplicity -1;
}

int Molecule::nbeta() const
{
	// The number of beta electrons is (nel - nunpairedalpha)/2
	return (nel - multiplicity + 1)/2;
}

// Calculate the nuclear energy, i.e. the sum of terms
// (Zi * Zj)/Rij for each distinct pair of nucleii, i, j
void Molecule::calcEnuc()
{
	double zi;
	enuc = 0.0;

	// Outer loop over all atoms
	for (int i = 0; i < natoms; i++) {
		zi = atoms[i].getEffectiveCharge();
		// Inner loop over all atoms with index
		// greater than i, to avoid double counting
		for (int j = i+1; j < natoms; j++){
			enuc += (zi*atoms[j].getEffectiveCharge())/dist(i, j);
		}
	}
}
  
// Compute the centre of mass of the molecule
Vector Molecule::com() const 
{
	// Coordinate vector of centre of mass, set to all
	// zeroes for summations
	Vector c = Vector::Zero(3);
  
	// Loop over all atoms
	double sum = 0.0;
	double m;
	for (int i = 0; i < log.getNatoms(); i++){
		m = atoms[i].getMass();
		c = c + m*atoms[i].getCoords();
		sum += m;
	}

	c = (1.0/sum)*c;
	return c;
}

// Calculate the principal moments of inertia;
// optionally, shift to the inertial coordinate system
// First the elements of the inertia tensor are
// calculated, and then this is diagonalised to 
// give the eigenvalues - the principal moments - 
// which are returned as a 3-vector
Vector Molecule::getInertia(bool shift) 
{
	Vector v(3); // To return the moments in
	Matrix inertia = Matrix::Zero(3, 3); // To contain the elements

	// We must translate to centre of mass coordinates
	Vector c = -1.0*com();
	translate(c(0), c(1), c(2));

	// Loop over all atoms
	double m;
	for (int i = 0; i < log.getNatoms(); i++){
		m = atoms[i].getMass();
		c = atoms[i].getCoords();

		// The diagonal elements are the mass times the sum
		// of the squares of the other coordinates
		// e.g. I(0, 0) = m(y^2 + z^2)
		inertia(0, 0) += m*(c(1)*c(1) + c(2)*c(2));
		inertia(1, 1) += m*(c(0)*c(0) + c(2)*c(2));
		inertia(2, 2) += m*(c(0)*c(0) + c(1)*c(1));

		// The off diagonal elements are minus the 
		// mass times the two coordinates
		// e.g. I(1, 0) = -m*x*y = I(0 , 1)
		inertia(0, 1) -= m*c(0)*c(1);
		inertia(0, 2) -= m*c(0)*c(2);
		inertia(1, 2) -= m*c(1)*c(2);
	}
	// Finally fill out the symmetric elements
	inertia(1, 0) = inertia(0, 1); 
	inertia(2, 0) = inertia(0, 2); 
	inertia(2, 1) = inertia(1, 2);
    
	// Now we diagonalise, only need eigenvalues
	EigenSolver es(inertia);
	v = es.eigenvalues();
	if(shift)
		rotate(es.eigenvectors());
	
	return v;
}


// Calculate the (Euclidean) distance between the atoms i and j
// dist2 gives the square of this distance
double Molecule::dist2(int i, int j) const
{
	Vector dvec = atoms[j].getCoords() - atoms[i].getCoords(); 
	return dvec.dot(dvec);
}

double Molecule::dist(int i, int j) const {
	return sqrt(dist2(i, j));
}

// Calculate bond, out of plane, and torsional angles
double Molecule::bondAngle(int i, int j, int k) const
{
	// Get the differences of the position vectors
	// and calculate the angle between them
	Vector rij = atoms[j].getCoords();
	Vector rjk = rij;
	rij = rij - atoms[i].getCoords();
	rjk = atoms[k].getCoords() - rjk;
	return angle(rij, rjk); // Vector class friend function
}

double Molecule::oopAngle(int i, int j, int k, int l) const
{
	// The out of plane angle is given by the formula:
	// [(elj x elk).eli]/sin(anglejlk)
	// and sin(anglejlk) = ||ejl x elk||
	// therefore start by forming the unit vectors
	Vector co = atoms[l].getCoords();
	Vector elk = atoms[k].getCoords() - co;
	Vector elj = atoms[j].getCoords() - co;
	Vector eli = atoms[i].getCoords() - co;
  
  	co = cross(elj, elk);
	double sinangle = co.norm();
  
	// Normalise
	elk.normalize();
	elj.normalize();
	eli.normalize();
  
	// Do the triple product and divide by the sine
	return (cross(elj, elk)).dot(eli)/sinangle;
}

double Molecule::torsionAngle(int i, int j, int k, int l) const
{
	// The torsional angle is given by the formula:
	// [(eij x ejk).(ejk x ekl)]/[sin(ijk)*sin(jkl)]

	// Form the identity vectors
	Vector jv = atoms[j].getCoords();
	Vector kv = atoms[k].getCoords();
	Vector eij = jv - atoms[i].getCoords();
	Vector ejk = kv - jv;
	Vector ekl = atoms[l].getCoords() - kv;
  
	// Normalise
	eij.normalize();
	ejk.normalize();
	ekl.normalize();

	// Get angles and cross products
	double sinijk = (cross(eij, ejk)).norm();
	double sinjkl = (cross(ejk, ekl)).norm();
	jv = cross(eij, ejk);
	kv = cross(ejk, ekl);
 
	// Return the inner product over the product of angles
	return jv.dot(kv)/(sinijk*sinjkl);
}

// Compute the rotational type and constants of the molecule
// using the getInertia method
std::string Molecule::rType()
{
	std::string rstring = "asymmetric"; // Asymm. is default case
	// Diatomic case
	if(log.getNatoms() == 2){
		rstring = "diatomic";
	} else {
		// Calculate the principal moments of inertia
		Vector I = getInertia(false);

		// All moments are >= 0, so no need to use fabs
		// Linear if Ic = Ib > Ia = 0
		// Prolate symm. if Ic = Ib > Ia /= 0
		// Oblate symm. if Ic > Ib = Ia /= 0
		// Spherical if Ic = Ib = Ia /= 0
		double CUTOFF = 0.01;

		if ( (I(2)-I(1)) < CUTOFF && I(0) < CUTOFF ){
			rstring = "linear";
		}   else if ( (I(2) - I(1)) < CUTOFF && (I(1)-I(0)) > CUTOFF &&  I(0) >= CUTOFF){
			rstring = "prolate";
		} else if ( (I(1) - I(0)) < CUTOFF && (I(2) - I(1)) > CUTOFF &&  I(0) >= CUTOFF){
			rstring = "oblate";
		} else if ( (I(1) - I(0)) < CUTOFF && (I(2) - I(1)) < CUTOFF){
			rstring = "spherical";
		}
	}
	return rstring;
}

Vector Molecule::rConsts(int units)
{
	// First get the principal moments of inertia
	Vector I = getInertia(false);
  
	double K;
	// Choose units
	if (units == 0){ // reciprocal centimetres
		K = Logger::RTOCM; 
	} else { // GHz
		K = Logger::RTOGHZ; // ~8.39214994 x 10^2 GHz . A^2 . amu
	}
	// Constant Bi = K/Ii 
	for (int i = 0; i < 3; i++){
		if (I(i) < log.precision()) { I[i] = 0.0; }  
		else { I[i] = K/I(i); }
	}
	return I;
}

void Molecule::buildShellBasis() {
	BasisReader b(log.bnames);
	for (int i = 0; i < natoms; i++) {
		Atom &a = atoms[i];
		double pos[3] = { a.getX(), a.getY(), a.getZ() };
		b.readShellBasis(bfset, a.getCharge(), pos, i);
	}
}

void Molecule::buildECPBasis() {
	if (bfset.hasECPS()) {
		BasisReader b(log.bnames); 
		for (int i = 0; i < natoms; i++) {
			Atom& a = atoms[i];
			auto it = log.bnames.find(-a.getCharge());
			if (it != log.bnames.end()) {
				ECP newECP = b.readECP(a.getCharge(), ecpset, a.getPos());
				ecpset.addECP(newECP, i);
			}
			atoms[i].setCore(ecpset);
			nel -= atoms[i].getCharge() - atoms[i].getEffectiveCharge();
		}
	}
}

void Molecule::updateBasisPositions() {
	std::vector<libint2::Shell>& shells = bfset.getIntShells();
	for (int i = 0; i < natoms; i++) {
		auto x = atoms[i].getX(); auto y = atoms[i].getY(); auto z = atoms[i].getZ();
		for (int j = 0; j < shells.size(); j++) {
			if (bfset.getShellAtom(j) == i) {
				shells[j].O[0] = x;
				shells[j].O[1] = y;
				shells[j].O[2] = z;
			}
		}
		
		if (bfset.hasECPS()) {
			for (int j = 0; j < ecpset.getN(); j++) {
				if (ecpset.getAtom(j) == i) ecpset.getECP(j).setPos(x, y, z);
			}
		}
	}
}
  


  
  
