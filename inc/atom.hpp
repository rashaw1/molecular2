/*
 *
 *       PURPOSE: To define a class Atom, which will contain the basic info
 *                and methods needed by the Molecule class, for each atom in
 *                a molecular calculation. 
 *    
 *       class Atom:
 *                owns: bfs - a set of contracted basis functions
 *                data: x, y, z - the cartesian coordinates of the atom
 *                      nbfs - the number of basis functions
 *                      nshells - the number of shells of basis functions
 *                      shells - a list corresponding to which bfs are in 
 *                               which shell
 *                      lnums - a list of what angular momentum type each
 *                              shell in shells has
 *                      charge - the atomic number (i.e. the charge in a.u.)
 *                      mass - the atomic mass (i.e. the mass in a.m.u)
 *                accessors: all of the above have get... routines
 *                           note that getCoords() returns an [x, y, z] vector
 *                           bfs has a setBasis(Basis) routine
 *                           getShellBF(shell, i) - return the ith bf from a given shell 
 *                           getShellPrim(shell, i) - same, but for primitives
 *                           getNSpherical() returns the number of cgbfs in spherical basis
 *                           getNSpherShellBF(shell) - returns the number of cgbfs in shell
 *                                     in the spherical basis.
 *                routines:
 *                      rotate(Matrix U) - rotate coords according to unitary
 *                                         transformation matrix, U
 *                      translate(x, y, z) - translate coords [+x, +y, +z]
 *                      
 *       DATE           AUTHOR             CHANGES 
 *       =====================================================================
 *       27/08/15       Robert Shaw        Original code.
 *
 */

#ifndef ATOMHEADERDEF
#define ATOMHEADERDEF

// Includes
#include "eigen_wrapper.hpp"
#include "basis.hpp"

// Declare forward dependencies
class ECPBasis;

// Begin class definition
class Atom
{
private:
  int charge, core;
  double x, y, z, mass;
  double *pos; 
public:
  // Constructors
  Atom() : charge(-1), core(0) { } // Default
  Atom(const Vector& coords, int q, double m); // q = charge, m = mass
  Atom(const Atom& other); // Copy constructor

  // Accessors
  int getCharge() const { return charge; }
  int getEffectiveCharge() const { return charge - core; }
  void setCore(ECPBasis& ecpset); 
  double getMass() const { return mass; }

  Vector getCoords() const; 
  double* getPos() const { return pos; }
  double getX() const { return x; }
  double getY() const { return y; }
  double getZ() const { return z; }

  // Routines
  void rotate(const Matrix& U); 
  void translate(double dx, double dy, double dz);
  // Overloaded operators
  Atom& operator=(const Atom& other);
};  
  
#endif
