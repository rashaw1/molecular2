/*
 *      PURPOSE: To define a class Molecule, which will be the work-horse of
 *               the entire suite of programs. It is defined as follows:
 *        
 *       class Molecule:
 *           owns:
 *                atoms: an array of the atoms contained in the molecule
 *                logger: an object for logging all and any information,
 *                        errors, global data, etc.
 *           data:
 *                charge: the overall charge of the molecule
 *                nel: the total number of electrons in the molecule
 *                multiplicity: the spin multiplicity of the molecule
 *                enuc: the nuclear energy of the molecule
 *                natoms: the number of atoms in atoms
 *           accessors: - all data has a get... accessor, e.g. getNel
 *           routines:
 *                rotate(Matrix U): rotate the coordinate system, according
 *                                  to the transformation in the unitary
 *                                  matrix, U
 *                translate(x, y, z): translate the coordinate system 
 *                nalpha: returns the number of alpha electrons
 *                nbeta: returns the number of beta electrons
 *                calcEnuc: calculates the nuclear energy, and stores in enuc
 *                com: calculates the centre of mass of the molecule
 *                getInertia(shift): calculates the principal moments 
 *                                  of inertia, optionally shifting to the
 *                                  inertial frame of reference if shift=true
 *                dist(i, j): calculates the distance between atoms i and j
 *                dist2(i, j): same as above, but squared
 *                bondAngle(i, j, k): calc. bond angle between atoms i,j,k
 *                oopAngle(i, j, k, l): calc out of plane angle
 *                torsionAngle(i, j, k, l): torsional angle
 *                rType: returns the molecules rotational type, giving one of
 *                       ['diatomic', 'linear', 'asymmetric', 'oblate',
 *                        'prolate', 'spherical']
 *                rConsts(units): returns rotational constants, A,B,C
 *                                units = 0 -> cm-1; 1 -> MHz
 *
 *     DATE            AUTHOR                CHANGES
 *     =======================================================================
 *     26/08/15        Robert Shaw           Original code.
 * 
 */

#ifndef MOLECULEHEADERDEF
#define MOLECULEHEADERDEF

// Includes
#include "basis.hpp"
#include "atom.hpp"
#include <string>
#include <map>
#include <vector>
#include "bf.hpp"
#include "ecp.hpp"
#include "eigen_wrapper.hpp"

// Declare forward dependcies
class Fragment;
class Molecule;
class ProgramController;
struct Construct;

using SharedMolecule = std::shared_ptr<Molecule>; 
using SharedFragment = std::shared_ptr<Fragment>;  

// Begin class definition
class Molecule : public std::enable_shared_from_this<Molecule> 
{
protected:
  ECPBasis ecpset;
  Basis bfset;
  Atom* atoms;
  int charge, nel, multiplicity, natoms;
  bool parent, angstrom, fragmented, has_ecps;
  std::vector<Fragment> fragments;
  std::map<int, std::string> bnames;
  double enuc;
public:
	
  std::shared_ptr<ProgramController> control;
  
  // Constructors and destructor
  void init(Construct& c); // An initialisation function
  Molecule(std::shared_ptr<ProgramController> control, Construct& c, int q = 0); // Need the log for input, q is charge
  Molecule(std::shared_ptr<ProgramController> control, int q); 
  Molecule(const Molecule& other); // Copy constructor
  ~Molecule(); // Deletes the atom array
  
  Atom parseGeomLine(std::string line);
  
  void buildShellBasis();
  void buildECPBasis();
  void updateBasisPositions();
  
  // Accessors
  int getNAtoms() const { return natoms; }
  int getCharge() const { return charge; }
  int getNel() const { return nel; }
  int getMultiplicity() const { return multiplicity; }
  double getEnuc() const { return enuc; }
  Atom& getAtom(int i) { return atoms[i]; } // Return atom i
  BF& getBF(int q, int i) { return bfset.getBF(q, i); } // Return basis func. i of atom q
  Basis& getBasis() { return bfset; }
  ECPBasis& getECPBasis() { return ecpset; }
  ECP& getECP(int i) { return ecpset.getECP(i); }
  bool hasECPS() { return has_ecps; }
  std::vector<Fragment>& getFragments() { return fragments; }
  
  // Routines
  void rotate(const Matrix& U);
  void translate(double x, double y, double z);
  int nalpha() const;
  int nbeta() const;
  void calcEnuc(); 
  Vector com() const;
  Vector getInertia(bool shift = false);
  double dist(int i, int j) const;
  double dist2(int i, int j) const;
  double bondAngle(int i, int j, int k) const;
  double oopAngle(int i, int j, int k, int l) const;
  double torsionAngle(int i, int j, int k, int l) const;
  std::string rType();
  Vector rConsts(int units);
  
  Molecule& operator=(const Molecule& other); 
};

class Fragment : public Molecule 
{
private:
	std::vector<Atom> frag_atoms; 
	SharedMolecule mol; 
public:
	Fragment(std::shared_ptr<ProgramController> control, SharedMolecule m, Atom* as, int nat, int q = 0, int mult = 1); 
	Fragment(const Fragment& other);
	~Fragment();
	void init(Atom* as, int nat, int q, int mult);
	
	Fragment& operator=(const Fragment& other); 
};

#endif
