/*
 *
 *     PURPOSE: To define the class Basis representing a
 *              basis set.
 *
 *     class Basis:
 *              owns: bfs - a set of BFs
 *              data: name - the name of the basis set, needed for file io
 *                    charges - a list of atomic numbers corresponding to
 *                            each BF in bfs.
 *                    shells - a list representing which bfs are in which
 *                             shells - i.e. if the first 3 are in one shell,
 *                             the next two in another, etc, it would be
 *                             3, 2, ...
 *                    lnums - a list with the angular quantum number 
 *                            of each shell in shells                 
 *              routines:
 *                    findPosition(q) - find the first position in the bfs 
 *                                      array of an atom of atomic number q
 *                    findShellPosition(q) - find the first position in the 
 *                                           shells array
 *                    getNBFs() - returns the total number of basis functions
 *                    getName() - returns the name
 *                    getCharges() - returns the vector of charges
 *                    getBF(charge, i) - return the ith BF corresponding to
 *                                       an atom of atomic number charge   
 *                    getSize(charge) - return how many basis functions an atom
 *                                      of atomic number charge has
 *                    getShellSize(charge) - return how many shells it has
 *                    getShells(charge) - return a vector of the correct subset
 *                                        of shells
 *                    getLnums(charge) - return vector of subset of lnums
 *                 
 *
 *     DATE        AUTHOR            CHANGES
 *     ====================================================================
 *     27/08/15    Robert Shaw       Original code.
 *
 */

#ifndef BASISHEADERDEF
#define BASISHEADERDEF

// Includes
#include "eigen_wrapper.hpp"
#include <string>
#include <libint2.hpp>
#include <vector>
#include <map>

// Forward declarations
class BF;

// Begin class definitions

class Basis
{
private:
  BF* bfs;
  std::string name;
  std::map<int, std::string> names;
  iVector charges;
  iVector shells;
  iVector lnums;
  int maxl;
  bool ecps; 
  std::vector<libint2::Shell> intShells;
  std::vector<int> shellAtomList;
public:
  // Constructors and destructor
  // Note - no copy constructor, as it doesn't really seem necessary
  // Need to specify the name of the basis, n, and a list of the 
  // distinct atoms that are needed (as a vector of atomic numbers)
  Basis() : name("Undefined") { } // Default constructor
  Basis(std::map<int, std::string> ns, iVector& atoms, bool _ecps = false);
  ~Basis(); // Destructor
  // Accessors
  int getNBFs() const { return charges.size(); }
  std::string getName() const { return name; }
  std::string getName(int q) const;
  iVector getCharges() const { return charges; }
  int findPosition(int q) const;
  int findShellPosition(int q) const;
  bool hasECPS() const { return ecps; }
  BF& getBF(int i);
  BF& getBF(int q, int i);
  int getSize(int q) const;
  int getMaxL() const { return maxl; }
  int getShellSize(int q) const;
  iVector getShells(int q) const;
  int getShellAtom(int i) const { return shellAtomList[i]; }
  std::vector<libint2::Shell>& getIntShells() { return intShells; }
  iVector getLnums(int q) const;
  
  void addShell(int l, std::vector<libint2::real_t> &exps, std::vector<std::vector <libint2::real_t>> &coeffs, double *pos, int atom);
  
  // Overloaded operators
  Basis& operator=(const Basis& other);
};

#endif
