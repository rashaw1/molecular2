/*
 *
 *   PURPOSE: To declare class BasisReader, for reading a basis set in from files.
 * 
 *
 *            class BasisReader - reads in basis specifications:
 *                    data: name, ifstream input
 *                    routines:
 *                       readNbfs(q) - read the number of basis functions an atom
 *                                     with charge q has in this basis set
 *                       readBF(q, i) - reads the ith basis function corresponding to
 *                                      atom q.
 *                       readShells(qs) - create a vector of the #shells that the atoms
 *                                        specified in qs would have
 *                       readLnums(qs) - same, but gives the ang. momentum. numbers
 *            
 *
 *   DATE         AUTHOR         CHANGES 
 *   ==========================================================================
 *   30/08/15     Robert Shaw    Original code.
 *
 */

#ifndef BASISREADERHEADERDEF
#define BASISREADERHEADERDEF

// Includes
#include <string>
#include <fstream>
#include <map>
#include "eigen_wrapper.hpp"
#include "gshell.hpp"
#include "pugixml.hpp"

// Declare forward dependencies
class BF;
class Basis;
class ECP;
class ECPBasis;

class BasisReader
{
private:
  std::map<int, std::string> names;
  std::map<std::string, int> obs_list, jk_list, mp2_list, ecp_list; 
  
  std::map<std::string, int> read_basis_list(std::string name); 
  std::vector<GaussianShell> read_basis(std::string atom, double* pos,
  	 									std::string basis, std::string name, 
  										std::map<std::string, int>& basis_list); 

public:
  BasisReader(std::map<int, std::string> ns); // Constructor

  ECP readECP(int q, ECPBasis& ecpset, double *pos); 
  void readShellBasis(Basis& b, int q, double *pos, int atom, int type = 0);
};

#endif
