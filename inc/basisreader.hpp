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

// Declare forward dependencies
class BF;
class Basis;
class ECP;
class ECPBasis;

class BasisReader
{
private:
  std::map<int, std::string> names;
  std::ifstream input;
  void openFile(int q);
  void closeFile();
public:
  BasisReader(std::map<int, std::string> ns) : names(ns) {} // Constructor
  int readNbfs(int q);
  BF readBF(int q, int i);
  ECP readECP(int q, ECPBasis& ecpset, double *pos);
  iVector readShells(int q);
  iVector readLnums(int q);
  iVector readShells(iVector& qs);
  iVector readLnums(iVector& qs);
 
  void readShellBasis(Basis& b, int q, double *pos, int atom, int type = 0);
};

#endif
