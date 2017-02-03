/*
 *    PURPOSE: To implement basis.hpp, defining class Basis, representing a basis set.
 *   
 *    DATE            AUTHOR              CHANGES
 *    =====================================================================
 *    27/08/15        Robert Shaw         Original code.
 *
 */

// Includes
#include "basis.hpp"
#include "bf.hpp"
#include "basisreader.hpp"
#include "ioutil.hpp"
#include <iostream>
#include <array>

// Constructors and destructors
Basis::Basis(std::map<int, std::string> ns, iVector& atoms, bool _ecps) : names(ns), maxl(0), ecps(_ecps)
{
  BasisReader input(ns); // Make a basis reading object
  auto it = ns.find(0);
  if (it != ns.end()) name = it->second;
  else name = "sto-3g";
  std::transform(name.begin(), name.end(), name.begin(), ::toupper);
  int natoms = atoms.size(); // Get how many different atoms there are
  // Now determine how many basis functions are needed
  int nbfs = 0;
  iVector qnbfs(natoms); // Store the nbfs for each atom
  for (int i = 0; i < natoms; i++){
    qnbfs[i] = input.readNbfs(atoms(i));
    nbfs += qnbfs(i);
  }
  
  // Allocate memory
  if(nbfs > 0){
    bfs = new BF[nbfs];
  } else {
    bfs = NULL;
  }
  charges.resize(nbfs);  

  // Fill the array, and charges
  int k = 0; // Index counter for bfs array
  for (int i = 0; i < natoms; i++){
    for (int j = 0; j < qnbfs(i); j++){
      bfs[k] = input.readBF(atoms(i), j);
      charges[k] = atoms(i);
      k++;
    }
  }

  // Read in the shell data
  shells = input.readShells(atoms);
  lnums = input.readLnums(atoms);
}

Basis::~Basis()
{
  if (charges.size() > 0){
    delete[] bfs;
  }
}

std::string Basis::getName(int q) const {
	auto it = names.find(q);
	std::string n = name;
	if (it != names.end()) n = it->second;
	return n; 
}

// Accessors
// Find the first position in bfs at which atom of atomic number q occurs
int Basis::findPosition(int q) const
{
  // Loop to find q
  int position = 0;
  bool found = false;
  while (!found && position < charges.size()) {
    found = (charges(position) == q ? true : false);
    position++;
  }
  return position-1;
}

// Same as above but for shells instead of bfs
int Basis::findShellPosition(int q) const
{
  int position = findPosition(q);
 
  int sum = 0;
  int i = 0;
  // Locate the entries in shells corresponding to q
  while(sum < position && i < shells.size()){
    sum+= shells(i);
    i++;
  }
  // i now points to the first element in shells corresponding to q
  return i;
}


// Return the ith basis function corresponding to an atom
// with atomic number q
BF& Basis::getBF(int i)
{
  return bfs[i];
}

BF& Basis::getBF(int q, int i)
{
  int position = findPosition(q);
  // No bounds checking
  return bfs[position+i];
}

// Return the number of basis functions that an atom of atomic number q has
int Basis::getSize(int q) const
{
  int position = findPosition(q);
  // Now find the last position of q
  int size = 1;
  bool found = false;
  while(!found && (size+position) < charges.size()){
    found = (charges(position+size) != q ? true : false);
    if (!found){
      size++;
    }
  }
  return size;
}

// Find the number of shells that atom of atomic number q has
int Basis::getShellSize(int q) const
{
  int nbfs = getSize(q);
  int i = findShellPosition(q);

  // Starting 
  int sum = 0;
  int size = 0;
  while(sum < nbfs && i < shells.size()){
      sum += shells(i);
      i++; size++;
  }
  //if (i < shells.size()){
  //  size--;
  //}
  return size;
}

// Return the subset of shells corresponding to q
iVector Basis::getShells(int q) const
{
  int position = findShellPosition(q);
  int size = getShellSize(q);
  iVector s(size);

  // Fill in the values
  for (int i = 0; i < size; i++){
    s[i] = shells[position+i];
  }
  return s;
}

// Find th subset of lnums corresponding to atom with atomic number q
iVector Basis::getLnums(int q) const
{
  int position = findShellPosition(q);
  int size = getShellSize(q);
  iVector l(size);
  // Fill in the values
  for (int i = 0; i < size; i++){
    l[i] = lnums[position+i];
  }
  return l;
}

void Basis::addShell(int l, std::vector<libint2::real_t> &exps, std::vector<std::vector <libint2::real_t>> &coeffs, double *pos, int atom) {
	using libint2::Shell;
	
	maxl = l > maxl ? l : maxl; 
	
	std::vector<libint2::Shell::Contraction> contr_arr; 
	for (auto &contr : coeffs) {
		libint2::Shell::Contraction newC;
		newC.l = l;
		for (auto c : contr) newC.coeff.push_back(c);
		contr_arr.push_back(newC);
	}
	
	std::array<libint2::real_t, 3> O = { {pos[0], pos[1], pos[2]} };

	intShells.push_back(Shell(exps, contr_arr, O));
	shellAtomList.push_back(atom);
}

  
// Overloaded operators
Basis& Basis::operator=(const Basis& other)
{
  // If basis functions already exist, deallocate memory
  if(charges.size() > 0){
    delete[] bfs;
  }

  // Assign attributes
  name = other.name;
  names = other.names;
  charges = other.charges;
  shells = other.shells;
  lnums = other.lnums;
  ecps = other.ecps;
  
  // Copy across bfs
  int nbfs = charges.size();
  if (nbfs > 0){
    bfs = new BF[nbfs];
    for (int i = 0; i < nbfs; i++){
      bfs[i] = other.bfs[i];
    }
  }
  return *this;
}

  


