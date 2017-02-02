/*
 *
 *   PURPOSE:  To define a class Fock, which contains the data and routines
 *             needed for methods based in Fock space. 
 *
 *       class Fock
 *              owns: hcore,  JK, orthog - the core hamiltonian, coulomb/ exchange, orthogonalising matrices
 *                          matrices.
 *                    dens, dens_1, dens_2 - the current and two previous density
 *                          matrices (previous needed for DIIS).
 *                    integrals - the integral engine
 *              data:
 * 
 *              routines:
 *   
 * 
 *     DATE        AUTHOR              CHANGES
 *    ==========================================================================
 *    14/09/15     Robert Shaw         Original code.
 * 
 */

#ifndef FOCKHEADERDEF
#define FOCKHEADERDEF

// Includes
#include "matrix.hpp"
#include "integrals.hpp"
#include "molecule.hpp"
#include "mvector.hpp"
#include <vector>

// Forward declarations
class Atom;

// Begin class definition
class Fock
{
protected:
  Matrix hcore;
  Matrix jkints;
  Matrix jints;
  Matrix kints;
  Matrix orthog;
  Matrix fockm;
  Matrix focka;
  Matrix CP;
  Matrix forces;
  Matrix hessian;
  Vector eps;
  std::vector<Matrix> focks;
  Matrix dens;
  IntegralEngine& integrals;
  Molecule& molecule;
  bool direct, twoints, fromfile, diis;
  int nbfs, iter, MAX;
public:
  Fock(IntegralEngine& ints, Molecule& m);
  Fock(const Fock& other);
  IntegralEngine& getIntegrals() { return integrals; }
  Molecule& getMolecule() { return molecule; }
  Matrix& getHCore() { return hcore; }
  Matrix& getFockAO() { return focka; }
  Matrix& getFockMO() { return fockm; }
  Matrix& getOrthog() { return orthog; }
  Matrix& getCP() { return CP; }
  Vector& getEps() { return eps; }
  Matrix& getJK() { return jkints; }
  Matrix& getJ() { return jints; } 
  Matrix& getK() { return kints; }
  virtual Matrix getS() { return integrals.getOverlap(); }
  Matrix& getDens() { return dens; }
  void setDIIS(bool d) { diis = d; } 
  void formHCore();
  void formOrthog();
  void transform(bool first = false);
  void diagonalise();
  virtual void makeJK();
  virtual void formJK(Matrix& P);
  void formJKdirect(const Matrix& Schwarz, Matrix& P);
  void formJKfile();
  void makeFock();
  virtual void makeFock(Matrix& jbints);
  void makeDens(int nocc);
  void average(Vector &w);
  void simpleAverage(Matrix& D0, double weight = 0.5);
  
  void compute_forces(const std::vector<Atom> &atoms, int nocc); 
  void compute_hessian(const std::vector<Atom> &atoms, int nocc);
  
  template <unsigned deriv_order>
  std::vector<EMatrix> compute_2body_fock_deriv(const std::vector<Atom> &atoms, const EMatrix& D);
};

class FockFragment : public Fock 
{
private:
	int start, end;
public:
	Eigen::MatrixXd Sxx; 
	FockFragment(IntegralEngine& ints, Molecule& m, int start, int end);
	FockFragment(const FockFragment& other);
	void buildFock(Eigen::MatrixXd& qfq, Eigen::MatrixXd& qfp, Eigen::MatrixXd& pfp); 
};

#endif
