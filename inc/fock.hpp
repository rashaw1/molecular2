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
#include "integrals.hpp"
#include "molecule.hpp"
#include "eigen_wrapper.hpp"
#include <vector>

// Forward declarations
class Atom;
class Command; 

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
  SharedMolecule molecule;
  bool direct, twoints, fromfile, diis;
  double precision; 
  int nbfs, iter, MAX, nocc;
public:
  Fock(Command& cmd, IntegralEngine& ints, SharedMolecule m);
  Fock(const Fock& other);
  IntegralEngine& getIntegrals() { return integrals; }
  SharedMolecule getMolecule() { return molecule; }
  Matrix& getHCore() { return hcore; }
  Matrix& getFockAO() { return focka; }
  Matrix& getFockMO() { return fockm; }
  Matrix& getOrthog() { return orthog; }
  Matrix& getCP() { return CP; }
  Vector& getEps() { return eps; }
  Matrix& getJK() { return jkints; }
  Matrix& getJ() { return jints; } 
  Matrix& getK() { return kints; }
  Matrix& getDens() { return dens; }
  Matrix& getForces() { return forces; }
  Matrix& getHessian() { return hessian; }
  virtual Matrix& getS() { return integrals.getOverlap(); }
  Matrix getHCore() const { return hcore; }
  Matrix getFockAO() const { return focka; }
  Matrix getFockMO() const { return fockm; }
  Matrix getOrthog() const { return orthog; }
  Matrix getCP() const { return CP; }
  Vector getEps() const { return eps; }
  Matrix getJK() const { return jkints; }
  Matrix getJ() const { return jints; } 
  Matrix getK() const { return kints; }
  Matrix getForces() const { return forces; }
  Matrix getHessian() const { return hessian; }
  virtual Matrix getS() const { return integrals.getOverlap(); }
  virtual Matrix getDens() const { return dens; }
  void setDIIS(bool d) { diis = d; } 
  void formHCore();
  void formOrthog();
  virtual void transform(bool first = false);
  virtual void diagonalise();
  virtual void makeJK();
  void formJK(Matrix& P, double multiplier = 1.0);
  void formJKdirect(const Matrix& Schwarz, Matrix& P1, double multiplier = 1.0);
  virtual void formJKfile();
  virtual void makeFock();
  virtual void makeDens();
  virtual void average(Vector &w);
  void simpleAverage(Matrix& D0, double weight = 0.5);
  
  virtual void clearDiis(); 
  
  virtual void compute_forces(const std::vector<Atom> &atoms, int nocc); 
  virtual void compute_hessian(const std::vector<Atom> &atoms, int nocc);
  
  template <unsigned deriv_order>
  std::vector<EMatrix> compute_2body_fock_deriv(const std::vector<Atom> &atoms, const EMatrix& D);
};

class UnrestrictedFock : public Fock
{
protected:
	Matrix fock_alpha_ao, fock_beta_ao; 
	Matrix fock_alpha_mo, fock_beta_mo;
	Matrix dens_alpha, dens_beta; 
	Matrix jkints_alpha, jkints_beta; 
	Matrix CP_alpha, CP_beta;
	Vector eps_alpha, eps_beta; 
	std::vector<Matrix> alpha_focks, beta_focks;
	int nalpha, nbeta; 
public:
    UnrestrictedFock(Command& cmd, IntegralEngine& ints, SharedMolecule m);
    UnrestrictedFock(const UnrestrictedFock& other);
	
	Matrix& getFockAlphaAO() { return fock_alpha_ao; }
	Matrix& getFockBetaAO() { return fock_beta_ao; }
	Matrix& getFockAlphaMO() { return fock_alpha_ao; }
	Matrix& getFockBetaMO() { return fock_beta_ao; }
	Matrix& getDensAlpha() { return dens_alpha; }
	Matrix& getDensBeta() { return dens_beta; }
	Matrix& getJKAlpha() { return jkints_alpha; }
	Matrix& getJKBeta() { return jkints_beta; }
	Matrix& getCPAlpha() { return CP_alpha; }
	Matrix& getCPBeta() { return CP_beta; }
	Vector& getEpsAlpha() { return eps_alpha; }
	Vector& getEpsBeta() { return eps_beta; }
	
    virtual void transform(bool first = false);
    virtual void diagonalise();
    virtual void makeJK();
   	void formJK(Matrix& P1, Matrix& P2, double multiplier = 1.0);
    void formJKdirect(const Matrix& Schwarz, Matrix& P1, Matrix& P2, double multiplier = 1.0);
    virtual void formJKfile();
    virtual void makeFock();
    virtual void makeDens();
    virtual void average(Vector &w);
	
	virtual void clearDiis(); 
};

class FockFragment : public Fock 
{
private:
	int start, end;
public:
	Matrix Sxx; 
	FockFragment(Command& cmd, IntegralEngine& ints, SharedMolecule m, int start, int end);
	FockFragment(const FockFragment& other);
	Vector buildFock(Matrix& qfq, Matrix& qfp, Matrix& pfp, bool alpha = false); 
	virtual void gensolve();
};


class UnrestrictedFockFragment : public UnrestrictedFock
{
private: 
	int start, end;
public:
	Matrix Sxx;
	UnrestrictedFockFragment(Command& cmd, IntegralEngine& ints, SharedMolecule m, int start, int end);
	UnrestrictedFockFragment(const UnrestrictedFockFragment& other);
	Vector buildFock(Matrix& qfq, Matrix& qfp, Matrix& pfp, bool alpha); 
	virtual void gensolve();
};

#endif
