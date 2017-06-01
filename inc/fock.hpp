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
#include <string>

// Forward declarations
class Atom;
class Command; 

struct Domain {
	std::vector<int> starts, sizes, centres;
	int totalsize; 
	Matrix G; 
	
	Domain() {}
	Domain(const Domain& other) {
		starts = other.starts;
		sizes = other.sizes;
		centres = other.centres;
		totalsize = other.totalsize; 
	}
	
	void sumsizes() {
		totalsize = 0; 
		for (auto s : sizes)
			totalsize += s;
	}
	
	void print() {
		std::cout << "STARTS:" << std::endl; 
		for (auto s : starts)
			std::cout << s << ", "; 
		std::cout << std::endl;
		
		std::cout << "SIZES:" << std::endl;
		for (auto s : sizes)
			std::cout << s << ", "; 
		std::cout << std::endl;
		
		std::cout << "CENTRES:" << std::endl;
		for (auto c : centres)
			std::cout << c << ", ";
		std::cout << std::endl << std::endl; 
	}
};

struct FragmentInfo {
	int occ, nbfs, naux, start, auxstart, ndfshells; 
	double radius, mo_thresh, fit_thresh, r_thresh; 
	Vector com; 
	
	FragmentInfo() : occ(0), nbfs(0), naux(0), start(0), auxstart(0), radius(0.0) {}
	FragmentInfo(const FragmentInfo& other) {
		occ = other.occ;
		nbfs = other.nbfs;
		naux = other.naux;
		start = other.start;
		auxstart = other.auxstart; 
		radius = other.radius;
		com = other.com; 
	}
};

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
  Matrix fockinc; 
  Matrix CP;
  Matrix forces;
  Matrix hessian;
  Matrix xyK; 
  Matrix Linv, Linv2; 
  Vector eps;
  std::vector<Matrix> focks;
  std::vector<Domain> lmo_domains, ao_domains, fit_domains;
  std::vector<std::vector<int> > Kblocks; 
  Matrix dens, dens_diff;
  IntegralEngine& integrals;
  SharedMolecule molecule;
  bool direct, twoints, fromfile, diis, density_fitted;
  double precision; 
  int nbfs, iter, MAX, nocc;
  std::string guess; 
public:
	
    bool reset_incremental, started_incremental; 
    double next_reset, rms_error, incremental_threshold; 
	int last_reset; 
  
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
  Matrix& getXYK() { return xyK; }
  Matrix& getLinv() { return Linv; }
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
  std::vector<Domain>& getLMODomains() { return lmo_domains; }
  std::vector<Domain>& getAODomains() { return ao_domains; }
  std::vector<Domain>& getFitDomains() { return fit_domains; }
  virtual Matrix getS() const { return integrals.getOverlap(); }
  virtual Matrix getDens() const { return dens; }
  void setDIIS(bool d) { diis = d; } 
  void formHCore();
  void formOrthog();
  virtual void transform(bool first = false);
  virtual void diagonalise();
  virtual void makeJK(Matrix& P, double multiplier = 1.0);
  void formJK(Matrix& P, double multiplier = 1.0);
  void formJKdirect(const Matrix& Schwarz, Matrix& P1, double multiplier = 1.0);
  virtual void formJKfile();
  virtual void formJKdf(Matrix& Cocc, double multiplier = 1.0); 
  virtual void makeFock(); 
  virtual void makeFock(Matrix& P, double multiplier = 1.0 );
  virtual void makeDens();
  virtual void average(Vector &w);
  void simpleAverage(Matrix& D0, double weight = 0.5);
  
  virtual void addDiis(); 
  virtual void clearDiis(); 
  
  virtual void compute_forces(const std::vector<Atom> &atoms, int nocc); 
  virtual void compute_hessian(const std::vector<Atom> &atoms, int nocc);
  virtual void compute_hessian_numeric(const std::vector<Atom> &atoms, int nocc, Command &cmd);
  
  Matrix compute_2body_fock_general(
  	const std::vector<libint2::Shell>& obs, const Matrix& D, const std::vector<libint2::Shell>& D_bs,
  bool D_is_sheldiagonal = false,  // set D_is_shelldiagonal if doing SOAD
  double precision = std::numeric_limits<double>::epsilon()  // discard contributions smaller than this
  		); 

  Matrix compute_2body_fock_df(const Matrix& Cocc);
  Matrix compute_2body_fock_df_local(Matrix& Cocc, const Matrix& sigmainv, Matrix& Pt, std::vector<FragmentInfo>& finfo); 
  Matrix compute_2body_fock_df_local_file(Matrix& Cocc, const Matrix& sigmainv, Matrix& Pt, std::vector<FragmentInfo>& finfo);
  void build_domains(Matrix& Cocc, Matrix& V, std::vector<FragmentInfo>& finfo);  
  
  void compute_soad_guess(); 
  
  virtual Vector compute_xgrad(double fx, Matrix& xhessian, std::vector<int>& activex, Command& cmd); 
  
  template <unsigned deriv_order>
  std::vector<EMatrix> compute_2body_fock_deriv(const std::vector<Atom> &atoms, const EMatrix& D);
};

class UnrestrictedFock : public Fock
{
protected:
	Matrix fock_alpha_ao, fock_beta_ao; 
	Matrix fock_alpha_mo, fock_beta_mo;
	Matrix dens_alpha, dens_beta; 
	Matrix kints_alpha, kints_beta, jints_alpha, jints_beta; 
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
	Matrix& getKAlpha() { return kints_alpha; }
	Matrix& getKBeta() { return kints_beta; }
	Matrix& getJAlpha() { return jints_alpha; }
	Matrix& getJBeta() { return jints_beta; }
	Matrix& getCPAlpha() { return CP_alpha; }
	Matrix& getCPBeta() { return CP_beta; }
	Vector& getEpsAlpha() { return eps_alpha; }
	Vector& getEpsBeta() { return eps_beta; }
	
    virtual void transform(bool first = false);
    virtual void diagonalise();
    virtual void makeJK();
   	void formJK(Matrix& P1, Matrix& P2, double multiplier = 1.0);
    void formJKdirect(const Matrix& Schwarz, Matrix& P1, Matrix& P2, double multiplier = 1.0);
	void formJKdf(Matrix& ca, Matrix& cb, double multiplier = 1.0); 
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
