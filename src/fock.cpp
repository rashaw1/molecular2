/*
 *
 *   PURPOSE: To implement class Fock, a class containing the data and routines
 *            needed for Fock space based methods.
 *
 *   DATE          AUTHOR              CHANGES
 *  ===========================================================================
 *  14/09/15       Robert Shaw         Original code.
 * 
 */

//#ifndef EIGEN_USE_MKL_ALL
//#define EIGEN_USE_MKL_ALL
//#endif

#include "fock.hpp"
#include "error.hpp"
#include "tensor4.hpp"
#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <libint2.hpp>
#include "logger.hpp"

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
	EMatrix;

// Constructor
Fock::Fock(IntegralEngine& ints, Molecule& m) : integrals(ints), molecule(m)
{
  Eigen::setNbThreads(m.getLog().getNThreads());
	
  // Make the core hamiltonian matrix
  formHCore();

  // Form the orthogonalising matrix and get initial density guess
  try {
    formOrthog();
  } catch (Error e) {
    molecule.getLog().error(e);
  }

  // Retrieve whether this calculation should be done direct, or whether
  // the twoints matrix has been formed, or if the 2e integrals need to be
  // read from file.
  direct = molecule.getLog().direct();
  diis = molecule.getLog().diis();
  iter = 0;
  MAX = 10;
  twoints = false;
  if (!direct){
    Vector ests = integrals.getEstimates();
    if (ests[3] < molecule.getLog().getMemory())
      twoints = true;
  }

  fromfile = false;
  if (!twoints && !direct)
    fromfile = true;

}

// Form the core hamiltonian matrix
void Fock::formHCore()
{
  // The core Hamiltonian matrix is defined to be
  // the sum of the kinetic and nuclear attraction
  // matrices
  hcore = integrals.getKinetic() + integrals.getNucAttract();
  nbfs = hcore.nrows();
}

void Fock::formOrthog()
{
  Matrix S = integrals.getOverlap();
  Eigen::MatrixXd temp(S.nrows(), S.nrows());
  for (int i = 0; i < S.nrows(); i++){
    for (int j = 0; j < S.nrows(); j++){
      temp(i, j) = S(i, j);
    }
  }
  Matrix U(S.nrows(), S.nrows(), 0.0); Vector lambda(S.nrows(), 0.0);
  // Diagonalise the overlap matrix into lambda and U,
  // so that U(T)SU = lambda
  // We can now form S^(-1/2) - the orthogonalising matrix
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(temp);
  temp = es.eigenvectors();
  for (int i = 0; i < S.nrows(); i++){
    for (int j = 0; j < S.nrows(); j++){
      U(i, j) = temp(i, j);
    }
  }
  temp = es.eigenvalues().asDiagonal();
  for (int i = 0; i < S.nrows(); i++)
    lambda[i] = temp(i, i);
  
  orthog.assign(nbfs, nbfs, 0.0);
  for (int i = 0; i < nbfs; i++) {
    orthog(i, i) = 1.0/(std::sqrt(lambda(i)));
  }
  
  // S^-1/2  = U(lambda^-1/2)U(T)
  orthog = U * orthog * U.transpose();
}
  


// Transform the AO fock matrix to the MO basis 
void Fock::transform(bool first)
{
  if (first) { 
    // Form the core Fock matrix as (S^-1/2)(T)H(S^-1/2)
    fockm = (orthog.transpose()) * ( hcore * orthog);
  } else {
    // Form the orthogonalised fock matrix
    if (diis) DIIS();
   fockm = orthog.transpose() * (focka * orthog);
  }
}

// Diagonalise the MO fock matrix to get CP and eps
void Fock::diagonalise() 
{
  Eigen::MatrixXd temp(fockm.nrows(), fockm.nrows());
  for (int i = 0; i < fockm.nrows(); i++){
    for (int j = 0; j < fockm.nrows(); j++){
      temp(i, j) = fockm(i, j);
    }
  }

  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(temp);
  temp = es.eigenvectors();
  CP.assign(fockm.nrows(), fockm.nrows(), 0.0); eps.assign(fockm.nrows(), 0.0);
  for(int i = 0; i < fockm.nrows(); i++){
    for (int j = 0; j < fockm.nrows(); j++){
      CP(i, j) = temp(i, j);
    }
  }
  
  temp = es.eigenvalues().asDiagonal();
  for (int i = 0; i < fockm.nrows(); i++)
    eps[i] = temp(i, i);

  // Sort the eigenvalues and eigenvectors
  // using a selection sort
  int k;
  for (int i = 0; i < nbfs; i++){
    k=i;
    // Find the smallest element
    for (int j = i+1; j < nbfs; j++)
      if (eps(j) < eps(k)) { k=j; }
    
    // Swap rows of eps and columns of CP
    eps.swap(i, k);
    CP.swapCols(i, k);
  }

  // Form the initial SCF eigenvector matrix
  CP = orthog*CP;
}

// Construct the density matrix from CP, 
// for nocc number of occupied orbitals
void Fock::makeDens(int nocc)
{
  // Form the density matrix
  dens.assign(nbfs, nbfs, 0.0);
  for (int u = 0; u < nbfs; u++){
    for (int v = 0; v < nbfs; v++){
      for (int t = 0; t < nocc; t++){
	dens(u, v) += CP(u, t)*CP(v,t);
      }
    }
  }
  dens = 2.0*dens;
}

// Make the JK matrix, depending on how two electron integrals are stored/needed
void Fock::makeJK()
{
  if (twoints){
    formJK(); 
  } else if (direct) {
    formJKdirect();
  } else {
    try {
      formJKfile();
    } catch (Error e) {
      molecule.getLog().error(e);
    }
  }
}

// Form the 2J-K matrix, given that twoints is stored in memory
void Fock::formJK()
{
  jints.assign(nbfs, nbfs, 0.0);
  kints.assign(nbfs, nbfs, 0.0);
  for (int u = 0; u < nbfs; u++){
    for (int v = 0; v < nbfs ; v++){
      for (int s = 0; s < nbfs; s++){
	for (int l = 0; l < nbfs; l++){
	  jints(u, v) += dens(s, l)*integrals.getERI(u, v, s, l);
	  kints(u, v) += dens(s, l)*integrals.getERI(u, l, s, v);
	}
      }
    }
  }
  jkints = jints - 0.5*kints;
}
// Form JK using integral direct methods
void Fock::formJKdirect()
{
    using libint2::Shell;
    using libint2::Engine;
    using libint2::Operator;
	
	std::vector<Shell>& shells = molecule.getBasis().getIntShells();
	const auto n = integrals.nbasis(shells);
	jints.assign(n, n, 0.0);
	kints.assign(n, n, 0.0);
	const Matrix& D = dens; 
	
	Engine engine(Operator::coulomb, integrals.max_nprim(shells), integrals.max_l(shells), 0);
	auto shell2bf = integrals.map_shell_to_basis_function(shells); 
	
	const auto& buf = engine.results();
	
    for(auto s1=0; s1!=shells.size(); ++s1) {

      auto bf1_first = shell2bf[s1]; // first basis function in this shell
      auto n1 = shells[s1].size();   // number of basis functions in this shell

      for(auto s2=0; s2<=s1; ++s2) {

        auto bf2_first = shell2bf[s2];
        auto n2 = shells[s2].size();

        for(auto s3=0; s3<=s1; ++s3) {

          auto bf3_first = shell2bf[s3];
          auto n3 = shells[s3].size();

          const auto s4_max = (s1 == s3) ? s2 : s3;
          for(auto s4=0; s4<=s4_max; ++s4) {

            auto bf4_first = shell2bf[s4];
            auto n4 = shells[s4].size();

            // compute the permutational degeneracy (i.e. # of equivalents) of the given shell set
            auto s12_deg = (s1 == s2) ? 1.0 : 2.0;
            auto s34_deg = (s3 == s4) ? 1.0 : 2.0;
            auto s12_34_deg = (s1 == s3) ? (s2 == s4 ? 1.0 : 2.0) : 2.0;
            auto s1234_deg = s12_deg * s34_deg * s12_34_deg;

            engine.compute(shells[s1], shells[s2], shells[s3], shells[s4]);
            const auto* buf_1234 = buf[0];
            if (buf_1234 == nullptr)
              continue; // if all integrals screened out, skip to next quartet

            // 1) each shell set of integrals contributes up to 6 shell sets of the Fock matrix:
            //    F(a,b) += (ab|cd) * D(c,d)
            //    F(c,d) += (ab|cd) * D(a,b)
            //    F(b,d) -= 1/4 * (ab|cd) * D(a,c)
            //    F(b,c) -= 1/4 * (ab|cd) * D(a,d)
            //    F(a,c) -= 1/4 * (ab|cd) * D(b,d)
            //    F(a,d) -= 1/4 * (ab|cd) * D(b,c)
            // 2) each permutationally-unique integral (shell set) must be scaled by its degeneracy,
            //    i.e. the number of the integrals/sets equivalent to it
            // 3) the end result must be symmetrized
            for(auto f1=0, f1234=0; f1!=n1; ++f1) {
              const auto bf1 = f1 + bf1_first;
              for(auto f2=0; f2!=n2; ++f2) {
                const auto bf2 = f2 + bf2_first;
                for(auto f3=0; f3!=n3; ++f3) {
                  const auto bf3 = f3 + bf3_first;
                  for(auto f4=0; f4!=n4; ++f4, ++f1234) {
                    const auto bf4 = f4 + bf4_first;

                    const auto value = buf_1234[f1234];

                    const auto value_scal_by_deg = value * s1234_deg;

                    jints(bf1,bf2) += D(bf3,bf4) * value_scal_by_deg;
                    jints(bf3,bf4) += D(bf1,bf2) * value_scal_by_deg;
                    kints(bf1,bf3) += 0.5 * D(bf2,bf4) * value_scal_by_deg;
                    kints(bf2,bf4) += 0.5 * D(bf1,bf3) * value_scal_by_deg;
                    kints(bf1,bf4) += 0.5 * D(bf2,bf3) * value_scal_by_deg;
                    kints(bf2,bf3) += 0.5 * D(bf1,bf4) * value_scal_by_deg;
                  }
                }
              }
            }

          }
        }
      }
    }

    // symmetrize the result
    jints = 0.25 * (jints + jints.transpose());
	kints = 0.25 * (kints + kints.transpose());
	jkints = jints - 0.5*kints;
}

// Form the JK matrix from two electron integrals stored on file
void Fock::formJKfile()
{
}

// Add an error vector
void Fock::addErr(Vector e)
{
  if (iter > MAX) {
    errs.erase(errs.begin()); // Remove first element
  }
  errs.push_back(e); // Push e onto the end of errs
}	
		

void Fock::makeFock()
{
  focka = hcore + jkints;
  if (diis) { // Archive for averaging
    if (iter > MAX) {
      focks.erase(focks.begin());
    }
    focks.push_back(focka);
  }		
}

void Fock::makeFock(Matrix& jbints)
{
  focka = hcore + 0.5*(jints + jbints - kints);
  if (diis) { // Archive for averaging
    if (iter > MAX) {
      focks.erase(focks.begin());
    }
    focks.push_back(focka);
  }
}

// Perform DIIS averaging
// Direct inversion of the iterative subspace
// Greatly improves convergence and numerical behaviour
// of the scf iterations.
void Fock::DIIS()
{
  if (iter > 1) {
    int lim = (iter < MAX ? iter : MAX);
    Matrix B(lim+1, lim+1, -1.0); // Error norm matrix
    B(lim, lim) = 0.0;
    
    // The elements of B are <e_i | e_j >
    for (int i = 0; i < lim; i++){
      for (int j = i; j < lim; j++){
	B(i, j) = inner(errs[i], errs[j]);
	B(j, i) = B(i, j);
      }
    }

    // Solve the linear system of equations for the weights
    Vector w(lim+1, 0.0); w[lim] = -1.0;
    
    Eigen::MatrixXd temp(lim+1, lim+1);
    for (int i = 0; i < lim+1; i++){
      for (int j = 0; j < lim+1; j++){
	temp(i, j) = B(i, j);
      }
    }

    Eigen::JacobiSVD<Eigen::MatrixXd> svd(temp, Eigen::ComputeFullU | Eigen::ComputeFullV);
    temp = svd.singularValues().asDiagonal();
    for (int i = 0; i < lim+1; i++){
      if (fabs(temp(i, i)) > molecule.getLog().precision())
	temp(i, i) = 1.0/temp(i, i);
    }
    temp = svd.matrixV() * temp * svd.matrixU().transpose();
    
    for (int i = 0; i < lim+1; i++){
      for (int j = 0; j < lim+1; j++){
	B(i, j) = temp(i, j);
      }
    }

    w = B*w;
    
    // Average the fock matrices according to the weights
    focka.assign(nbfs, nbfs, 0.0);
    for (int i = 0; i < lim; i++) {
      focka = focka + w(i)*focks[i];   
    } 
  }
  iter++;
}

void Fock::simpleAverage(Matrix& D0, double weight)
{
  dens = weight*dens + (1.0-weight)*D0;
}
