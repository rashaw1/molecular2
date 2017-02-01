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
#include "atom.hpp"
#include "tensor4.hpp"
#include <iostream>
#include <cmath>
#include <Eigen/Eigenvalues>
#include <libint2.hpp>
#include "logger.hpp"

typedef Eigen::DiagonalMatrix<double, Eigen::Dynamic, Eigen::Dynamic>
	DiagonalMatrix;

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
	MAX = 8;
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

Fock::Fock(const Fock& other) : integrals(other.integrals), molecule(other.molecule) {
	hcore = other.hcore;
	jkints = other.jkints;
	jints = other.jints;
	kints = other.kints;
	orthog = other.orthog;
	fockm = other.fockm;
	focka = other.focka;
	CP = other.CP;
	forces = other.forces;
	hessian = other.hessian;
	eps = other.eps;
	focks = other.focks;
	dens = other.dens;
	direct = other.direct;
	twoints = other.twoints;
	fromfile = other.fromfile;
	diis = other.diis;
	nbfs = other.nbfs;
	iter = other.iter;
	MAX = other.MAX; 
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
  
void Fock::average(Vector &w) {
	if (diis && iter > 2) {
		// Average the fock matrices according to the weights
		focka.assign(nbfs, nbfs, 0.0);
		int offset = focks.size() - w.size();
		for (int i = offset; i < focks.size(); i++) {
			focka = focka + w[i-offset]*focks[i]; 
		} 
	}
}

// Transform the AO fock matrix to the MO basis 
void Fock::transform(bool first)
{
	if (first) { 
		// Form the core Fock matrix as (S^-1/2)(T)H(S^-1/2)
		fockm = (orthog.transpose()) * ( hcore * orthog);
	} else {
		// Form the orthogonalised fock matrix
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
		formJK(dens); 
	} else if (direct) {
		formJKdirect(integrals.getPrescreen(), dens);
	} else {
		try {
			formJKfile();
		} catch (Error e) {
			molecule.getLog().error(e);
		}
	}
}

// Form the 2J-K matrix, given that twoints is stored in memory
void Fock::formJK(Matrix& P)
{
	jints.assign(nbfs, nbfs, 0.0);
	kints.assign(nbfs, nbfs, 0.0);
	for (int u = 0; u < nbfs; u++){
		for (int v = 0; v < nbfs ; v++){
			for (int s = 0; s < nbfs; s++){
				for (int l = 0; l < nbfs; l++){
					jints(u, v) += P(s, l)*integrals.getERI(u, v, s, l);
					kints(u, v) += P(s, l)*integrals.getERI(u, l, s, v);
				}
			}
		}
	}
	jkints = jints - 0.5*kints;
}
// Form JK using integral direct methods
void Fock::formJKdirect(const Matrix& Schwarz, Matrix& P)
{
	using libint2::Shell;
	using libint2::Engine;
	using libint2::Operator;
	
	std::vector<Shell>& shells = molecule.getBasis().getIntShells();
	const auto n = integrals.nbasis(shells);
	jints.assign(n, n, 0.0);
	kints.assign(n, n, 0.0);
	const Matrix& D = P; 
	
	const auto do_schwarz_screen = Schwarz.ncols() != 0 && Schwarz.nrows() != 0;
	EMatrix D_shblk_norm = integrals.compute_shellblock_norm(shells, D);
	auto fock_precision = molecule.getLog().precision();
	auto maxnprim =  integrals.max_nprim(shells);
	auto maxnprim4 = maxnprim * maxnprim * maxnprim * maxnprim;
	auto engine_precision = std::min(fock_precision / D_shblk_norm.maxCoeff(), std::numeric_limits<double>::epsilon()) / maxnprim4;
	
	Engine engine(Operator::coulomb, integrals.max_nprim(shells), integrals.max_l(shells), 0);
	engine.set_precision(engine_precision);
	auto shell2bf = integrals.map_shell_to_basis_function(shells); 
	
	const auto& buf = engine.results();
	
	for(auto s1=0; s1!=shells.size(); ++s1) {

		auto bf1_first = shell2bf[s1]; // first basis function in this shell
		auto n1 = shells[s1].size();   // number of basis functions in this shell

		for(auto s2=0; s2<=s1; ++s2) {

			auto bf2_first = shell2bf[s2];
			auto n2 = shells[s2].size();
		
			const auto Dnorm12 = do_schwarz_screen ? D_shblk_norm(s1, s2) : 0.0;

			for(auto s3=0; s3<=s1; ++s3) {

				auto bf3_first = shell2bf[s3];
				auto n3 = shells[s3].size();
		  
				const auto Dnorm123 = do_schwarz_screen ? std::max(D_shblk_norm(s1, s3),
				std::max(D_shblk_norm(s2, s3), Dnorm12)) : 0.0;

				const auto s4_max = (s1 == s3) ? s2 : s3;
				for(auto s4=0; s4<=s4_max; ++s4) {

					const auto Dnorm1234 = do_schwarz_screen ? std::max(D_shblk_norm(s1, s4),
					std::max(D_shblk_norm(s2, s4), std::max(D_shblk_norm(s3, s4), Dnorm123))) : 0.0;
						
					if(do_schwarz_screen && Dnorm1234 * Schwarz(s1, s2) * Schwarz(s3, s4) < fock_precision) 
						continue;
					
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
		

void Fock::makeFock()
{
	focka = hcore + jkints;
	if (diis) { // Archive for averaging
		if (iter >= MAX) {
			focks.erase(focks.begin());
		}
		focks.push_back(focka);
		iter++;
	}		
}

void Fock::makeFock(Matrix& jbints)
{
	focka = hcore + 0.5*(jints + jbints - kints);
	if (diis) { // Archive for averaging
		if (iter >= MAX) {
			focks.erase(focks.begin());
		}
		focks.push_back(focka);
		iter++;
	}
}

void Fock::simpleAverage(Matrix& D0, double weight)
{
	dens = weight*dens + (1.0-weight)*D0;
}

void Fock::compute_forces(const std::vector<Atom> &atoms, int nocc) {
	std::cout << "Here" << std::endl << std::flush;
	using libint2::Shell;
	using libint2::Operator;

	std::vector<Shell>& shells = molecule.getBasis().getIntShells();
	EMatrix D(dens.nrows(), dens.ncols());
	for (int i = 0; i < dens.nrows(); ++i) {
		for (int j = 0; j < dens.ncols(); ++j) D(i, j) = dens(i, j);
	}
	
	Matrix F1(atoms.size(), 3, 0.0);
	Matrix F_Pulay(atoms.size(), 3, 0.0);
	
	// One-body contributions to forces
	std::cout << "Kinetic and potential.. " << std::endl << std::flush;
	auto T1 = integrals.compute_1body_ints_deriv<Operator::kinetic>(1, shells, atoms);
	auto V1 = integrals.compute_1body_ints_deriv<Operator::nuclear>(1, shells, atoms);
	for (auto atom = 0, i = 0; atom != atoms.size(); ++atom) {
		for (auto xyz = 0; xyz != 3; ++xyz, ++i) {
			auto force = 2 * (T1[i] + V1[i]).cwiseProduct(D).sum();
			F1(atom, xyz) += force;
		}
	}
	std::cout << "done." << std::endl << std::flush;
	
	// Pulay force
	EMatrix evals_occ = EMatrix::Zero(nocc, nocc);
	for (int i = 0; i < nocc; ++i) evals_occ(i, i) = eps[i];
	EMatrix C_occ(nocc, nocc);
	for (int i = 0; i < nocc; ++i){ 
		for (int j = 0; j < nocc; ++j) C_occ(i, j) = CP(i, j);
	}
	EMatrix W = C_occ * evals_occ * C_occ.transpose();
	std::cout << "Overlap.. " << std::endl << std::flush;
	auto S1 = integrals.compute_1body_ints_deriv<Operator::overlap>(1, shells, atoms);
	std::cout << "done." << std::endl << std::flush;
	for (auto atom = 0, i = 0; atom != atoms.size(); ++atom) {
		for (auto xyz = 0; xyz != 3; ++xyz, ++i) {
			auto force = 2 * S1[i].cwiseProduct(W).sum();
			F_Pulay(atom, xyz) -= force;
		}
	}
	
	// Two-body contrtibutions to forces
	Matrix F2(atoms.size(), 3, 0.0);
	auto G1 = compute_2body_fock_deriv<1>(atoms, D);
	for (auto atom = 0, i = 0; atom != atoms.size(); ++atom) {
		for (auto xyz = 0; xyz != 3; ++xyz, ++i) {
			// identity prefactor since E(HF) = trace(H + F, D) = trace(2H + G, D)
			auto force = G1[i].cwiseProduct(D).sum();
			F2(atom, xyz) += force;
		}
	}
	
	// Compute nuclear repulsion forces
	Matrix FN(atoms.size(), 3, 0.0);
	for (auto a1 = 1; a1 != atoms.size(); ++a1) {
		const auto& atom1 = atoms[a1];
		for (auto a2 = 0; a2 < a1; ++a2) {
			const auto& atom2 = atoms[a2];

			auto x12 = atom1.getX() - atom2.getX();
			auto y12 = atom1.getY() - atom2.getY();
			auto z12 = atom1.getZ() - atom2.getZ();
			auto r12_2 = x12 * x12 + y12 * y12 + z12 * z12;
			auto r12 = sqrt(r12_2);
			auto r12_3 = r12 * r12_2;

			auto z1z2_over_r12_3 =
				atom1.getEffectiveCharge() * atom2.getEffectiveCharge() / r12_3;

			auto fx = -x12 * z1z2_over_r12_3;
			auto fy = -y12 * z1z2_over_r12_3;
			auto fz = -z12 * z1z2_over_r12_3;
			FN(a1, 0) += fx;
			FN(a1, 1) += fy;
			FN(a1, 2) += fz;
			FN(a2, 0) -= fx;
			FN(a2, 1) -= fy;
			FN(a2, 2) -= fz;
		}
	}
	
	forces = F1 + F_Pulay + F2 + FN;
	forces.print();
}

void Fock::compute_hessian(const std::vector<Atom> &atoms, int nocc) {
	const auto ncoords = 3 * atoms.size();
	const auto nelem = ncoords * (ncoords +1) / 2;
	
	using libint2::Shell;
	using libint2::Operator;

	std::vector<Shell>& shells = molecule.getBasis().getIntShells();
	EMatrix D(dens.nrows(), dens.ncols());
	for (int i = 0; i < dens.nrows(); ++i) {
		for (int j = 0; j < dens.ncols(); ++j) D(i, j) = dens(i, j);
	}
	
	Matrix H1(ncoords, ncoords, 0.0);
	Matrix H_Pulay(ncoords, ncoords, 0.0);
	
	// One-body contributions to the hessian
	auto T2 = integrals.compute_1body_ints_deriv<Operator::kinetic>(2, shells, atoms);
	auto V2 = integrals.compute_1body_ints_deriv<Operator::nuclear>(2, shells, atoms);
	for (auto row = 0, i = 0; row != ncoords; ++row) {
		for (auto col = row; col != ncoords; ++col, ++i) {
			auto hess = 2 * (T2[i] + V2[i]).cwiseProduct(D).sum();
			H1(row, col) += hess;
		}
	}
	
	// Pulay hessian
	EMatrix evals_occ = EMatrix::Zero(nocc, nocc);
	for (int i = 0; i < nocc; ++i) evals_occ(i, i) = eps[i];
	EMatrix C_occ(nocc, nocc);
	for (int i = 0; i < nocc; ++i){ 
		for (int j = 0; j < nocc; ++j) C_occ(i, j) = CP(i, j);
	}
	EMatrix W = C_occ * evals_occ * C_occ.transpose();
	auto S2 = integrals.compute_1body_ints_deriv<Operator::overlap>(2, shells, atoms);
	for (auto row = 0, i = 0; row != ncoords; ++row) {
		for (auto col = row; col != ncoords; ++col, ++i) {
			auto hess = 2 * S2[i].cwiseProduct(W).sum();
			H_Pulay(row, col) -= hess;
		}
	}
	
	// Two-body contributions to the hessian
	Matrix H2(ncoords, ncoords, 0.0);
	auto G2 = compute_2body_fock_deriv<2>(atoms, D);
	for (auto row = 0, i = 0; row != ncoords; ++row) {
		for (auto col = row; col != ncoords; ++col, ++i) {
			// identity prefactor since E(HF) = trace(H + F, D) = trace(2H + G, D)
			auto hess = G2[i].cwiseProduct(D).sum();
			H2(row, col) += hess;
		}
	}
	
	// Nuclear repulsion hessian
	Matrix HN(ncoords, ncoords, 0.0);
	for (auto a1 = 1; a1 != atoms.size(); ++a1) {
		const auto& atom1 = atoms[a1];
		for (auto a2 = 0; a2 < a1; ++a2) {
			const auto& atom2 = atoms[a2];

			auto x12 = atom1.getX() - atom2.getX();
			auto y12 = atom1.getY() - atom2.getY();
			auto z12 = atom1.getZ() - atom2.getZ();
			auto x12_2 = x12 * x12;
			auto y12_2 = y12 * y12;
			auto z12_2 = z12 * z12;
			auto r12_2 = x12 * x12 + y12 * y12 + z12 * z12;
			auto r12 = sqrt(r12_2);
			auto r12_5 = r12 * r12_2 * r12_2;

			auto z1z2_over_r12_5 = atom1.getEffectiveCharge() * atom2.getEffectiveCharge() / r12_5;

			HN(3*a1 + 0, 3*a1 + 0) += z1z2_over_r12_5 * (3*x12_2 - r12_2);
			HN(3*a1 + 1, 3*a1 + 1) += z1z2_over_r12_5 * (3*y12_2 - r12_2);
			HN(3*a1 + 2, 3*a1 + 2) += z1z2_over_r12_5 * (3*z12_2 - r12_2);
			HN(3*a1 + 0, 3*a1 + 1) += z1z2_over_r12_5 * (3*x12*y12);
			HN(3*a1 + 0, 3*a1 + 2) += z1z2_over_r12_5 * (3*x12*z12);
			HN(3*a1 + 1, 3*a1 + 2) += z1z2_over_r12_5 * (3*y12*z12);

			HN(3*a2 + 0, 3*a2 + 0) += z1z2_over_r12_5 * (3*x12_2 - r12_2);
			HN(3*a2 + 1, 3*a2 + 1) += z1z2_over_r12_5 * (3*y12_2 - r12_2);
			HN(3*a2 + 2, 3*a2 + 2) += z1z2_over_r12_5 * (3*z12_2 - r12_2);
			HN(3*a2 + 0, 3*a2 + 1) += z1z2_over_r12_5 * (3*x12*y12);
			HN(3*a2 + 0, 3*a2 + 2) += z1z2_over_r12_5 * (3*x12*z12);
			HN(3*a2 + 1, 3*a2 + 2) += z1z2_over_r12_5 * (3*y12*z12);

			HN(3*a2 + 0, 3*a1 + 0) -= z1z2_over_r12_5 * (3*x12_2 - r12_2);
			HN(3*a2 + 1, 3*a1 + 1) -= z1z2_over_r12_5 * (3*y12_2 - r12_2);
			HN(3*a2 + 2, 3*a1 + 2) -= z1z2_over_r12_5 * (3*z12_2 - r12_2);
			HN(3*a2 + 1, 3*a1 + 0) -= z1z2_over_r12_5 * (3*y12*x12);
			HN(3*a2 + 2, 3*a1 + 0) -= z1z2_over_r12_5 * (3*z12*x12);
			HN(3*a2 + 2, 3*a1 + 1) -= z1z2_over_r12_5 * (3*z12*y12);
			HN(3*a2 + 0, 3*a1 + 1) -= z1z2_over_r12_5 * (3*x12*y12);
			HN(3*a2 + 0, 3*a1 + 2) -= z1z2_over_r12_5 * (3*x12*z12);
			HN(3*a2 + 1, 3*a1 + 2) -= z1z2_over_r12_5 * (3*y12*z12);
		}
	}
	hessian = H1 + H_Pulay + H2 + HN;
	std::cout << "\n\n";
	hessian.print();
	std::cout << "\n\n";
}

template <unsigned deriv_order>
std::vector<EMatrix> Fock::compute_2body_fock_deriv(const std::vector<Atom> &atoms, const EMatrix& D)
{
	using libint2::Shell;	
	using libint2::Engine;	
	using libint2::Operator;		
	
	std::vector<Shell>& shells = molecule.getBasis().getIntShells();
	const Matrix& Schwarz = integrals.getPrescreen();
	double precision = molecule.getLog().precision();
	const auto n = integrals.nbasis(shells);
	const auto nshells = shells.size();
	const auto nderiv = libint2::num_geometrical_derivatives(atoms.size(), deriv_order);
	const auto nderiv_shellset = libint2::num_geometrical_derivatives(atoms.size(), deriv_order);
	const auto ncoords_times_2 = (atoms.size() * 3) * 2;
	
	std::vector<EMatrix> G(nderiv, EMatrix::Zero(n, n));
	
	const auto do_schwarz_screen = Schwarz.ncols() != 0 && Schwarz.nrows() != 0;
	EMatrix D_shblk_norm = integrals.compute_shellblock_norm(shells, dens);
	auto fock_precision = molecule.getLog().precision();
	auto maxnprim =  integrals.max_nprim(shells);
	auto maxnprim4 = maxnprim * maxnprim * maxnprim * maxnprim;
	auto engine_precision = std::min(fock_precision / D_shblk_norm.maxCoeff(), std::numeric_limits<double>::epsilon()) / maxnprim4;
	
	Engine engine(Operator::coulomb, maxnprim, integrals.max_l(shells), deriv_order);
	engine.set_precision(engine_precision);
	
	auto shell2bf = integrals.map_shell_to_basis_function(shells); 
	auto shell2atom = integrals.map_shell_to_atom(atoms, shells);
	
	const auto& buf = engine.results();
	
	size_t shell_atoms[4];
	auto obs_shellpair_list = integrals.compute_shellpair_list(shells);
	
	// loop over permutationally-unique set of shells
	for (auto s1 = 0l, s1234 = 0l; s1 != nshells; ++s1) {
		auto bf1_first = shell2bf[s1];  // first basis function in this shell
		auto n1 = shells[s1].size();       // number of basis functions in this shell
		shell_atoms[0] = shell2atom[s1];

		for (const auto& s2 : obs_shellpair_list[s1]) {
			auto bf2_first = shell2bf[s2];
			auto n2 = shells[s2].size();
			shell_atoms[1] = shell2atom[s2];

			const auto Dnorm12 = do_schwarz_screen ? D_shblk_norm(s1, s2) : 0.0;

			for (auto s3 = 0; s3 <= s1; ++s3) {
				auto bf3_first = shell2bf[s3];
				auto n3 = shells[s3].size();
				shell_atoms[2] = shell2atom[s3];

				const auto Dnorm123 = do_schwarz_screen ? std::max(D_shblk_norm(s1, s3),
				std::max(D_shblk_norm(s2, s3), Dnorm12)) : 0.0;

				const auto s4_max = (s1 == s3) ? s2 : s3;
				for (const auto& s4 : obs_shellpair_list[s3]) {
					if (s4 > s4_max)
						break;

					const auto Dnorm1234 = do_schwarz_screen ? std::max(D_shblk_norm(s1, s4),
					std::max(D_shblk_norm(s2, s4), std::max(D_shblk_norm(s3, s4), Dnorm123))) : 0.0;

					if (do_schwarz_screen && Dnorm1234 * Schwarz(s1, s2) * Schwarz(s3, s4) < fock_precision)
						continue;

					auto bf4_first = shell2bf[s4];
					auto n4 = shells[s4].size();
					shell_atoms[3] = shell2atom[s4];

					const auto n1234 = n1 * n2 * n3 * n4;

					// compute the permutational degeneracy (i.e. # of equivalents) of
					// the given shell set
					auto s12_deg = (s1 == s2) ? 1.0 : 2.0;
					auto s34_deg = (s3 == s4) ? 1.0 : 2.0;
					auto s12_34_deg = (s1 == s3) ? (s2 == s4 ? 1.0 : 2.0) : 2.0;
					auto s1234_deg = s12_deg * s34_deg * s12_34_deg;

					// computes contribution from shell set \c idx to the operator matrix with
					// index \c op
					auto add_shellset_to_dest = [&](std::size_t op, std::size_t idx, int coord1, int coord2, double scale = 1.0) {
						auto& g = G[op];
						auto shset = buf[idx];
						const auto weight = scale * s1234_deg;

						for (auto f1 = 0, f1234 = 0; f1 != n1; ++f1) {
							const auto bf1 = f1 + bf1_first;
							for (auto f2 = 0; f2 != n2; ++f2) {
								const auto bf2 = f2 + bf2_first;
								for (auto f3 = 0; f3 != n3; ++f3) {
									const auto bf3 = f3 + bf3_first;
									for (auto f4 = 0; f4 != n4; ++f4, ++f1234) {
										const auto bf4 = f4 + bf4_first;

										const auto value = shset[f1234];
										const auto wvalue = value * weight;

										g(bf1, bf2) += D(bf3, bf4) * wvalue;
										g(bf3, bf4) += D(bf1, bf2) * wvalue;
										g(bf1, bf3) -= 0.25 * D(bf2, bf4) * wvalue;
										g(bf2, bf4) -= 0.25 * D(bf1, bf3) * wvalue;
										g(bf1, bf4) -= 0.25 * D(bf2, bf3) * wvalue;
										g(bf2, bf3) -= 0.25 * D(bf1, bf4) * wvalue;
									}
								}
							}
						}
					};

					engine.compute2<Operator::coulomb, libint2::BraKet::xx_xx, deriv_order>(shells[s1], shells[s2], shells[s3], shells[s4]);
					if (buf[0] == nullptr)
						continue; // if all integrals screened out, skip to next quartet

					switch (deriv_order) {
						case 0: {
							int coord1 = 0, coord2 = 0;
							add_shellset_to_dest(0, 0, coord1, coord2);
						} break;

						case 1: {
							for (auto d = 0; d != 12; ++d) {
								const int a = d / 3;
								const int xyz = d % 3;

								auto coord = shell_atoms[a] * 3 + xyz;
								auto& g = G[coord];

								int coord1 = 0, coord2 = 0;

								add_shellset_to_dest(coord, d, coord1, coord2);

							}  // d \in [0,12)
						} break;

						case 2: {
							// computes upper triangle index
							// n2 = matrix size times 2
							// i,j = (unordered) indices
#define upper_triangle_index(n2, i, j) (i < j ? i : j) * (n2 - (i < j ? i : j) - 1) / 2 + (i > j ? i : j)
							// look over shellsets in the order in which they appear
							std::size_t shellset_idx = 0;
							for (auto c1 = 0; c1 != 4; ++c1) {
								auto a1 = shell_atoms[c1];
								auto coord1 = 3 * a1;
								for (auto xyz1 = 0; xyz1 != 3; ++xyz1, ++coord1) {
									for (auto c2 = c1; c2 != 4; ++c2) {
										auto a2 = shell_atoms[c2];
										auto xyz2_start = (c1 == c2) ? xyz1 : 0;
										auto coord2 = 3 * a2 + xyz2_start;
										for (auto xyz2 = xyz2_start; xyz2 != 3; ++xyz2, ++coord2) {
											double scale = (coord1 == coord2 && c1 != c2) ? 2.0 : 1.0;

											const auto coord12 = upper_triangle_index(ncoords_times_2, coord1, coord2);
											add_shellset_to_dest(coord12, shellset_idx, coord1, coord2, scale);
											++shellset_idx;
										}
									}
								}
							}
						} break;
#undef upper_triangle_index

						default:
						assert(deriv_order <= 2 &&
							"support for 3rd and higher derivatives of the Fock "
								"matrix not yet implemented");
					}
				}
			}
		}
	}

	std::vector<EMatrix> GG(nderiv);
	for (auto d = 0; d != nderiv; ++d) {
		GG[d] = 0.5 * (G[d] + G[d].transpose());
	}

	return GG;
}

FockFragment::FockFragment(IntegralEngine& ints, Molecule& m) : Fock(ints, m)
{
}

FockFragment::FockFragment(const FockFragment& other) : Fock(other) {

}

