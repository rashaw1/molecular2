/*
*
*   PURPOSE: To implement class IntegralEngine, which calculates, stores,
*            and processes the molecular integrals needed in ab initio
*            quantum chemistry calculations.
*
*   DATE         AUTHOR           CHANGES
*   ======================================================================
*   02/09/15     Robert Shaw      Original code.
*   03/09/15     Robert Shaw      Merged overlap and kinetic integrals.
*   04/09/15     Robert Shaw      Now supports general contracted.
*   06/09/15     Robert Shaw      Nuclear attraction ints for prims.
*   07/09/15     Robert Shaw      formNucAttract() now works, as does
*                                 makeSpherical(ints, lnums)
*   08/09/15     Robert Shaw      Auxiliary 2e- integrals, twoe
*   09/09/15     Robert Shaw      Shell 2e- integrals, up to (m0|pq)
*                                 - need to sphericalise, increment, 
*                                 then sphericalise again.
*   10/09/15     Robert Shaw      Finished shell twoe. Got rid of 
*   							   makeContracted and makeSpherical
*                                 for 2e ints, as absorbed into twoe.
*/

#include "error.hpp"
#include "integrals.hpp"
#include "mathutil.hpp"
#include "basis.hpp"
#include "logger.hpp"
#include "tensor4.hpp"
#include "ecp.hpp"
#include "ecpint.hpp"
#include "gshell.hpp"

#include <cmath>
#include <iomanip>
#include <string>
#include <thread>   

// Constructor
IntegralEngine::IntegralEngine(Molecule& m) : molecule(m)
{
	// Calculate sizes
	int natoms = molecule.getNAtoms();
	int N = 0; // No. of cartesian basis functions
	int M = 0; // No. of spherical basis functions
	for (int i = 0; i < natoms; i++){
		N += m.getAtom(i).getNbfs();
		M += m.getAtom(i).getNSpherical();
	}
	// Cartesian is easy - there are (N^2+N)/2
	// unique 1e integrals and ([(N^2+N)/2]^2 + (N^2+N)/2)/2
	// unique 2e integrals
	int ones = (N*(N+1));
	sizes.resize(4);
	sizes[0] = ones;
	sizes[1] = (ones*(ones+1));
  
	ones = (M*(M+1));
	sizes[2] = ones;
	sizes[3] = (ones*(ones+1));

	molecule.getLog().title("INTEGRAL GENERATION");
  
	molecule.getLog().print("Forming the one electron integrals\n");
  	
	auto &shells = molecule.getBasis().getIntShells();
	std::vector<Atom> atoms;
	for (int i = 0; i < molecule.getNAtoms(); i++) atoms.push_back(molecule.getAtom(i));
	
	sints = compute_1body_ints(shells, libint2::Operator::overlap);
	tints = compute_1body_ints(shells, libint2::Operator::kinetic);
	naints = compute_1body_ints(shells, libint2::Operator::nuclear, atoms);
	
	if(molecule.getBasis().hasECPS()) {
		naints = naints + compute_ecp_ints(shells);
	}
	
	Vector lnums(naints.nrows());
	int pos = 0;
	for (auto a : atoms) {
		Vector ltemp = molecule.getBasis().getLnums(a.getCharge());
		for (int i = 0; i < ltemp.size(); i++) lnums[pos++] = ltemp[i];
	}
	Matrix newnaints = makeSpherical(naints, lnums);
	newnaints.print(); std::cout << std::endl;
	//naints.print(); std::cout << std::endl << std::endl; 
  
	molecule.getLog().print("One electron integrals complete\n");
	molecule.getLog().localTime();
    
	Vector ests = getEstimates();
	if ( molecule.getLog().getMemory() < ests(3) && !molecule.getLog().direct() ) {
		Error e("MEMERR", "Not enough memory for ERIs.");
		molecule.getLog().error(e);
		molecule.getLog().setDirect(true);
	}
	
	prescreen = compute_schwarz_ints<>(shells);
	if (prescreen.nrows() < 10) { 
		molecule.getLog().print("Forming the two electron repulsion integrals.\n");
		molecule.getLog().print("PRESCREENING MATRIX:\n");
		molecule.getLog().print(prescreen);
		molecule.getLog().print("\n\n");
	}
		
	if ( molecule.getLog().direct() ){
		molecule.getLog().print("Two electron integrals to be calculated on the fly.\n");
	} else { // Check memory requirements
		twoints = compute_eris(shells);
		
		molecule.getLog().print("Two electron integrals completed.\n");
		std::string mem = "Approximate memory usage = ";
		mem += std::to_string(M*M*M*M*sizeof(double)/(1024.0*1024.0));
		mem += " MB\n";
		molecule.getLog().print(mem);
		molecule.getLog().localTime();
		
		if (molecule.getLog().twoprint()) {
			molecule.getLog().print("Writing ERIs to file.\n");
			printERI(molecule.getLog().getIntFile(), M);
		}	
	} 
	
	//std::cout << "\n\n";
	//sints.print();
	//std::cout << "\n\n";
	//tints.print(); std::cout << "\n";
	//naints.print(); std::cout << "\n";
}

IntegralEngine::~IntegralEngine() {

}

// Accessors

// Return estimates of the memory that will be needed by the 
// one and two electron integrals. Returns as:
// [1e cart, 2e cart, 1e spher, 2e spher]
Vector IntegralEngine::getEstimates() const
{
	Vector estimates(4);
	// The amount of memory is roughly the number of integrals times the size
	// of a double in memory
	double TOMB = 1.0/(1024.0 * 1024.0);
	estimates[0] = TOMB*sizeof(double)*sizes(0);
	estimates[1] = TOMB*sizeof(double)*sizes(1);
	estimates[2] = TOMB*sizeof(double)*sizes(2);
	estimates[3] = TOMB*sizeof(double)*sizes(3);
	return estimates;
}

// Return a particular integral from twoints (taking into account symmetries)
double IntegralEngine::getERI(int i, int j, int k, int l) const
{
	return twoints(i, j, k, l);
} 

// Print a sorted list of ERIs to ostream output
void IntegralEngine::printERI(std::ostream& output, int NSpher) const
{
	output << " TWO ELECTRON INTEGRALS: \n";
	// print them out
	int icount = 0;
	int scount = 0;
	for (int c1 = 0; c1 < NSpher; c1++){
		for (int c2 = 0; c2 < c1+1; c2++){
			for (int c3 = 0; c3 < c1+1; c3++){
				for(int c4 = 0; c4 < c3+1; c4++){
					icount++;
					double multiplier = 0.125;
					if (c1!=c2) { multiplier *= 2.0; }
					if (c3!=c4) { multiplier *= 2.0; }
					if ((c1+c2) != (c3+c4)) { multiplier*=2.0; }
					else if ((c1!=c3) && (c1 != c4)) { multiplier *=2.0; }
					else if ((c2!=c3) && (c2 != c4)) { multiplier *=2.0; }
					if (fabs(getERI(c4, c3, c2, c1)) < molecule.getLog().thrint()) { scount++; multiplier = 0; }
					output << std::setw(6) << c1+1;
					output << std::setw(6) << c2+1;
					output << std::setw(6) << c3+1;
					output << std::setw(6) << c4+1;
					output << std::setw(20) << multiplier*getERI(c4, c3, c2, c1);
					output << "\n";
				}
			}
		}
	}			
	output << "N 2e Ints: " << icount << "\n";
	output << "N insig. Ints: " << scount << "\n";
}

// Contract a set of 1e- integrals
// Assumes that integrals are ordered as: 00, 01, 02, ..., 10, 11, 12, ...,
// and so on, where the first number refers to the index of c1, and the second, 
// that of c2.
double IntegralEngine::makeContracted(Vector& c1, Vector& c2, Vector& ints) const
{
	double integral = 0.0;
	int N1 = c1.size();
	int N2 = c2.size();
	// Loop over contraction coefficients
	for (int i = 0; i < N1; i++){
		for (int j = 0; j < N2; j++){
			integral += c1(i)*c2(j)*ints(i*N2+j);
		}
	}
	return integral;
}

// Sphericalise a matrix of 1e- integrals (ints)
// where the cols have angular momenta lnums.
// Returns matrix of integrals in canonical order
Matrix IntegralEngine::makeSpherical(const Matrix& ints, const Vector& lnums) const
{
	// Calculate the size of matrix needed
	int scount = 0, pcount = 0, dcount = 0, fcount = 0, gcount = 0; 
	for (int i = 0; i < lnums.size(); i++){
		switch((int)(lnums(i))){
			case 1: { pcount++; break; } // p-type
			case 2: { dcount++; break; } // d-type
			case 3: { fcount++; break; } // f-type
			case 4: { gcount++; break; } // g-type
			default: { scount++; } // Assume s-type
		}
	}

	// Number of spherical basis functions
	int M = scount + pcount + 5*(dcount/6) + 7*(fcount/10) + 9*(gcount/15); 
	// Number of cartesian basis functions.
	int N = lnums.size();

	// Declare the matrix to return the transformed integrals in
	Matrix retInts;

	// Construct a reduced list of lnums, and
	// corresponding m-quantum numbers
	Vector slnums(M); Vector smnums(M);
	int j=0, k = 0; // Counter for slnums
	while(j < N){
		switch((int)(lnums(j))){
			case 1: { // p-type 
				slnums[k] = 1; smnums[k] = 1;
				slnums[++k] = 1; smnums[k] = -1;
				slnums[++k] = 1; smnums[k++] = 0; 
				j += 3;
				break;
			} 
			case 2: { // d-type
				slnums[k] = 2; smnums[k] = 0;
				slnums[++k] = 2; smnums[k] = -2;
				slnums[++k] = 2; smnums[k] = 1;
				slnums[++k] = 2; smnums[k] = 2;
				slnums[++k] = 2; smnums[k++] = -1;
				j += 6;
				break;
			}
			case 3: { // f-type
				slnums[k] = 3; smnums[k] = 0;
				slnums[++k] = 3; smnums[k] = -1;
				slnums[++k] = 3; smnums[k] = 1;
				slnums[++k] = 3; smnums[k] = -2;
				slnums[++k] = 3; smnums[k] = 2;
				slnums[++k] = 3; smnums[k] = -3;
				slnums[++k] = 3; smnums[k++] = 3;
				j+=10;
				break;
			}
			case 4: { // g-type
				slnums[k] = 4; smnums[k] = 0;
				slnums[++k] = 4; smnums[k] = -1;
				slnums[++k] = 4; smnums[k] = 1;
				slnums[++k] = 4; smnums[k] = -2;
				slnums[++k] = 4; smnums[k] = 2;
				slnums[++k] = 4; smnums[k] = -3;
				slnums[++k] = 4; smnums[k] = 3;
				slnums[++k] = 4; smnums[k] = -4;
				slnums[++k] = 4; smnums[k++] = 4;
				j += 15;
				break;
			}
			default: { // assume s-type
				slnums[k] = 0; smnums[k++] = 0;
				j++;
			}
		}
	}

	// Make the transformation matrix
	Matrix trans(M, N, 0.0);
	// Loop through all the slnums, looking up the coefficients
	// as needed. 
	j = 0; // Column counter
	int m = 0; // m-counter
	for(int i = 0; i < M; i++){
		// Get the coefficients
		formTransMat(trans, i, j, (int)(slnums(i)), (int)(smnums(i)));
		if (m == 2*slnums(i)){ // Increment j by a suitable amount
			switch((int)(slnums(i))){
				case 1: { j+=3; break;}
				case 2: { j+=6; break;}
				case 3: { j+=10; break;}
				case 4: { j+=15; break;}
				default: { j+=1; }
			}
			m = 0;
		} else { m++; }
	}

	// Now transform the integral matrix
	retInts = trans*ints;
	retInts = retInts*(trans.transpose());
	return retInts;
}

size_t IntegralEngine::nbasis(const std::vector<libint2::Shell>& shells) {
	size_t n = 0;
	for (const auto& shell: shells)
		n += shell.size();
	return n;
}

size_t IntegralEngine::max_nprim(const std::vector<libint2::Shell>& shells) {
	size_t n = 0;
	for (auto shell: shells)
		n = std::max(shell.nprim(), n);
	return n;
}

int IntegralEngine::max_l(const std::vector<libint2::Shell>& shells) {
	int l = 0;
	for (auto shell: shells)
		for (auto c: shell.contr)
			l = std::max(c.l, l);
	return l;
}

std::vector<size_t> IntegralEngine::map_shell_to_basis_function(const std::vector<libint2::Shell>& shells) {
	std::vector<size_t> result;
	result.reserve(shells.size());

	size_t n = 0;
	for (auto shell: shells) {
		result.push_back(n);
		n += shell.size();
	}

	return result;
}

std::vector<long> IntegralEngine::map_shell_to_atom(const std::vector<Atom>& atoms, const std::vector<libint2::Shell>& shells) {
	std::vector<long> result;
	result.reserve(shells.size());
	for(const auto& s: shells) {
		auto a = std::find_if(atoms.begin(), atoms.end(), [&s](const Atom& a){ return s.O[0] == a.getX() && s.O[1] == a.getY() && s.O[2] == a.getZ(); } );
		result.push_back( a != atoms.end() ? a - atoms.begin() : -1);
	}
	return result;
}

Matrix IntegralEngine::compute_1body_ints(const std::vector<libint2::Shell>& shells,
libint2::Operator obtype,
const std::vector<Atom>& atoms)
{
	using libint2::Shell;
	using libint2::Engine;
	using libint2::Operator;

	const auto n = nbasis(shells);
	EMatrix result(n,n);

	// construct the overlap integrals engine                                                                                                                             
	Engine engine(obtype, max_nprim(shells), max_l(shells), 0);

	// nuclear attraction ints engine needs to know where the charges sit ...                                                                                                                  
	// the nuclei are charges in this case; in QM/MM there will also be classical charges                                                                                                      
	if (obtype == Operator::nuclear) {
		std::vector<std::pair<double,std::array<double,3>>> q;
		for(const auto& atom : atoms) {
			q.push_back( {static_cast<double>(atom.getEffectiveCharge()), {{atom.getX(), atom.getY(), atom.getZ()}}} );
		}
			
		engine.set_params(q);
	}

	auto shell2bf = map_shell_to_basis_function(shells);

	// buf[0] points to the target shell set after every call  to engine.compute()                                                                                                             
	const auto& buf = engine.results();

	// loop over unique shell pairs, {s1,s2} such that s1 >= s2                                                                                                                                
	// this is due to the permutational symmetry of the real integrals over Hermitian operators: (1|2) = (2|1)                                                                                 
	for(auto s1=0; s1!=shells.size(); ++s1) {

		auto bf1 = shell2bf[s1]; // first basis function in this shell                                                                                                                           
		auto n1 = shells[s1].size();

		for(auto s2=0; s2<=s1; ++s2) {

			auto bf2 = shell2bf[s2];
			auto n2 = shells[s2].size();

			// compute shell pair; return is the pointer to the buffer                                                                                                                             
			engine.compute(shells[s1], shells[s2]);

			// "map" buffer to a const Eigen Matrix, and copy it to the corresponding blocks of the result                                                                                         
			Eigen::Map<const EMatrix> buf_mat(buf[0], n1, n2);
			result.block(bf1, bf2, n1, n2) = buf_mat;
			if (s1 != s2) // if s1 >= s2, copy {s1,s2} to the corresponding {s2,s1} block, note the transpose!                                                                                     
				result.block(bf2, bf1, n2, n1) = buf_mat.transpose();

		}
	}
  
	Matrix ret(n, n);
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j)
			ret(i, j) = result(i, j);
	}
	return ret;
}

S8EvenTensor4 IntegralEngine::compute_eris(const std::vector<libint2::Shell>& shells) {
	using libint2::Shell;
	using libint2::Engine;
	using libint2::Operator;

	const auto n = nbasis(shells);
	S8EvenTensor4 eris(n, 0.0);

	Engine engine(Operator::coulomb, max_nprim(shells), max_l(shells), 0);

	auto shell2bf = map_shell_to_basis_function(shells);

	const auto& buf = engine.results();

	// loop over shell quartets
	for (auto s1=0; s1 != shells.size(); ++s1) {
    
		auto bf1_first = shell2bf[s1];
		auto n1 = shells[s1].size();

		for (auto s2=0; s2 <= s1; ++s2) {
      
			auto bf2_first = shell2bf[s2];
			auto n2 = shells[s2].size();
      
			for (auto s3=0; s3 <= s1; ++s3) {
	
				auto bf3_first = shell2bf[s3];
				auto n3 = shells[s3].size();

				const auto s4_max = (s1 == s3) ? s2 : s3;
				for(auto s4=0; s4<=s4_max; ++s4) {
	  
					auto bf4_first = shell2bf[s4];
					auto n4 = shells[s4].size();

					engine.compute(shells[s1], shells[s2], shells[s3], shells[s4]);

					const auto* buf_1234 = buf[0];
					if (buf_1234 == nullptr)
						continue;

					for(auto f1=0, f1234=0; f1!=n1; ++f1) {
            
						const auto bf1 = f1 + bf1_first;
            
						for(auto f2=0; f2!=n2; ++f2) {
              
							const auto bf2 = f2 + bf2_first;
             
							for(auto f3=0; f3!=n3; ++f3) {
                
								const auto bf3 = f3 + bf3_first;
             
								for(auto f4=0; f4!=n4; ++f4, ++f1234) {
                  
									const auto bf4 = f4 + bf4_first;
									eris(bf1, bf2, bf3, bf4) = buf_1234[f1234];
									eris(bf3, bf4, bf1, bf2) = buf_1234[f1234];
		  
								}
							}
						}
					}
				}
			}
		}
	}
	return eris;
}

EMatrix IntegralEngine::compute_shellblock_norm(const std::vector<libint2::Shell> &shells, const Matrix& A) {
	const auto nsh = shells.size();
	EMatrix Ash(nsh, nsh);
	EMatrix EA(A.nrows(), A.ncols()); 
	for (int i = 0; i < A.nrows(); ++i) {
		for (int j = 0; j < A.ncols(); ++j) EA(i, j) = A(i, j);
	}
	
	auto shell2bf = map_shell_to_basis_function(shells);
	for (size_t s1 = 0; s1 != nsh; ++s1) {
		const auto& s1_first = shell2bf[s1];
		const auto& s1_size = shells[s1].size();
		for (size_t s2 = 0; s2 != nsh; ++s2) {
			const auto& s2_first = shell2bf[s2];
			const auto& s2_size = shells[s2].size();

			Ash(s1, s2) = EA.block(s1_first, s2_first, s1_size, s2_size).lpNorm<Eigen::Infinity>();
		}
	}

	return Ash;
}

template <libint2::Operator Kernel>
Matrix IntegralEngine::compute_schwarz_ints( const std::vector<libint2::Shell> &bs1, const std::vector<libint2::Shell> _bs2, 
 bool use_2norm, typename libint2::operator_traits<Kernel>::oper_params_type params) {
     const std::vector<libint2::Shell>& bs2 = (_bs2.empty() ? bs1 : _bs2);
     const auto nsh1 = bs1.size();
     const auto nsh2 = bs2.size();
     const auto bs1_equiv_bs2 = (&bs1 == &bs2);

     Matrix K(nsh1, nsh2, 0.0);

     // construct the 2-electron repulsion integrals engine
     using libint2::Engine;

     // !!! very important: cannot screen primitives in Schwarz computation !!!
     auto epsilon = 0.;
     Engine engine(Kernel, std::max(max_nprim(bs1), max_nprim(bs2)), 
	 	std::max(max_l(bs1), max_l(bs2)), 0, epsilon, params);

	 const auto& buf = engine.results();

	 // loop over permutationally-unique set of shells
	 for (auto s1 = 0l, s12 = 0l; s1 != nsh1; ++s1) {
		 auto n1 = bs1[s1].size();  // number of basis functions in this shell

		 auto s2_max = bs1_equiv_bs2 ? s1 : nsh2 - 1;
		 for (auto s2 = 0; s2 <= s2_max; ++s2, ++s12) {
			 auto n2 = bs2[s2].size();
			 auto n12 = n1 * n2;

			 engine.compute2<Kernel, libint2::BraKet::xx_xx, 0>(bs1[s1], bs2[s2],
			 bs1[s1], bs2[s2]);

			 // the diagonal elements are the Schwarz ints ... use Map.diagonal()
			 Eigen::Map<const EMatrix> buf_mat(buf[0], n12, n12);
			 auto norm2 = use_2norm ? buf_mat.diagonal().norm()
				 : buf_mat.diagonal().lpNorm<Eigen::Infinity>();
			 K(s1, s2) = std::sqrt(norm2);
			 if (bs1_equiv_bs2) K(s2, s1) = K(s1, s2);
		 }
	 }

     return K;
 }
 
 shellpair_list_t IntegralEngine::compute_shellpair_list(const std::vector<libint2::Shell>& bs1, 
 const std::vector<libint2::Shell> _bs2, double threshold) {
 		
	 const std::vector<libint2::Shell> &bs2 = (_bs2.empty() ? bs1 : _bs2);
	 const auto nsh1 = bs1.size();
	 const auto nsh2 = bs2.size();
	 const auto bs1_equiv_bs2 = (&bs1 == &bs2);

	 // construct the 2-electron repulsion integrals engine
	 using libint2::Engine;
	 using libint2::Operator;
	 Engine engine(Operator::overlap, std::max(max_nprim(bs1), max_nprim(bs2)),
	 std::max(max_l(bs1), max_l(bs2)), 0);

	 shellpair_list_t result;
	 const auto& buf = engine.results();

	 // loop over permutationally-unique set of shells
	 for (auto s1 = 0l, s12 = 0l; s1 != nsh1; ++s1) {
		 if (result.find(s1) == result.end())
			 result.insert(std::make_pair(s1, std::vector<size_t>()));

		 auto n1 = bs1[s1].size();  // number of basis functions in this shell

		 auto s2_max = bs1_equiv_bs2 ? s1 : nsh2 - 1;
		 for (auto s2 = 0; s2 <= s2_max; ++s2, ++s12) {
	
			 auto on_same_center = (bs1[s1].O == bs2[s2].O);
			 bool significant = on_same_center;
			 if (not on_same_center) {
				 auto n2 = bs2[s2].size();
				 engine.compute(bs1[s1], bs2[s2]);
				 Eigen::Map<const EMatrix> buf_mat(buf[0], n1, n2);
				 auto norm = buf_mat.norm();
				 significant = (norm >= threshold);
			 }

			 if (significant) {
				 result[s1].emplace_back(s2);
			 }
		 }
	 }

	 // resort shell list in increasing order, i.e. result[s][s1] < result[s][s2]
	 // if s1 < s2
	 for (auto s1 = 0l; s1 != nsh1; ++s1) {
		auto& list = result[s1];
		std::sort(list.begin(), list.end());
	 }

	 return result;
 }
 
 Matrix IntegralEngine::compute_ecp_ints(const std::vector<libint2::Shell>& shells, int deriv_order) {
 	const auto n = nbasis(shells);
	Matrix ecps(n, n, 0.0);
	std::cout << n << std::endl;
	
	// Initialise ecp integral engine
	molecule.getLog().print("\nIntialising ECP integral calculations...\n");
	ECPIntegral ecpint(molecule.getECPBasis(), molecule.getBasis().getMaxL(), molecule.getECPBasis().getMaxL(), deriv_order);
	molecule.getLog().localTime();
	
	ECPBasis& ecpset = molecule.getECPBasis();
	auto shell2bf = map_shell_to_basis_function(shells);
	
	// loop over shells
	for(auto s1=0; s1!=shells.size(); ++s1) {

		auto bf1 = shell2bf[s1];
		auto n1 = shells[s1].size();
		
		double A[3] = { shells[s1].O[0], shells[s1].O[1], shells[s1].O[2] };
		GaussianShell shellA(A, shells[s1].contr[0].l);
		for (auto c : shells[s1].contr)
			for (int i = 0; i < c.coeff.size(); i++) 
				shellA.addPrim(shells[s1].alpha[i], c.coeff[i]);

		for(auto s2=0; s2<=s1; ++s2) {

			auto bf2 = shell2bf[s2];
			auto n2 = shells[s2].size();
			
			double B[3] = { shells[s2].O[0], shells[s2].O[1], shells[s2].O[2] };
			GaussianShell shellB(B, shells[s2].contr[0].l); 
			for (auto c : shells[s2].contr)
				for (int i = 0; i < c.coeff.size(); i++)
					shellB.addPrim(shells[s2].alpha[i], c.coeff[i]);
			
			Matrix shellPairInts = ecpint.compute_pair(shellA, shellB);
			//shellPairInts.print(); std::cout << "\n\n";
			for (int i = bf1; i < bf1 + shellPairInts.nrows(); i++) {
				for (int j = bf2; j < bf2 + shellPairInts.ncols(); j++) {
					ecps(i, j) = shellPairInts(i-bf1, j-bf2);
					ecps(j, i) = ecps(i, j);
				}
			}
			
			bf2 += n2; 
		}
		
		bf1 += n1; 
	}
	
	ecps.print();
	std::cout << std::endl << std::endl;
	return ecps;
 }
 
