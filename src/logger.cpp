/*
*
*   PURPOSE: To implement the file logger.hpp, defining the class Logger, a communication and
*            storage device used throughout the molecular suite.
* 
*   DATE          AUTHOR           CHANGES
*   =========================================================================================
*   29/08/15      Robert Shaw      Original code.
*   30/08/15      Robert Shaw      Implemented the many, many print functions
*/

// Includes
#include <algorithm>
#include <iomanip>
#include <ctime>
#include <map>
#include "logger.hpp"
#include "ProgramController.hpp"
#include "molecule.hpp"
#include "bf.hpp"
#include "fock.hpp"
#include "pbf.hpp"
#include "error.hpp"
#include "ioutil.hpp"
#include <libint2.hpp>
#include "ecpint.hpp"

// Define static constants
const double Logger::RTOCM = 60.19969093279;
const double Logger::RTOGHZ = 1804.74133155814269;
const double Logger::TOKCAL = 627.509469;
const double Logger::TOEV = 27.21138505;
const double Logger::TOKJ = 2625.49962;
const double Logger::TOANG = 0.52917721092;
const double Logger::TOBOHR = 1.889726124565;
const double Logger::TOWAVENUMBERS = 32298.63331; 

// Constructor
Logger::Logger(ProgramController& _control, std::ofstream& out, std::ostream& e) : control(_control), outfile(out), errstream(e)
{
	// Timer is started on initialisation of logger.
	first_time = std::chrono::steady_clock::now(); 
	last_time = first_time;

	// Initialise error array - fixed maximum of 20 errors, because any more than that and the 
	// whole idea is pointless really! Dynamic memory seems an unecessary overhead for just storing
	// errors.
	errs = new Error[20];
	nerr = 0; // No errors as of yet (hopefully)!
}

void Logger::init_intfile() {
	if (control.get_option<bool>("printeris")) { 
		std::string intfilename = control.get_option<std::string>("intfile"); 
		intfile.open(intfilename);
		if (!intfile.is_open()){
			Error e1("FILEIO", "Unable to open integral file.");
			error(e1);
		}
	}
}

void Logger::init_optfile() {
	if (!optfile.is_open()) {
		std::string optfilename = control.getName() + ".info";
		optfile.open(optfilename); 
		if (!optfile.is_open()) { 
			Error e1("FILEIO", "Unable to open optimisation file."); 
			error(e1); 
		} 
	}
}

Logger& Logger::operator=(const Logger& other) {
	control = other.control;
	nerr = other.nerr; 
	if (nerr > 0) {
		errs = new Error[nerr];
		for (int i = 0; i < nerr; i++) errs[i] = other.errs[i]; 
	}
	first_time = other.first_time;
	last_time = other.last_time;
	return *this;
}

// Destructor - get rid of the atoms and errors arrays
// and close intfile if necessary
Logger::~Logger()
{
	delete[] errs;
	if (intfile.is_open()){
		intfile.close();
	}
}

// Overloaded print functions

// Print out a basic string message
void Logger::print(const std::string& msg) const
{
	outfile << msg << "\n";
}

// Print out a vector to a given precision, either horizontally or vertically
void Logger::print(const Vector& v, int digits, bool vertical) const
{
	std::string ender = "";
	if(vertical){
		ender = "\n";
	}
	for (int i = 0; i < v.size(); i++){
		outfile << std::fixed << std::setprecision(digits) << std::setw(digits+5) << v[i] << ender;
		if (!vertical && i % 10 == 9) outfile << std::endl; 
	}
	outfile << std::endl;
}

// Print out a matrix, row by row, to a given precision
void Logger::print(const Matrix& m, int digits) const
{
	// Print out row by row
	for (int i = 0; i < m.rows(); i++){
		Vector temp = m.row(i);
		print (temp, digits, false);
	}
}

// Print out the basis set specification in the format:
// (example is C atom in 3-21G basis set)
// BASIS: 6-31G
// Total no. of cgbfs: 9
// Total no. of primitives: 9
// ===============
//  Specification
// ===============
// Atom      Shell     #CGBFs    #Prims
// ......................................
// C         s          3         6
//           p          6         3
// (if full = true, then continue with this)
// ===============
// Basis Functions
// ===============
// Atom      Shell     BF    Coeff        Exponent
// .................................................
// C         s         1     0.0617669    172.2560
//                           0.358794     25.91090
//                           0.700713     5.533350
// etc...
void Logger::print(Basis& b, bool full) const
{
	// Collect the data needed for printing
	int nbfs = b.getNBFs(); // Store number of cgbfs and prims
	int nprims = 0;

	iVector qs = b.getCharges();
	// Sort the qs and get rid of duplicates
	std::sort(qs.data(), qs.data()+qs.size(), [](int lhs, int rhs){ return rhs > lhs; });
	iVector qtemp(qs.size());
	qtemp[0] = qs(0);
	int k = 1;
	for (int i = 1; i < qs.size(); i++){
		if (qs(i) != qtemp[k-1]){
			qtemp[k] = qs(i);
			k++;
		}
	}
	qs = qtemp;
	qs.conservativeResize(k);

	// Now sum over all basis functions to get the number of prims
	for (int i = 0; i < nbfs; i++)
		nprims += b.getBF(i).getNPrims();


	// Start printing
	title("Basis Set");
	outfile << "DEFAULT: " << b.getName() << "\n";
	for (int i = 0; i < qs.size(); i++)
		outfile << getAtomName(qs[i]) << " " << b.getName(qs[i]) << "\n";
	outfile << "Total no. of cgbfs: " << nbfs << "\n";
	outfile << "Total no. of prims: " << nprims << "\n";
	title("Specification");
	outfile << std::setw(8) << "Atom";
	outfile << std::setw(8) << "Shell";
	outfile << std::setw(8) << "#CGBFs";
	outfile << std::setw(8) << "#Prims\n";
	outfile << std::string(35, '.') << "\n";
  
	// loop over the atom types
	outfile << std::setprecision(2);
	iVector subshells; iVector sublnums; 
	for (int i = 0; i < k; i++){
		int nc = 0; int np = 0;
		outfile << std::setw(8) << getAtomName(qs(i));
		subshells = b.getShells(qs[i]);

		sublnums = b.getLnums(qs[i]);
		outfile << std::setw(8) << getShellName(sublnums[0]);
		outfile << std::setw(8) << subshells[0];

		for (int j = 0; j < subshells[0]; j++){
			np += b.getBF(qs[i], j).getNPrims();
		}

		nc += subshells[0];
		outfile << std::setw(8) << np << "\n";
		for (int j = 1; j < subshells.size(); j++){
			outfile << std::setw(8) << "";
			outfile << std::setw(8) << getShellName(sublnums[j]);
			outfile << std::setw(8) << subshells[j];
			np = 0;
			for (int l = 0; l < subshells[j]; l++){
				np += b.getBF(qs[i], nc + l).getNPrims();	
			}
			nc += subshells[j];
			outfile << std::setw(8) << np << "\n";
		}

	}
	outfile << std::setprecision(8);
	// Now print out basis functions if required
	if (full) {
		title("Basis Functions");
		outfile << std::setw(8) << "Atom";
		outfile << std::setw(8) << "Shell";
		outfile << std::setw(5) << "BF";
		outfile << std::setw(18) << "Coeff";
		outfile << std::setw(18) << "Exponent\n";
		outfile << std::string(58, '.') << "\n";
		// Loop over all the basis functions
		iVector subshell; iVector sublnums; 
		Vector coeffs; Vector exps;
		std::string filler = "";
		for (int i = 0; i < k; i++){
			subshell = b.getShells(qs(i));
			sublnums = b.getLnums(qs(i));

			// Loop over shells
			int sum = 0;
			for (int r = 0; r < subshell.size(); r++){ 
				// Loop over bfs
				for (int s = 0; s < subshell[r]; s++){
					BF& bftemp = b.getBF(qs(i), s+sum);
					coeffs = bftemp.getCoeffs();
					exps = bftemp.getExps();

					// Loop over coeffs/exps
					for (int t = 0; t < coeffs.size(); t++){
						filler = ((r == 0 && s==0 && t==0) ? getAtomName(qs[i]) : "");
						outfile << std::setw(8) << filler;
						filler = ((s == 0 && t == 0) ? getShellName(sublnums[r]) : "");
						outfile << std::setw(8) << filler;
						filler = (t == 0 ? std::to_string(s+1) : "");
						outfile << std::setw(5) << filler;
						outfile << std::setw(18) << std::setprecision(8) << coeffs(t);
						outfile << std::setw(18) << std::setprecision(8) << exps(t) << "\n";
					}
				}
				sum += subshell[r];
			}
		}
	}
}      

// Print out the details of an atom, taking the form:
// Atom Type   Atomic Number    Atomic Mass(amu)  #CGBFS   x    y    z
void Logger::print(const Atom& a) const
{
	int q = a.getCharge();
	Vector c = a.getCoords();
	outfile << std::fixed << std::setprecision(6);
	outfile << std::setw(10) << getAtomName(q);
	outfile << std::setw(10) << a.getEffectiveCharge();
	outfile << std::setw(10) << a.getMass();
	outfile << std::setw(10) << a.getNbfs();
	outfile << std::setw(4) << "(" << std::setw(8) << c(0);
	outfile << ", " << std::setw(8) << c(1);
	outfile << ", " << std::setw(8) << c(2) << ")\n";
}

// Print out details of Molecule, taking the form:
// ========
// MOLECULE
// ========
// No. of e- = nel,  charge = q,  singlet/doublet/triplet etc.
// Enuc = nuclear energy
// (if inertia = true then also prints the section:
// ............................
// Principal Moments of Inertia
// ............................
// Ia = ...,  Ib = ...,  Ic = ...
// Rotational type: ...
// ....................
// Rotational Constants
// ....................
// (the coordinate system is also transformed to inertial coords)
// =====
// ATOMS
// =====
// Atom     z     Mass (a.m.u)     Coordinates
// ...............................................................
// C        6     12.014           (0.0, 3.53, -1.24)
// etc...
void Logger::print(Molecule& mol, bool inertia) const
{
	title("Molecule");
	// Print out basic details
	outfile << "# electrons = " << mol.getNel() << ",  ";
	outfile << "charge = " << mol.getCharge() << ",  ";
	std::string temp;
	// Get the state type
	switch (mol.getMultiplicity()) {
		case 1: { temp = "Singlet"; break; }
		case 2: { temp = "Doublet"; break; }
		case 3: { temp = "Triplet"; break; }
		case 4: { temp = "Quartet"; break; }
		case 5: { temp = "Quintet"; break; }
		default: { temp = "Spintacular"; break; } // Are hextets even a thing?!
	}
	outfile << temp << "\n";
	outfile << "ENUC = " << mol.getEnuc() << " Hartree\n";
  
	// Print inertial details if needed, and rotate
	// into the inertial coordinate system
	if (inertia) {
		// Get it
		temp = mol.rType();
		Vector rconsts(3);
		rconsts = mol.rConsts(1); // MHz
		Vector inert(3);
		inert = mol.getInertia(false);
		// Print it out
		outfile << std::string(30, '.') << "\n";
		outfile << "Principal Moments of Inertia\n";
		outfile << std::string(30, '.') << "\n";
		outfile << "Ia = " << std::setw(12) << inert(0);
		outfile << ",  Ib = " << std::setw(12) << inert(1);
		outfile << ",  Ic = " << std::setw(12) << inert(2) << "\n";
		outfile << "Rotational type: " << temp << "\n";
		outfile << std::string(29, '.') << "\n";
		outfile << "Rotational Constants / GHz\n";
		outfile << std::string(29, '.') << "\n";
		outfile << "A = " << std::setw(12) << rconsts(0);
		outfile << ",  B = " << std::setw(12) << rconsts(1);
		outfile << ",  C = " << std::setw(12) << rconsts(2) << "\n";
    
	}
  
	// Finally, print out all the atoms
	title("Atoms");
	outfile << std::setw(10) << "Atom" << std::setw(10) << "z";
	outfile << std::setw(10) << "Mass" << std::setw(10) << "#CGBFs"; 
	outfile << std::setw(30) << "Coordinates" << "\n";
	outfile << std::string(70, '.') << "\n";
	for (int i = 0; i < mol.getNAtoms(); i++)
		print(mol.getAtom(i));
}

// Print out a contracted gaussian basis function, in the form:
// #Prims = ..., s/px/py/pz/dxy/... type orbital
// Coefficients: ......
// Exponents: ......
void Logger::print(BF& bf) const
{
	// Get information
	Vector coef = bf.getCoeffs();
	Vector exp = bf.getExps();
	int lx = bf.getLx(); int ly = bf.getLy(); int lz = bf.getLz();
  
	// Work out the orbital type
	std::string temp;
	switch(lx + ly + lz){
		case 0: { temp = "s"; break; }
		case 1: { 
			if (lx == 1){
				temp = "px";
			} else if (ly == 1){
				temp = "py";
			} else {
				temp = "pz";
			}
			break;
		}
		case 2: { 
			if (lx == 1 && ly == 1){
				temp = "dxy";
			} else if (lx == 1 && lz == 1){
				temp = "dxz";
			} else if (ly == 1 && lz == 1){
				temp = "dyz";
			} else if (lz == 2) {
				temp = "dz2";
			} else {
				temp = "d(x2 - y2)";
			}
			break;
		}
		case 3: { temp = "f"; break; } 
		case 4: { temp = "g"; break; }
		case 5: { temp = "h"; break; }
		case 6: { temp = "Very angular"; break; }
	}

	// Print it out
	outfile << "#Primitives = " << coef.size();
	outfile << ",  " << temp << " type basis function\n";
	outfile << std::setw(15) << "Coefficients: ";
	for (int i = 0; i < coef.size(); i++){
		outfile << std::setw(12) << std::setprecision(9) << coef(i);
		if ((i % 4) == 0 && i != 0){ // Print 4 per line
			outfile << "\n" << std::setw(15) << "";
		}
	}
	outfile << "\n";
	outfile << std::setw(15) <<  "Exponents: ";
	for (int i = 0; i < exp.size(); i++){
		outfile << std::setw(12) << std::setprecision(9) << exp(i);
		if ( (i%4) == 0 && i!=0) {
			outfile << "\n" << std::setw(15) << "";
		}
	}
}

// Print out a primitive gaussian basis function in the form:
// Lx = ..., Ly = ..., Lz = ...
// norm = ..., exponent = ...
void Logger::print(const PBF& pbf) const
{
	outfile << "Lx = " << pbf.getLx();
	outfile << ",  Ly = " << pbf.getLy();
	outfile << ",  Lz = " << pbf.getLz() << "\n";
	outfile << "Norm = " << pbf.getNorm();
	outfile << ",  Exponent = " << pbf.getExponent() << "\n";
}
 

// Specialised printing routines

// Print a title like so:
//
// 
// =====
// TITLE
// =====
//
void Logger::title(const std::string& msg) const 
{
	std::string temp = msg;
	// Make upper case
	std::transform(temp.begin(), temp.end(), temp.begin(), ::toupper);
  
	// Print
	outfile << "\n\n";
	outfile << std::string(temp.length(), '=') << "\n";
	outfile << temp << "\n";
	outfile << std::string(temp.length(), '=') << "\n";
	outfile << "\n";
}

// Print a result like so:
//
// **********************
// result goes here
// **********************
// 
void Logger::result(const std::string& msg) const
{
	outfile << "\n";
	outfile << std::string(msg.length(), '*') << "\n";
	outfile << msg << "\n";
	outfile << std::string(msg.length(), '*') << "\n\n";
}

void Logger::result(const std::string& name, const double value, const std::string& units) const
{
	int length = name.length() + 25 + units.length();
	outfile << "\n";
	outfile << std::string(length, '*') << "\n";
	outfile << name << " = " <<  std::setprecision(12) << value << " " << units << "\n";
	outfile << std::string(length, '*') << "\n\n";
}

// Log an error
void Logger::error(Error& e)
{
	nerr++;
	errs[nerr%20] = e;
	errstream << "ERROR: " << e.getCode() << "\n";
	errstream << "Message: " << e.getMsg() << "\n";
	errTime();
}


// Print out an iteration table
// Initialise table header first
void Logger::initIteration()
{
	outfile << "\n";
	outfile << std::setw(12) << "Iteration";
	outfile << std::setw(24) << "Energy";
	outfile << std::setw(24) << "Delta E";
	outfile << std::setw(24) << "Delta D";
	outfile << std::setw(20) << "Time elapsed\n";
	outfile << std::string(107, '-') << "\n";
}

void Logger::initIterationCC()
{
	outfile << "\n";
	outfile << std::setw(12) << "Iteration";
	outfile << std::setw(24) << "Energy";
	outfile << std::setw(24) << "Delta E";
	outfile << std::setw(20) << "Delta T1";
	outfile << std::setw(20) << "Delta T2";
	outfile << std::setw(15) << "Time Interm";
	outfile << std::setw(15) << "Time Amps";
	outfile << std::setw(15) << "Time total\n";
	outfile << std::string(150, '-') << "\n";
}

void Logger::initALMOTable()
{
	outfile << "\n";
	outfile << std::setw(12) << "Frag. A";
	outfile << std::setw(12) << "Frag. B";
	outfile << std::setw(20) << "RCOM (Bohr)"; 
	outfile << std::setw(24) << "Disp. (kcal/mol)"; 
	outfile << std::setw(24) << "Disp. exch.";
	outfile << std::setw(20) << "Time elapsed\n";
	outfile << std::string(115, '-') << "\n";
}

// Print a single iteration
void Logger::iteration(int iter, double energy, double delta, double dd)
{
	outfile << std::setw(12) << iter;
	outfile << std::setw(24) << std::setprecision(12) << energy;
	outfile << std::setw(24) << delta;
	outfile << std::setw(24) << dd;
	outfile << std::setw(20) << std::setprecision(6) << getLocalTime();
	outfile << "\n";
  
	flush();
}

void Logger::ALMORow(int f1, int f2, double sep, double edisp, double edispexch)
{
	outfile << std::setw(12) << f1+1;
	outfile << std::setw(12) << f2+1;
	outfile << std::setw(20) << std::setprecision(3) << sep;
	outfile << std::setw(24) << std::setprecision(8) << edisp * TOKCAL; 
	outfile << std::setw(24) << edispexch * TOKCAL; 
	outfile << std::setw(20) << std::setprecision(6) << getLocalTime();
	outfile << "\n";
  
	flush();
}

// Print a single iteration
void Logger::iterationCC(int iter, double energy, double delta_e, double delta_s, double delta_d, double t_interm, double t_amps, double t_iter)
{
	outfile << std::setw(12) << iter;
	outfile << std::setw(24) << std::setprecision(12) << energy;
	outfile << std::setw(24) << delta_e;
	outfile << std::setw(20) << std::setprecision(9) << delta_s;
	outfile << std::setw(20) << delta_d;
	outfile << std::setw(15) << std::setprecision(6) << t_interm;
	outfile << std::setw(15) << t_amps;
	outfile << std::setw(15) << t_iter;
	outfile << "\n";
  
	flush();
}

// Print out the orbitals from an SCF calculation
void Logger::orbitals(const Vector& eps, int nel, bool one)
{
	if (!one) {
		outfile << "ORBITALS (Energies in Hartree)\n\n";
	}
	int nlines = (eps.size()+1)/2;
	int i = 0;
	while (i < nlines-1){
		outfile << std::setw(12) << i+1;
		outfile << std::setw(15) << eps(i);
		outfile <<  std::setw(12) << nlines+i+1;
		outfile << std::setw(15) << eps(nlines+i);
		outfile << "\n";
		i++;
	}
	outfile << std::setw(12) << i+1;
	outfile << std::setw(15) << eps(i);
	if (eps.size()%2 == 0) { 
		outfile << std::setw(12) << nlines+i+1;
		outfile << std::setw(15) << eps(nlines+i);
	}
	outfile << "\n\n";
	if (nel > 0) {
		i = (one ? nel : nel/2);
		outfile << std::setw(12) <<"HOMO:";
		outfile << std::setw(12) <<  i;
		outfile << std::setw(15) << eps(i-1)*TOEV << " eV\n";  
		if (i < eps.size()) { 	 
			outfile << std::setw(12) << "LUMO:";
			outfile << std::setw(12) << i+1;
			outfile << std::setw(15) << eps(i)*TOEV << " eV\n";
		} else {
			outfile << "MINIMAL BASIS USED, SO NO LUMO" << std::endl; 
		}
	}	
}

void Logger::coefficient_matrix(const Vector& eps, int nel, const Matrix& coeffs, bool one) 
{
	
	int nlines = nel; 
	if (!one) {
		outfile << "ORBITALS (Energies in Hartree)\n\n";
		nlines /= 2; 
	}
	
	int nsublines = coeffs.rows(); 
	nsublines = nsublines % 10 == 0 ? nsublines / 10 : nsublines / 10 + 1; 
	
	outfile << std::setw(12) << "Orbital";
	outfile << std::setw(15) << "Energy";
	outfile << std::setw(12) << "Occ."; 
	outfile << std::setw(15) << "Coefficients" << std::endl; 
	
	for (int i = 0; i < nlines; i++) {
		outfile << std::setw(12) <<  i+1;
		outfile << std::setw(15) << std::setprecision(6) << eps(i);
		outfile << std::setw(12) << (one ? 1 : 2); 
		
		int ctr = 0; 
		for (int j = 0; j < nsublines; j++) {
			int maxk = coeffs.rows() - ctr; 
			maxk = maxk < 10 ? maxk : 10; 
			for (int k = 0; k < maxk; k++)
				outfile << std::setw(15) << coeffs(k+ctr, i); 
			outfile << std::endl; 
			outfile << std::setw(12) << " "; 
			outfile << std::setw(15) << " ";
			outfile << std::setw(12) << " ";
			ctr+= 10; 
		}
		outfile << std::endl; 
	}
	
	outfile << "\n\n";
	if (nel > 0) {
		int i = (one ? nel : nel/2);
		outfile << std::setw(12) <<"HOMO:";
		outfile << std::setw(12) <<  i;
		outfile << std::setw(15) << eps(i-1)*TOEV << " eV\n";  
		if (i < eps.size()) { 	 
			outfile << std::setw(12) << "LUMO:";
			outfile << std::setw(12) << i+1;
			outfile << std::setw(15) << eps(i)*TOEV << " eV\n";
		} else {
			outfile << "MINIMAL BASIS USED, SO NO LUMO" << std::endl; 
		}
	}	
} 

void Logger::frequencies(const Vector& freqs, const Matrix& modes, bool printmodes) {
	int nrows = modes.rows(); 
	
	outfile << "\nFREQUENCIES (in wavenumbers)\n\n"; 
	outfile << std::setw(10) << "Mode" << std::setw(20) << "Frequency" << std::endl; 
	for (int i = 0; i < nrows; i++) {
		outfile << std::setw(10) << i+1;
		outfile << std::setw(20) << std::setprecision(6) << freqs[i] << std::endl; 
	}
		
	if (printmodes) {
		outfile << "\nNORMAL MODES\n\n"; 
		
		int ntrips = nrows / 3; 
		int i = 0; 
		for (int row = 0; row < ntrips; row++) {
			outfile << std::setw(20) << "Coordinate"; 
			outfile << std::setw(20) << i+1; 
			outfile << std::setw(20) << i+2;
			outfile << std::setw(20) << i+3 << std::endl;
			
			int j = 0; 
			for (int atom = 0; atom < ntrips; atom++) {
				outfile << std::setw(20) << "X" + std::to_string(atom+1); 
				outfile << std::setw(20) << std::setprecision(8) << modes(j, i); 
				outfile << std::setw(20) << modes(j, i+1); 
				outfile << std::setw(20) << modes(j, i+2); 
				outfile << std::endl; 
				
				outfile << std::setw(20) << "Y" + std::to_string(atom+1); 
				outfile << std::setw(20) << std::setprecision(8) << modes(j+1, i); 
				outfile << std::setw(20) << modes(j+1, i+1); 
				outfile << std::setw(20) << modes(j+1, i+2); 
				outfile << std::endl; 
				
				outfile << std::setw(20) << "Z" + std::to_string(atom+1); 
				outfile << std::setw(20) << std::setprecision(8) << modes(j+2, i); 
				outfile << std::setw(20) << modes(j+2, i+1); 
				outfile << std::setw(20) << modes(j+2, i+2); 
				outfile << std::endl; 
					
				j+=3;
			}
			
			i += 3;  
			outfile << std::endl << std::endl;  
		}
		
	}
}

// Optimisation printing routines

void Logger::initIterationOpt() {
	
	outfile << "\n";
	outfile << std::setw(12) << "Iteration";
	outfile << std::setw(24) << "Energy";
	outfile << std::setw(24) << "Delta E";
	outfile << std::setw(24) << "Grad. Norm";
	outfile << std::setw(24) << "Step. Norm"; 
	outfile << std::setw(20) << "Time elapsed\n";
	outfile << std::string(131, '-') << "\n";
	
	init_optfile(); 
	
	optfile << "\nSTARTING OPTIMIZATION\n";
	
}

void Logger::optg_dump(int iter, Vector& grad, Vector& step, SharedMolecule m, Matrix& hessian,
				double trust, double delta_e, double grad_norm, double step_norm, 
			double energy, double expect) 
{

	outfile << std::setw(12) << iter-1;
	outfile << std::setw(24) << std::setprecision(12) << energy;
	outfile << std::setw(24) << delta_e;
	outfile << std::setw(24) << grad_norm;
	outfile << std::setw(24) << step_norm; 
	outfile << std::setw(20) << std::setprecision(6) << getLocalTime();
	outfile << "\n";	
	
	flush(); 
	
	optfile << "\nITERATION " << iter-1 << std::endl << std::endl; 
	
	optfile << std::setw(10) << "Coord."; 
	optfile << std::setw(20) << "Current"; 
	optfile << std::setw(20) << "Next"; 
	optfile << std::setw(20) << "Gradient"; 
	optfile << std::setw(20) << "HNorm" << std::endl; 
	 
	std::vector<int> activeAtoms = m->getActiveList();
	int ctr = 0;
	for (int i : activeAtoms) {
		Atom& a = m->getAtom(i); 
		
		optfile << std::setw(10) << "X" << i; 
		optfile << std::setw(20) << a.getX() - step[ctr];
		optfile << std::setw(20) << a.getX(); 
		optfile << std::setw(20) << grad[ctr]; 
		optfile << std::setw(20) << hessian.row(ctr++).norm() << std::endl; 
		
		optfile << std::setw(10) << "Y" << i; 
		optfile << std::setw(20) << a.getY() - step[ctr];
		optfile << std::setw(20) << a.getY(); 
		optfile << std::setw(20) << grad[ctr]; 
		optfile << std::setw(20) << hessian.row(ctr++).norm() << std::endl;
		
		optfile << std::setw(10) << "Z" << i; 
		optfile << std::setw(20) << a.getZ() - step[ctr];
		optfile << std::setw(20) << a.getZ(); 
		optfile << std::setw(20) << grad[ctr]; 
		optfile << std::setw(20) << hessian.row(ctr++).norm() << std::endl;  
	}
	
	optfile << std::endl << "Expected change in energy (Ha): " << -expect << std::endl;
	optfile << std::endl << "Actual change in energy (Ha): " << delta_e << std::endl;
	optfile << std::endl << "Gradient norm: " << grad_norm << std::endl;
	optfile << std::endl << "Step norm: " << step_norm << std::endl;
	optfile << std::endl << "Current trust radius: " << trust << std::endl; 
 				
}
void Logger::optx_dump(int iter, Vector& grad, Vector& step, SharedMolecule m, Matrix& hessian,
		    double trust, double delta_e, double grad_norm, double step_norm, 
		  	double energy, double expect, std::vector<int>& activex) 
{
	outfile << std::setw(12) << iter;
	outfile << std::setw(24) << std::setprecision(12) << energy;
	outfile << std::setw(24) << delta_e;
	outfile << std::setw(24) << grad_norm;
	outfile << std::setw(24) << step_norm; 
	outfile << std::setw(20) << std::setprecision(6) << getLocalTime();
	outfile << "\n";	
	
	flush(); 
	
	optfile << "\nITERATION " << iter << std::endl << std::endl; 
	
	optfile << std::setw(10) << "Exponent"; 
	optfile << std::setw(20) << "Current"; 
	optfile << std::setw(20) << "Next"; 
	optfile << std::setw(20) << "Gradient"; 
	optfile << std::setw(20) << "HNorm" << std::endl; 
	 
	int nexps = m->getBasis().getNExps(); 
	double currexp; 
	int ctr = 0; 
	for (int i : activex) {
		currexp = m->getBasis().getExp(i); 
		
		optfile << std::setw(10) << i;  
		optfile << std::setw(20) << currexp - step[ctr];
		optfile << std::setw(20) << currexp; 
		optfile << std::setw(20) << grad[ctr]; 
		optfile << std::setw(20) << hessian.row(ctr++).norm() << std::endl; 
  
	}
	
	optfile << std::endl << "Expected change in energy (Ha): " << -expect << std::endl;
	optfile << std::endl << "Actual change in energy (Ha): " << delta_e << std::endl;
	optfile << std::endl << "Gradient norm: " << grad_norm << std::endl;
	optfile << std::endl << "Step norm: " << step_norm << std::endl;
	optfile << std::endl << "Current trust radius: " << trust << std::endl; 				
}

void Logger::mo_map(Vector& coeffs, SharedMolecule m, int fineness, std::string& filename) {
	
	std::ofstream mofile;
	mofile.open(filename); 
	
	if (mofile.is_open()) {
		
		double maxx = 0.0;
		double maxy = 0.0;
		double minx = 0.0;
		double miny = 0.0;
		
		double currx, curry;
		for (int i = 0; i < m->getNAtoms(); i++) {
			currx = m->getAtom(i).getX();
			curry = m->getAtom(i).getY();
			maxx = currx > maxx ? currx : maxx;
			maxy = curry > maxy ? curry : maxy;
			minx = currx < minx ? currx : minx;
			miny = curry < miny ? curry : miny; 
		}
		
		miny -= 0.5;
		minx -= 0.5;
		maxx += 0.5;
		maxy += 0.5; 
		
		double xgap = (maxx - minx) / ((double) fineness); 
		double ygap = (maxy - miny) / ((double) fineness);
		
		std::vector<libint2::Shell>& shells = m->getBasis().getIntShells(); 
		
		currx = minx;
		double x, y, z, r, r2; 
		int l; 
		for (int i = 0 ; i <= fineness; i++) {
			
			curry = miny;
			for (int j = 0; j <= fineness; j++) {
				
				double value = 0.0; 
				int ctr = 0; 
				for (auto& s : shells) {
					x = currx - s.O[0];
					y = curry - s.O[1];
					z = s.O[2]; 
					
					r2 = x*x + y*y + z*z;
					r = std::sqrt(r2); 
					
					std::vector<libint2::real_t>& exps = s.alpha; 
					std::vector<libint2::Shell::Contraction>& contr = s.contr;
					for (auto& c : contr) {
						l = c.l; 
						double tempval = 0.0;
						for (int k = 0; k < c.coeff.size(); k++)
							tempval += c.coeff[k] * std::exp(-exps[k] * r2); 
						tempval *= std::pow(r, l); 
						
						double cost = r > 0 ? z / r : 0.0; 
						double phi = atan2(y, x); 
						std::vector<double> fac = facArray(2*l + 1);
						std::vector<double> dfac = dfacArray(2*l + 1);
						TwoIndex<double> spherharms = realSphericalHarmonics(l, cost, phi, fac, dfac); 
						
						for (int m = -l; m <= l; m++)
							value += coeffs[ctr++] * tempval * spherharms(l, m+l); 
					}
				}
				
				mofile << std::setw(15) << std::setprecision(7) << currx; 
				mofile << std::setw(15) << curry; 
				mofile << std::setw(15) << value << std::endl; 
				
				curry += ygap; 
			}
			
			currx += xgap; 
		}
		
		
		
		
	} else {
		Error e1("FILEIO", "Could not open MO mapping file.");
		error(e1); 
	}
	
}

// Print orbital domains in local DFJK
void Logger::printDomains(Fock& f) {
	std::vector<Domain>& lmos = f.getLMODomains();
	std::vector<Domain>& aos = f.getAODomains(); 
	std::vector<Domain>& fits = f.getFitDomains(); 
	
	int nmos = lmos.size();
	assert(aos.size() == nmos && fits.size() == nmos); 
	
	outfile << "LOCAL DFJK DOMAINS" << std::endl << std::endl; 
	
	for (int i = 0; i < nmos; i++) {
		
		outfile << "ORBITAL " << i+1 << std::endl; 
		outfile << "LMO Domain (" << std::accumulate(lmos[i].sizes.begin(), lmos[i].sizes.end(), 0) << " functions): "; 
		for (auto d : lmos[i].centres)
			outfile << d+1 << ", ";
		outfile << std::endl; 
		
		outfile << "AO Domain (" << std::accumulate(aos[i].sizes.begin(), aos[i].sizes.end(), 0) << " functions): ";
		for (auto d : aos[i].centres)
			outfile << d+1 << ", ";
		outfile << std::endl;
		
		outfile << "DF Domain (" << std::accumulate(fits[i].sizes.begin(), fits[i].sizes.end(), 0) << " functions): ";
		for (auto d : fits[i].centres)
			outfile << d+1 << ", ";
		outfile << std::endl << std::endl;
	}
}


// Timing functions

// Print time elapsed since last call
void Logger::localTime()
{
	auto temp = std::chrono::steady_clock::now(); 
	outfile << "Time taken: " <<  std::chrono::duration<double, std::deci>(temp - last_time).count()/10.0
		<< " seconds\n";
	last_time = temp;
}
 
// Print total time taken
void Logger::globalTime()
{
	auto temp = std::chrono::steady_clock::now() - first_time; 
	outfile << "Total time: " << std::chrono::duration<double, std::deci>(temp).count()/10.0
		<< " seconds\n";
}

// Print time at which error occured
void Logger::errTime()
{
	auto temp = std::chrono::steady_clock::now() - first_time; 
	errstream << "Error after " << std::chrono::duration<double, std::deci>(temp).count()/10.0
		<< " seconds\n";
}

// Return the localTime/globalTime, instead of printing
double Logger::getLocalTime() 
{
	auto temp = std::chrono::steady_clock::now(); 
	auto rval = temp - last_time;
	last_time = temp;
	return std::chrono::duration<double, std::deci>(rval).count()/10.0; 
}

double Logger::getGlobalTime()
{
	auto temp = std::chrono::steady_clock::now() - first_time;
	return std::chrono::duration<double, std::deci>(temp).count()/10.0; 
}

// Flush the output streams
void Logger::flush() {
	outfile.flush();
	errstream.flush();
}

// Initialise the output
void Logger::init()
{
	outfile << "MOLECULAR 2015  (alpha version)\n";
	outfile << "A suite of ab initio quantum chemistry programs\n";
	std::time_t t = time(0);
	struct std::tm* now = std::localtime(&t);
	char buf[80];
	std::strftime(buf, sizeof(buf), "%d-%m-%Y %X", now);
	outfile << "Program called at date/time: " << buf << "\n";
}
	  
// Finalise the output
// Currently just prints the time and the number of errors
// that occurred
void Logger::finalise()
{
	outfile << std::string(30, '-') << "\n";
	globalTime(); // Print time elapsed
	outfile << "Number of errors: " << nerr << "\n";
}
  
