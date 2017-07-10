/* 
 *   
 *     PURPOSE: To define a class Logger, which will perform the majority of IO.
 *				It esentially acts as the communicator 
 *              between all the different units, 'logging' everything that happens. 
 *              For convenience, fundamental constants are stored here as well. 
 *
 *     class Logger:
 *              owns: outfile - an ofstream, for all primary logging functions to permanent
 *                              record.
 *                    infile - an ifstream for the input file
 *                    errstream - an ostream for logging all error messages - could be a file,
 *                                the console, any ostream
 *                    errs - an array of Error messages that have been thrown
 *                    intfile - the file to print two electron integrals to if wanted
 *              data: nerr - the number of errors accumulated
 *                    last_time - the last time that timer.elapsed was called
 *              fundamental constants:
 *                    M_PI - a definition of pi, in case one is not available for some reason
 *              conversion factors:
 *                    RTOCM, RTOGHZ - convert the rotation constants calculated by molecule into
 *                                    units of cm-1 or GHz, respectively.
 *                    TOKCAL, TOKJ - convert energies in Hartrees to kcal/mol or kj/mol
 *                    TOBOHR, TOANG - convert distances either from angstrom to bohr, or bohr 
 *                                    to angstrom, respectively
 *              routines:
 *                    OUTPUT:
 *                         print() - overloaded function that will just print a general message
 *                                   or object in a generic way
 *                         title() - will print a message as a title
 *                         result() - will print a message specifically formatted to stand out
 *                                    as a result
 *                         error() - will log an Error passed to it, in errs, and print it to
 *                                   the errstream.
 *                         localTime() - will print the time elapsed since last call to outfile,
 *                                       and set last_time to current time
 *                         globalTime() - will print the total time elapsed since beginning to
 *                                       outfile, without changing last_time.
 *                         getLocalTime(), getGlobalTime() - get values in secs rather than print.
 *                         errTime() - will print the total time elapsed to the errstream,
 *                                     without changing last_time.
 *                         finalise() - will close the output streams, flushing the buffers, and
 *                                      adding the customary final parts of the output.
 *
 *        DATE            AUTHOR              CHANGES    
 *        ===========================================================================
 *        27/08/15        Robert Shaw         Original code.
 *        28/08/15        Robert Shaw         Changed how timing works.
 */

#ifndef LOGGERHEADERDEF
#define LOGGERHEADERDEF

// Exceptional use of preprocessor here for fundamental constants
// purely because they aren't always defined by every compiler
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Includes
#include <string>
#include <fstream>
#include <iostream>
#include "basis.hpp"
#include "atom.hpp"
#include "molecule.hpp"
#include "eigen_wrapper.hpp"
#include "ecp.hpp"
#include <vector>
#include <chrono>
#include <libint2.hpp>

// Declare forward dependencies
class Error;
class Fock; 

// Begin class declaration
class Logger
{
private:
  ProgramController& control; 
  std::ofstream& outfile;
  std::ofstream intfile, optfile;
  std::ostream& errstream;
  Error* errs;
  int nerr;
  std::chrono::steady_clock::time_point last_time, first_time;
  
public:
  // Conversion factors
  static const double RTOCM;
  static const double RTOGHZ;
  static const double TOKCAL;
  static const double TOKJ;
  static const double TOEV;
  static const double TOBOHR;
  static const double TOANG;
  static const double TOWAVENUMBERS; 
  
  // Constructor/destructor
  Logger(ProgramController &control, std::ofstream& out, std::ostream& e);
  ~Logger(); // Delete the various arrays
  // Accessors
  Logger& operator=(const Logger& other);
  
  void init_intfile();
  void init_optfile(); 
  
  std::ofstream& getIntFile() { return intfile; }
  std::ofstream& getOptFile() { return optfile; }
  
  // Overloaded print functions
  void print(const std::string& msg) const; // Print a string message
  // Print out a vector with precision digits, either horizontally or vertically
  void print(const Vector& v, int digits = 6, bool vertical = false) const; 
  void print(const Matrix& m, int digits = 6) const; // Matrix with precision digits
  void print(Basis& b, SharedMolecule m, bool full = false) const; // Basis set - spec., no. of bfs, etc.
  void print(ECPBasis& b, bool full = false) const; // ECP Basis set - spec., no. of bfs, etc.
  void print(ECP& ecp, bool full = false) const; 
  void print(std::vector<libint2::Shell>& basis, std::vector<int>& shellAtomList, 
					SharedMolecule m, bool full = false) const; 
  void print(const Atom& a) const; // Atom - i.e coords, etc.
  // Print out the details of the molecule, including inertial data if wanted.
  void print(Molecule& mol, bool inertia = false) const; 

  // Print out an iteration
  void iteration(int iter, double energy, double delta, double dd);
  void iterationCC(int iter, double energy, double delta_e, double delta_s, double delta_d, double t_interm, double t_amps, double t_iter);
  void initIteration();
  void initIterationCC();
  void initALMOTable(); 
  void initCTTable(); 
  void ALMORow(int f1, int f2, double sep, double edisp, double edispexch, double eionic, double ebsse, double eintra); 
  void CTRow(int f1, int f2, double en); 
  void orbitals(const Vector& eps, int nel, bool one = false);
  void coefficient_matrix(const Vector& eps, int nel, const Matrix& coeffs, bool one = false); 
  void frequencies(const Vector& freqs, const Matrix& modes, bool printmodes); 
  void printDomains(Fock& f); 
  
  void initIterationOpt(); 
  void optg_dump(int iter, Vector& grad, Vector& step, SharedMolecule m, Matrix& hessian,
  				double trust, double delta_e, double grad_norm, double step_norm, 
				double energy, double expect);
  void optx_dump(int iter, Vector& grad, Vector& step, SharedMolecule m, Matrix& hessian,
			    double trust, double delta_e, double grad_norm, double step_norm, 
			  	double energy, double expect, std::vector<int>& activex);
				
  void mo_map(Vector& coeffs, SharedMolecule m, int fineness, std::string& filename);  
  
  // Specific logging formats
  void title(const std::string& msg) const;
  void result(const std::string& msg) const;
  void result(const std::string& name, const double value, const std::string& units) const;
  void error(Error& e);
  void localTime();
  void globalTime();
  void errTime();
  double getLocalTime();
  double getGlobalTime();
  void flush();
  void init();
  void finalise();
};
 
#endif

