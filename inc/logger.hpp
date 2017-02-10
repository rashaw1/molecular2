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
 *                    timer - a boost::timer::cpu_timer for keeping track of time elapsed, and the time
 *                            that the log was instantiated at
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
#include <boost/timer/timer.hpp>
#include <string>
#include <fstream>
#include <iostream>
#include "basis.hpp"
#include "atom.hpp"
#include "molecule.hpp"
#include "eigen_wrapper.hpp"
#include <vector>

// Declare forward dependencies
class BF;
class PBF;
class Error;

// Begin class declaration
class Logger
{
private:
  ProgramController& control; 
  std::ofstream& outfile, intfile;
  std::ostream& errstream;
  Error* errs;
  int nerr;
  boost::timer::cpu_timer timer;
  boost::timer::nanosecond_type last_time;
  
public:
  // Conversion factors
  static const double RTOCM;
  static const double RTOGHZ;
  static const double TOKCAL;
  static const double TOKJ;
  static const double TOEV;
  static const double TOBOHR;
  static const double TOANG;
  
  // Constructor/destructor
  Logger(ProgramController &control, std::ofstream& out, std::ostream& e);
  ~Logger(); // Delete the various arrays
  // Accessors
  Logger& operator=(const Logger& other);
  
  std::ofstream& getIntFile() { return intfile; }
  
  // Overloaded print functions
  void print(const std::string& msg) const; // Print a string message
  // Print out a vector with precision digits, either horizontally or vertically
  void print(const Vector& v, int digits = 6, bool vertical = false) const; 
  void print(const Matrix& m, int digits = 6) const; // Matrix with precision digits
  void print(Basis& b, bool full = false) const; // Basis set - spec., no. of bfs, etc.
  void print(const Atom& a) const; // Atom - i.e coords, etc.
  // Print out the details of the molecule, including inertial data if wanted.
  void print(Molecule& mol, bool inertia = false) const; 
  void print(BF& bf) const; // Basis function - coeffs and each pbf
  void print(const PBF& pbf) const; // Primitive gaussian - exponent, norm, ang. momenta
  // Print out an iteration
  void iteration(int iter, double energy, double delta, double dd);
  void iterationCC(int iter, double energy, double delta_e, double delta_s, double delta_d, double t_interm, double t_amps, double t_iter);
  void initIteration();
  void initIterationCC();
  void orbitals(const Vector& eps, int nel, bool one = false);
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

