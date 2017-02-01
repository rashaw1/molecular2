/*
*
*   PURPOSE: Implements filereader.hpp, a class for reading input files.
*
*   DATE           AUTHOR           CHANGES
*   ==================================================================
*   30/08/15       Robert Shaw      Original code.
*
*/
 
#include "filereader.hpp"
#include "ioutil.hpp"
#include "error.hpp"
#include <algorithm>
#include <iostream> 

// Class FileReader implementation

// Destructor
FileReader::~FileReader()
{
}

// Match a token to a particular command
int FileReader::findToken(std::string t)
{
	// Get rid of extraneous spaces
	t.erase(std::remove(t.begin(), t.end(), ' '), t.end());
  
	// make t lowercase
	std::transform(t.begin(), t.end(), t.begin(), ::tolower);

	// Match it to a command
	int rval = 0;
	if (t == "charge") { rval = 1; }
	else if (t == "multiplicity") { rval = 2; }
	else if (t == "maxiter") { rval = 3; }
	else if (t == "precision") { rval = 4; }
	else if (t == "basis") { rval = 5; }
	else if (t == "geom") { rval = 6; }
	else if (t == "thrint") { rval = 7; }
	else if (t == "memory") { rval = 8; }
	else if (t == "integral") { rval = 9; }  
	else if (t == "print") { rval = 10; }
	else if (t == "direct") { rval = 11; }
	else if (t == "scf") { rval = 12; }
	else if (t == "nodiis") { rval = 13; }
	else if (t == "converge") { rval = 14; }
	else if (t == "hf") { rval = 15; }
	else if (t == "rhf") { rval = 16; }
	else if (t == "uhf") { rval = 17; }
	else if (t == "angstrom") { rval = 18; }
	else if (t == "nthreads") { rval = 19; }
	else if (t == "mp2") { rval = 20; }
	else if (t == "ccsd") { rval = 21; }
	else if (t == "ccsd(t)") { rval = 22; }
	else if (t == "ecp") { rval = 23; }
	else if (t == "fragments") { rval = 24; }
	else if (t == "almo") { rval = 25; }
	return rval;
}

// Read in the parameters from the input file
void FileReader::readParameters()
{
	// Set default values, in case none specified
	charge = 0;
	multiplicity = 1;
	maxiter = 50;
	precision = 1e-12;
	geomstart = 0; geomend = 0;
	thrint = 1e-12;
	converge = 1e-5;
	memory = 100;
	nthreads = 1;
	direct = false;
	twoprint = false;
	bprint = false;
	diis = true;
	angstrom = false;
	fragments = false;
	ecp = false;
	basis[0] = "sto-3g";

	// Read line by line and parse
	std::string line, token;
	size_t pos;
	int linecount = 0;
	while (std::getline(input, line)) {
		// Erase any comments
		pos = line.find('!');
		if (pos != std::string::npos){
			line.erase(pos, line.length());
		}

		// Tokenise
		pos = line.find(',');
		if(pos != std::string::npos){
			token = line.substr(0, pos);

			// Match token to command
			switch(findToken(token)){
				case 1: { // Charge
					charge = std::stoi(line.substr(pos+1, line.length()));
					break;
				}
				case 2: { // Multiplicity
					multiplicity = std::stoi(line.substr(pos+1, line.length()));
					break;
				}
				case 3: { // MaxIter
					maxiter = std::stoi(line.substr(pos+1, line.length()));
					break;
				}
				case 4: { // Precision
					precision = std::stod(line.substr(pos+1, line.length()));
					break;
				}
				case 5: { // Basis
					line.erase(0, pos+1);
					// Get rid of excess spaces
					pos = line.find(',');
					if (pos != std::string::npos) {
						token = line.substr(0, pos);
						token.erase(std::remove(token.begin(), token.end(), ' '), token.end());
						std::transform(token.begin(), token.end(), token.begin(), ::tolower);
						std::string token2 = line.substr(pos+1, line.length());
						token2.erase(std::remove(token2.begin(), token2.end(), ' '), token2.end());
						if (token == "default")
							basis[0] = token2;
						else {
							int q = getAtomCharge(token);
							basis[q] = token2;
						}	
					}
					break;
				}
				case 6: { // GeomStart
					geomstart = linecount + 1;
					line.erase(0, pos+1);
					line.erase(std::remove(line.begin(), line.end(), ' '), line.end());
					if (line.length() > 0) {
						switch(findToken(line)){
							case 18: { // Angstrom
								angstrom = true;
								break;
							}
							default: {
							}
						}
					}
					geomend = geomstart-1;
					while(line != "geomend"){
						if (std::getline(input, line))
							geomend++;
					}
					break;
				}
				case 7: { // Thrint
					thrint = std::stod(line.substr(pos+1, line.length()));
					break;
				}
				case 8: { // Memory
					memory = std::stod(line.substr(pos+1, line.length()));
					break;
				}
				case 9: { // Integral directive
					line.erase(0, pos+1);
					// Tokenise the next bit
					pos = line.find(',');
					if(pos != std::string::npos){
						token = line.substr(0, pos);
						switch(findToken(token)){
							case 10: { // Print the integrals, file specified
								twoprint = true;
								line.erase(0, pos+1);
								// Get rid of excess spaces
								line.erase(std::remove(line.begin(), line.end(), ' '), line.end());
								intfile  = line;
								break;
							}
							default: { 
								throw(Error("READIN", "Command " + token + " not found."));
							}
						}
					} else {
						// Get rid of excess spaces
						line.erase(std::remove(line.begin(), line.end(), ' '), line.end());
						switch(findToken(line)){
							case 10: { // Print the integrals to default file name
								twoprint = true;
								intfile = "twoints.out";
								break;
							}
							case 11: { // Direct
								direct = true;
								break;
							}
							case 5: { // print basis details
								bprint = true;
								break;
							}
							default: {
								throw(Error("READIN", "Command " + line + " not found."));
							}
						}
					}
					break;
				}
				case 12: { // SCF directive
					line.erase(0, pos+1);
					// Tokenise the next bit
					pos = line.find(',');
					if (pos != std::string::npos) {
						token = line.substr(0, pos);
						switch(findToken(token)){
							case 14: { // Convergence criterion specified
								converge = std::stod(line.substr(pos+1, line.length()));
								break;
							}
							default: {
								throw(Error("READIN", "Command " + token + " not found."));
							}
						}
					} else {
						line.erase(std::remove(line.begin(), line.end(), ' '), line.end());
						switch(findToken(line)){
							case 13: { // No DIIS wanted
								diis = false;
								break;
							}
							default: {
								throw(Error("READIN", "Command " + line + " not found."));
							}
						}
					}
					break;
				}
				case 15: { // HF directive
					commands.push_back("HF");
					break;
				}
				case 16: { // RHF directive
					commands.push_back("RHF");
					break;
				}
				case 17: { // UHF directive
					commands.push_back("UHF");
					break;
				}
				case 19: { // Nthreads
					nthreads = std::stoi(line.substr(pos+1, line.length()));
					break;
				}
				case 20: { // MP2 directive
					commands.push_back("MP2");
					break;
				}
				case 21: { // CCSD directive 
					commands.push_back("CCSD");
					break;
				}
				case 22: { // CCSD(T) directive
					commands.push_back("CCSD(T)");
					break;
				}
				case 23: { // ECP basis 
					line.erase(0, pos+1);
					ecp = true;
					// Get rid of excess spaces
					pos = line.find(',');
					if (pos != std::string::npos) {
						token = line.substr(0, pos);
						token.erase(std::remove(token.begin(), token.end(), ' '), token.end());
						std::transform(token.begin(), token.end(), token.begin(), ::tolower);
						std::string token2 = line.substr(pos+1, line.length());
						token2.erase(std::remove(token2.begin(), token2.end(), ' '), token2.end());
						int q = getAtomCharge(token);
						basis[-q] = token2;
					}
					break;
				}
				case 24: { // Fragment definitions
					fragments = true;
					std::getline(input, line);
					line.erase(std::remove(line.begin(), line.end(), ' '), line.end());
					while(line != "fragmentsend") {
						pos = line.find(',');
						std::vector<int> f;
						while (pos != std::string::npos) {
							token = line.substr(0, pos);
							f.push_back(std::stoi(token));
							line.erase(0, pos+1);
							pos = line.find(',');
						}
						f.push_back(std::stoi(line));
						if (f.size() > 0) {
							f[0] -= 1;
							frags.push_back(f);
						}
						std::getline(input, line);
						line.erase(std::remove(line.begin(), line.end(), ' '), line.end());	
					}
					break;
				}
				case 25: { //ALMO SCF
					commands.push_back("ALMO"); 
					break;
				}
				default: { // Unkown command issued
					throw(Error("READIN", "Command " + token + " not found."));
				}
			}
		}
		linecount++;
	}
	natoms = geomend - geomstart;
}

void FileReader::readGeometry()
{
	// The array size will be geomend-geomstart
	if (natoms != 0){
		// Rewind to beginning of file
		input.clear();
		input.seekg(0, std::ios::beg);
		std::string line;
		// Copy the geometry into the array
		for (int l = 0; l < geomend; l++){
			std::getline(input, line);
			if (l > geomstart - 1){
				geometry.push_back(line);
			}
		}
	} else {
		throw(Error("READGM", "Geometry not found.")); 
	}
}
