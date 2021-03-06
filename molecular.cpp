/*
*
*   PURPOSE: The main program loop for the molecular suite of ab initio quantum
*            chemistry routines.
* 
*   DATE              AUTHOR                CHANGES
*   ===========================================================================
*   23/09/15          Robert Shaw           Original code.
*
*/

#include <iostream>
#include <fstream>
#include <string>
#include "mpi.h"
#include "ProgramController.hpp"

int main (int argc, char* argv[])
{
	MPI_Init(&argc, &argv);
	if (argc == 1) { 
		std::cerr << "You must provide an input file.\nUSAGE: ./molecular inputfile.mol\n";
	} else {
		// Get the input file name, and create the error and output file names
		std::string ifname = argv[1];
		std::string name = ifname;
		std::size_t pos = name.find('.');
		if (pos != std::string::npos) { 
			name.erase(pos, name.length());
		}    
		std::string efname = name + ".log";
		std::string ofname = name + ".out";

		// Open the file streams
		std::ifstream input(ifname);
		if (!input.is_open()) {
			std::cerr << "Could not open input file.\n";
		} else {
			std::ofstream output(ofname);
			std::ofstream err(efname);
      
			// Create the program controller
			std::shared_ptr<ProgramController> control = std::make_shared<ProgramController>(input, output, err, name);
			control->run(); 
		   
			// Close file streams
			input.close();
			output.close();
			err.close();
		}
	}
	MPI_Finalize();
	return 0;
}

