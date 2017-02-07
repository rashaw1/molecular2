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
#include "logger.hpp"
#include "molecule.hpp"
#include "error.hpp"
#include "integrals.hpp"
#include "fock.hpp"
#include "scf.hpp"
#include "almoscf.hpp"
#include "mp2.hpp"
#include "cc.hpp"
#include "gaussquad.hpp"
#include "gshell.hpp"
#include "ecp.hpp"
#include <functional>
#include <cmath>
#include "ecpint.hpp"
#include "multiarr.hpp"
#include <libint2.hpp>
#include "optimiser.hpp"

void runmp2(MP2& mp2obj, Logger& log, SCF& hf, bool calc); 

int main (int argc, char* argv[])
{
	if (argc == 1) { 
		std::cerr << "You must provide an input file.\nUSAGE: ./molecular inputfile.mol\n";
	} else {
		// Get the input file name, and create the error and output file names
		std::string ifname = argv[1];
		std::string ofname = ifname;
		std::size_t pos = ofname.find('.');
		if (pos != std::string::npos) { 
			ofname.erase(pos, ofname.length());
		}    
		std::string efname = ofname + ".log";
		ofname += ".out";

		// Open the file streams
		std::ifstream input(ifname);
		if (!input.is_open()) {
			std::cerr << "Could not open input file.\n";
		} else {
			std::ofstream output(ofname);
			std::ofstream err(efname);
      
			// Create the logger
			Logger log(input, output, err);
			log.init();
			log.flush();

			try{
				// Create the molecule and print initial output
				Molecule mol(log, 1);
				mol.buildShellBasis();
				mol.buildECPBasis();
				mol.calcEnuc();
				log.print(mol, true);
				mol.updateBasisPositions();
				log.print(log.getBasis(), log.bprint());
				log.print("\nPRELIMINARIES FINISHED");
				log.localTime();
				log.flush();
				
				// Make an integral engine and set up Eigen
				libint2::initialize();
				Eigen::setNbThreads(mol.getLog().getNThreads());
				
				IntegralEngine ints(mol);

				// All calculations will need some form of HF
				Fock* focker;
				SCF* hf; 
								
				// Run the commands given
				int cmd = log.nextCmd();
				while (cmd!=0) {
					switch(cmd) {
						case 1: { // HF
							if(mol.getMultiplicity() > 1){
								focker = new UnrestrictedFock(ints, mol);
								hf = new SCF(mol, *focker); 
								hf->uhf();
							} else if ((mol.getNel()%2) != 0) {
								focker = new UnrestrictedFock(ints, mol);
								hf = new SCF(mol, *focker);
								hf->uhf();
							} else {
								focker = new Fock(ints, mol);
								hf = new SCF(mol, *focker); 
								hf->rhf();
							}
							break;
						}
						case 2: { // RHF
							focker = new Fock(ints, mol);
							hf = new SCF(mol, *focker);
							hf->rhf();
							break;
						}
						case 3: { // UHF
							focker = new UnrestrictedFock(ints, mol);
							hf = new SCF(mol, *focker);
							hf->uhf();
							break;
						}
						case 4: { // Just MP2
							MP2 mp2obj(*focker);
							runmp2(mp2obj, log, *hf, true);
							break;
						}
						case 5: { // MP2 + CCSD 
							MP2 mp2obj(*focker);
							runmp2(mp2obj, log, *hf, true);
							mp2obj.spatialToSpin();
							CCSD ccobj(mp2obj, false, log.diis());
							ccobj.compute();
							log.result("Total Energy = ", hf->getEnergy() + ccobj.getEnergy(), "Hartree");
							break;
						}
						case 6: { // MP2 + CCSD(T) 
							MP2 mp2obj(*focker);
							runmp2(mp2obj, log, *hf, true);
							mp2obj.spatialToSpin();
							CCSD ccobj(mp2obj, true, log.diis());
							ccobj.compute();
							log.result("Total Energy = ", hf->getEnergy() + ccobj.getEnergy() + ccobj.getETriples(), "Hartree");
							break;
						}
						case 7: { // Just CCSD
							MP2 mp2obj(*focker);
							runmp2(mp2obj, log, *hf, false);
							mp2obj.spatialToSpin();
							CCSD ccobj(mp2obj, false, log.diis());
							ccobj.compute();
							log.result("Total Energy = ", hf->getEnergy() + ccobj.getEnergy(), "Hartree");
							break;
						}
						case 8: { // Just CCSD(T) 
							MP2 mp2obj(*focker);
							runmp2(mp2obj, log, *hf, false);
							mp2obj.spatialToSpin();
							CCSD ccobj(mp2obj, true, log.diis());
							ccobj.compute();
							log.result("Total Energy = ", hf->getEnergy() + ccobj.getEnergy() + ccobj.getETriples(), "Hartree");
							break;
						}
						case 9: { // ALMO SCF 
							ALMOSCF almo(mol, *focker); 
							almo.rscf(); 
							break;
						}
						default: { }
					}
					log.flush();
					cmd = log.nextCmd();
				}
				
				// Finalise the run
				libint2::finalize();
				log.finalise();
			} catch (Error e){
				log.error(e);
			}
		   
			// Close file streams
			input.close();
			output.close();
			err.close();
		}
	}

	return 0;
}

void runmp2(MP2& mp2obj, Logger& log, SCF& hf, bool calc) {
	if (calc)
		log.title("MP2 CALCULATION");
	else
		log.title("INTEGRAL TRANSFORMATION");
	mp2obj.transformIntegrals();
	log.print("Integral transformation complete.\n");
	log.localTime();
	if (calc) {
		mp2obj.calculateEnergy();
		log.result("MP2 Energy Correction", mp2obj.getEnergy(), "Hartree");
		log.result("Total Energy = ", hf.getEnergy() + mp2obj.getEnergy(), "Hartree");
	}
}

