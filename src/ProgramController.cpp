#include "ProgramController.hpp"
#include "error.hpp"
#include "almoscf.hpp"
#include "cc.hpp"
#include "optimiser.hpp"
#include <libint2.hpp>
#include "eigen_wrapper.hpp"

#include <algorithm>

Option::Option(std::string& line) {
	parse(line); 
}

Option::Option(const Option& other) {
	name = other.name; 
	_value = other._value; 
} 

void Option::parse(std::string& line) {
	size_t pos = line.find(','); 
	if (pos != std::string::npos) {
		name = line.substr(0, pos); 
		_value = line.substr(pos+1, line.length());
	}
}

Command::Command(std::vector<std::string>& lines, int _molecule_id, std::string _name) : molecule_id(_molecule_id), name(_name) {
	parse(lines); 
}

Command::Command(const Command& other) {
	options = other.options; 
	molecule_id = other.molecule_id; 
	name = other.name; 
}

void Command::parse(std::vector<std::string>& lines) {
	for (auto& line : lines)
		options.push_back(Option(line)); 
}

Construct::Construct(std::vector<std::string>& lines, std::string& _name) : name(_name) {
	parse(lines); 
}

Construct::Construct(const Construct& other) {
	subconstructs = other.subconstructs;
	content = other.content; 
	name = other.name; 
}

void Construct::parse(std::vector<std::string>& lines) {
	size_t pos; 
	std::string token; 
	for (int i = 0; i < lines.size(); i++) {
		pos = lines[i].find('=');
		if (pos != std::string::npos) { // Subconstruct
			token = lines[i].substr(0, pos); 
			
			std::vector<std::string> sublines; 
			pos = lines[++i].find('}'); 
			bool end_found = (pos != std::string::npos);
			int nsubs = 0; 
			
			while(i < lines.size() && !end_found) {
				sublines.push_back(lines[i]);
				pos = lines[i].find('=');
				if (pos != std::string::npos) {
					nsubs++;
					i++;  
				} else {
					pos = lines[i++].find('}');
					if (pos != std::string::npos && nsubs == 0) {
						end_found = true;
						sublines.pop_back();
						i--;
					} else if (pos != std::string::npos) nsubs--; 
				}
			} 
			
			subconstructs.push_back(Construct(sublines, token)); 
		} else content.push_back(lines[i]); 
	}
}

ProgramController::ProgramController(std::ifstream& input, std::ofstream& output, std::ostream& err) : log(output, err) {
	log.init();
	log.flush();
	
	command_list["hf"] = call_hf;
	command_list["rhf"] = call_rhf; 
	command_list["uhf"] = call_uhf;
	command_list["mp2"] = call_mp2;
	command_list["ccsd"] = call_ccsd;
	command_list["almo"] = call_almo;
	command_list["optg"] = call_optg;
	
	parse(input); 
}

void ProgramController::cleanLine(std::string& line) {
	size_t pos = line.find('!');
	if (pos != std::string::npos){
		line.erase(pos, line.length());
	}
	
	line.erase(std::remove(line.begin(), line.end(), ' '), line.end()); 
	line.erase(std::remove(line.begin(), line.end(), '\t'), line.end()); 
	std::transform(line.begin(), line.end(), line.begin(), ::tolower); 
}

void ProgramController::parse(std::ifstream& input) {
	std::string line, token;
	size_t pos;
	int curr_molecule = -1; 
	
	while(std::getline(input, line)) {
		// Erase any comments
		cleanLine(line);
		
		// Tokenise
		pos = line.find('=');
		if (pos != std::string::npos) { // Construct
			token = line.substr(0, pos);
			std::vector<std::string> sublines; 
			if (std::getline(input, line)) {
				pos = line.find('}');
				bool end_found = (pos != std::string::npos);
				int nsubs = 0;  
			
				while(!end_found) {
					cleanLine(line);
					sublines.push_back(line);
					pos = line.find('=');
					if (pos != std::string::npos) { 
						nsubs++; 
						if (!std::getline(input, line)) end_found = true;
					} else {
						pos = line.find('}');
						if (pos != std::string::npos && nsubs == 0) {
							end_found = true;
							sublines.pop_back(); 
						} else if (pos != std::string::npos) {
							nsubs--;
							if (!std::getline(input, line)) end_found = true;
						}  else if (!std::getline(input, line)) end_found = true; 
					}
				}  
			}
			constructs.push_back(Construct(sublines, token));
			curr_molecule++; 
		} else {
			pos = line.find('('); 
			if (pos != std::string::npos) { // Command
				token = line.substr(0, pos);   
			 
				auto it = command_list.find(token);
				if (it == command_list.end()) {
					Error e("IO", "Could not find command " + token);
					log.error(e);
				} else { // Read in options
					bool options_end = false; 
					std::vector<std::string> lines; 
					if (std::getline(input, line)) {
						while(!options_end) {
							pos = line.find(')'); 
							if (pos != std::string::npos) options_end = true; 
							else {
								cleanLine(line); 
								lines.push_back(line);
								if (!std::getline(input, line)) options_end = true; 
							}
						}
					}
				
					commands.push_back(Command(lines, curr_molecule, token)); 
				}
			} else {
			
				pos = line.find(','); 
				if (pos != std::string::npos) // Global option
					global_options.push_back(Option(line)); 
			
			}	
		}
	}
}

void ProgramController::run() {
	
	libint2::initialize(); 
	Eigen::setNbThreads(get_option<int>("nthreads")); 
	
	focker = nullptr;
	hf = nullptr;
	mp2obj = nullptr;
	ints = nullptr; 
	
	// Build molecules
	for (auto& c : constructs)
		molecules.push_back(Molecule(*this, c));
	
	int curr_molecule = 0; 
	for (auto& m : molecules) {
		m.buildShellBasis();
		m.buildECPBasis();
		m.calcEnuc();
		log.print(m, true); 
		m.updateBasisPositions(); 
		log.print(m.getBasis(), get_option<bool>("bprint")); 
		log.localTime(); 
		log.flush(); 
		
		ints = new IntegralEngine(m); 
		
		for (auto& c : commands) {
			if(c.molecule_id != curr_molecule) continue; 
			else {
				auto it = command_list.find(c.name);
				if (it != command_list.end())  {
					std::function<void(Command&)> f = it->second;
					f(c); 
				}
			}
		}
		
		curr_molecule++; 
	}
	
	libint2::finalize();  
	
}

void ProgramController::call_hf(Command& c) {
	Molecule& m = molecules[c.molecule_id];
	if (m.getMultiplicity() > 1 || mol.getNel()%2 != 0)
		call_uhf(c);
	else 
		call_rhf(c); 
}

void ProgramController::call_rhf(Command& c) {
	Molecule& m = molecules[c.molecule_id]; 
	focker = new Fock(c, ints, m);
	hf = new SCF(m, *focker);
	hf->rhf(); 
}

void ProgramController::call_uhf(Command& c) {
	Molecule& m = molecules[c.molecule_id];
	focker = new UnrestrictedFock(c, ints, m);
	hf = new SCF(m, *focker);
	hf->uhf(); 
}
