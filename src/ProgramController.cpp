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

Option::Option(std::string _name, std::string val) : name(_name), _value(val) {} 

Option::Option(const Option& other) {
	name = other.name; 
	_value = other._value; 
} 

void Option::parse(std::string& line) {
	size_t pos = line.find(','); 
	if (pos != std::string::npos) {
		name = line.substr(0, pos); 
		_value = line.substr(pos+1, line.length());
		if (_value == "true" || _value == "yes" || _value == "on") _value = std::to_string(true);
		else if (_value == "false" || _value == "no" || _value == "off") _value = std::to_string(false);
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

ProgramController::ProgramController(std::ifstream& input, std::ofstream& output, std::ostream& err) : log(*this, output, err) {
	log.init();
	log.flush();
	
	using namespace std::placeholders;
	command_list["hf"] = std::bind(&ProgramController::call_hf, this, _1, _2); 
	command_list["rhf"] = std::bind(&ProgramController::call_rhf, this, _1, _2);  
	command_list["uhf"] = std::bind(&ProgramController::call_uhf, this, _1, _2); 
	command_list["mp2"] = std::bind(&ProgramController::call_mp2, this, _1, _2); 
	command_list["ccsd"] = std::bind(&ProgramController::call_ccsd, this, _1, _2); 
	command_list["ralmo"] = std::bind(&ProgramController::call_ralmo, this, _1, _2); 
	command_list["ualmo"] = std::bind(&ProgramController::call_ualmo, this, _1, _2); 
	command_list["optg"] = std::bind(&ProgramController::call_optg, this, _1, _2); 
	
	parse(input); 
	
	log.init_intfile();
	
	if (!is_option_set("memory")) set_option<double>("memory", 100.0);
	if (!is_option_set("nthreads")) set_option<int>("nthreads", 1);
	if (!is_option_set("printeris")) set_option<bool>("printeris", false);
	if (!is_option_set("direct")) set_option<bool>("direct", false);
	if (!is_option_set("intfile")) set_option<std::string>("intfile", "eris.ints"); 
	if (!is_option_set("thrint")) set_option<double>("thrint", 1e-12);
	if (!is_option_set("bprint")) set_option<bool>("bprint", false); 
	
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
					pos = line.find(')');
					if (pos == std::string::npos) {
						if(std::getline(input, line)) {
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
	
	if (get_option<bool>("debug")) {
		std::cout << "OPTIONS: " << std::endl;
		for (auto& op : global_options) std::cout << op.name << " " << op._value << std::endl;
		std::cout << "\nCOMMANDS: " << std::endl;
		for (auto& c : commands) std::cout << c.name << " " << c.molecule_id << std::endl; 
		std::cout << "\nCONSTRUCTS: " << std::endl; 
		for (auto& c : constructs) {
			std::cout << c.name << std::endl; 
			std::cout << "CONTENT" << std::endl;
			for (auto& line : c.content) std::cout << line << std::endl; 
			std::cout << "SUBCONSTRUCTS" << std::endl; 
			for (auto& sc : c.subconstructs) {
				std::cout << sc.name << std::endl;
				std::cout << "SUBCONTENT" << std::endl;
				for (auto& line : sc.content) std::cout << line << std::endl; 
				std::cout << std::endl;
			} 
			std::cout << std::endl;
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
	
	try {
		// Build molecules
		int curr_molecule = 0;
		for (auto& c : constructs) {
			std::shared_ptr<Molecule> m = std::make_shared<Molecule>(shared_from_this(), c);
			m->buildShellBasis();
			m->buildECPBasis();
			m->calcEnuc();
			log.print(*m, true); 
			m->updateBasisPositions(); 
			log.print(m->getBasis(), get_option<bool>("bprint")); 
			log.localTime(); 
			log.flush(); 

			ints = std::make_shared<IntegralEngine>(m); 
		
			done_hf = false;
			done_transform = false; 
		
			for (auto& cmd : commands) {
				if(cmd.molecule_id != curr_molecule) continue; 
				else {
					auto it = command_list.find(cmd.name);
					if (it != command_list.end())  {
						std::function<void(Command&, SharedMolecule)> f = it->second;
						f(cmd, m); 
					}
				}
			}
		
			curr_molecule++; 
		}
	} catch (Error e) {
		log.error(e); 
	}
	libint2::finalize(); 
	log.finalise(); 
	
}

void ProgramController::call_hf(Command& c, SharedMolecule m) {
	if (m->getMultiplicity() > 1 || m->getNel()%2 != 0)
		call_uhf(c, m);
	else 
		call_rhf(c, m); 
}

void ProgramController::call_rhf(Command& c, SharedMolecule m) {
	if(!c.is_option_set("diis")) c.set_option<bool>("diis", true); 
	if(!c.is_option_set("maxdiis")) c.set_option<int>("maxdiis", 8);
	if(!c.is_option_set("maxiter")) c.set_option<int>("maxiter", 40);
	if(!c.is_option_set("converge")) c.set_option<double>("converge", 1e-5); 
	if(!c.is_option_set("precision")) c.set_option<double>("precision", 1e-12); 
	
	focker = std::make_shared<Fock>(c, *ints, m);
	hf = std::make_shared<SCF>(c, m, *focker); 
	hf->rhf(); 
	done_hf = true;
}

void ProgramController::call_uhf(Command& c, SharedMolecule m) {
	if(!c.is_option_set("diis")) c.set_option<bool>("diis", true); 
	if(!c.is_option_set("maxdiis")) c.set_option<int>("maxdiis", 8);
	if(!c.is_option_set("maxiter")) c.set_option<int>("maxiter", 40);
	if(!c.is_option_set("converge")) c.set_option<double>("converge", 1e-5); 
	if(!c.is_option_set("precision")) c.set_option<double>("precision", 1e-12); 

	focker = std::make_shared<UnrestrictedFock>(c, *ints, m);
	hf = std::make_shared<SCF>(c, m, *focker); 
	hf->uhf(); 
	done_hf = true;
}

void ProgramController::call_mp2(Command& c, SharedMolecule m) {
	if(done_hf && !done_transform) { 
		mp2obj = std::make_shared<MP2>(*focker); 
		runmp2(*mp2obj, *hf, true); 
		done_transform = true;
	} else {
		Error e("MP2", "HF is required before MP2 can be done.");
		log.error(e);
	}
}

void ProgramController::call_ccsd(Command& c, SharedMolecule m) {
	if(!c.is_option_set("diis")) c.set_option<bool>("diis", true);
	if(!c.is_option_set("triples")) c.set_option<bool>("triples", false);
	if(!c.is_option_set("maxdiis")) c.set_option<int>("maxdiis", 5);
	if(!c.is_option_set("maxiter")) c.set_option<int>("maxiter", 30); 
	if(!c.is_option_set("converge")) c.set_option<double>("converge", 1e-5); 
	
	if (done_hf) {
		if(!done_transform) {
			mp2obj = std::make_shared<MP2>(*focker); 
			runmp2(*mp2obj, *hf, false);
			done_transform = true;
		}
		
		mp2obj->spatialToSpin();
		CCSD ccobj(c, *mp2obj); 
		ccobj.compute();
		log.result("Total Energy = ", hf->getEnergy() + ccobj.getEnergy() + ccobj.getETriples(), "Hartree");
	} 
}

void ProgramController::call_ralmo(Command& c, SharedMolecule m) {
	if(!c.is_option_set("diis")) c.set_option<bool>("diis", true);
	if(!c.is_option_set("maxdiis")) c.set_option<int>("maxdiis", 6);
	if(!c.is_option_set("converge")) c.set_option<double>("converge", 1e-5); 
	if(!c.is_option_set("maxiter")) c.set_option<int>("maxiter", 40);
	if(!c.is_option_set("perturb")) c.set_option<int>("perturb", 0);
	if(!c.is_option_set("precision")) c.set_option<double>("precision", 1e-12); 

	focker = std::make_shared<Fock>(c, *ints, m); 
	ALMOSCF almo(c, m, *focker); 
	almo.rscf(); 
}

void ProgramController::call_ualmo(Command& c, SharedMolecule m) {
	if(!c.is_option_set("diis")) c.set_option<bool>("diis", true);
	if(!c.is_option_set("maxdiis")) c.set_option<int>("maxdiis", 6);
	if(!c.is_option_set("converge")) c.set_option<double>("converge", 1e-5); 
	if(!c.is_option_set("maxiter")) c.set_option<int>("maxiter", 40);
	if(!c.is_option_set("perturb")) c.set_option<int>("perturb", 0);
	if(!c.is_option_set("precision")) c.set_option<double>("precision", 1e-12); 

	focker = std::make_shared<UnrestrictedFock>(c, *ints, m);
	ALMOSCF almo(c, m, *focker); 
	almo.uscf(); 
}

void ProgramController::call_optg(Command& c, SharedMolecule m) {
	// To be implemented
}

void ProgramController::runmp2(MP2& mp2obj, SCF& hf, bool calc) {
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

ProgramController& ProgramController::operator=(const ProgramController& other) {
	global_options = other.global_options;
	commands = other.commands;
	constructs = other.constructs;
	command_list = other.command_list;
	focker = other.focker;
	hf = other.hf;
	mp2obj = other.mp2obj;
	ints = other.ints;
	done_hf = other.done_hf;
	done_transform = other.done_transform; 
	log = other.log;
	return *this;
}
