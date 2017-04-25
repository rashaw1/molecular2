#ifndef PROGRAM_CONTROLLER_HEAD
#define PROGRAM_CONTROLLER_HEAD

#include <vector>
#include <fstream>
#include <string>
#include <iostream>
#include <functional>
#include <sstream>
#include "logger.hpp"
#include "molecule.hpp"
#include "fock.hpp"
#include "scf.hpp"
#include "mp2.hpp"
#include "integrals.hpp"

struct Option {	
	Option() { }
	Option(std::string& line);
	Option(std::string _name, std::string val); 
	Option(const Option& other); 
	
	std::string name; 
	std::string _value; 
	
	template <typename T>
	T get_value() {
		std::stringstream ss(_value);
		T convertedValue;
		if ( ss >> convertedValue ) return convertedValue;
		else throw std::runtime_error("conversion failed");
	}
	
	template <typename T>
	void set_value(T val) {
		std::ostringstream os;
		os << val; 
		_value = os.str(); 
	}
	
	void parse(std::string& line); 
};

class Command {
private:
	std::vector<Option> options; 
	
public:
	Command() : molecule_id(-1), name("default") { } 
	Command(std::vector<std::string>& lines, int molecule_id, std::string name);
	Command(const Command& other);  
	
	virtual void parse(std::vector<std::string>& lines);
	
	template <typename T>
	T get_option(std::string name) {
		T value = T();
		for (auto& op : options) {
			if (op.name == name) {
				value = op.get_value<T>();
				break;
			}
		}
		return value;  
	}
	
	template <typename T>
	void set_option(std::string name, T val) {
		bool found = false;
		for (auto& op : options) {
			if (op.name == name) {
				op.set_value<T>(val);
				found = true;
				break;
			}
		}
		if (!found) {
			std::ostringstream os;
			os << val; 
			options.push_back(Option(name, os.str()));
		} 
	} 
	
	bool is_option_set(std::string name) {
		bool found = false;
		for (auto& op :options) {
			if (op.name == name) {
				found = true;
				break;
			}
		}
		return found; 
	}
	
	int molecule_id; 
	std::string name;  
};

struct Construct {
	std::vector<Construct> subconstructs;
	std::vector<std::string> content;  
	std::string name; 
	
	Construct() : name("default") { }
	Construct(std::vector<std::string>& lines, std::string& name); 
	Construct(const Construct& other);
	
	virtual void parse(std::vector<std::string>& lines);
};

class ProgramController : public std::enable_shared_from_this<ProgramController> {
private:
	std::vector<Option> global_options;
	std::vector<Command> commands; 
	std::vector<Construct> constructs;  
	std::map<std::string, std::function<void(Command&, SharedMolecule)>> command_list;
	
	std::shared_ptr<Fock> focker;
	std::shared_ptr<SCF> hf; 
	std::shared_ptr<MP2> mp2obj; 
	std::shared_ptr<IntegralEngine> ints; 
	
	std::string jobname; 
	
	bool done_hf;
	bool done_transform; 
	
	void runmp2(MP2& mp2obj, SCF& hf, bool calc); 
	
public:
	
	Logger log; 
	ProgramController(std::ifstream& input, std::ofstream& output, std::ostream& err, std::string& name); 
	ProgramController& operator=(const ProgramController& other);
	
	std::string& getName() { return jobname; }
	
	template <typename T>
	T get_option(std::string name) {
		T value = T(); 
		for (auto& op : global_options) {
			if (op.name == name) {
				value = op.get_value<T>();
				break;
			}
		}
		return value;
	}
	
	template <typename T>
	void set_option(std::string name, T val) {
		bool found = false;
		for (auto& op : global_options) {
			if (op.name == name) {
				op.set_value<T>(val); 
				found = true;
				break;
			}
		}
		if (!found) {
			std::ostringstream os;
			os << val; 
			global_options.push_back(Option(name, os.str()));
		}
	}
	
	bool is_option_set(std::string name) {
		bool found = false;
		for (auto& op : global_options) {
			if (op.name == name) {
				found = true;
				break;
			}
		}
		return found; 
	}
	 
	void parse(std::ifstream& input); 
	void run(); 
	void cleanLine(std::string& line); 
	
	void call_hf(Command& c, SharedMolecule m);
	void call_rhf(Command& c, SharedMolecule m);
	void call_uhf(Command& c, SharedMolecule m);
	void call_mp2(Command& c, SharedMolecule m);
	void call_ccsd(Command& c, SharedMolecule m);
	void call_ralmo(Command& c, SharedMolecule m);
	void call_ualmo(Command& c, SharedMolecule m);
	void call_optg(Command& c, SharedMolecule m);
	void call_optx(Command& c, SharedMolecule m); 
	void call_nuctest(Command& c, SharedMolecule m);
	
};

using SharedPC = std::shared_ptr<ProgramController>; 

#endif
