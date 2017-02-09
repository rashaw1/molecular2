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
	T Command::getOption(std::string name) {
		T value = T();
		for (auto& op : options)
			if (op.name == name) value = op.get_option<T>();
		return value;  
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

class ProgramController {
private:
	std::vector<Option> global_options;
	std::vector<Command> commands; 
	std::vector<Construct> constructs;  
	std::map<std::string, std::function<void(Command&)> command_list;
	std::vector<Molecule> molecules; 
	
	std::shared_ptr<Fock> focker;
	std::shared_ptr<SCF> hf; 
	std::shared_ptr<MP2> mp2obj; 
	std::shared_ptr<IntegralEngine> ints; 
	
public:
	
	Logger log; 
	ProgramController() { }
	ProgramController(std::ifstream& input, std::ofstream& output, std::ostream& err); 
	
	template <typename T>
	T get_option(std::string name) {
		T value = T(); 
		for (auto& op : global_options) 
			if (op.name == name) value = op.get_value<T>();
		return value;
	}
	
	template <typename T>
	void set_option(std::string name, T val) {
		bool found = false;
		for (auto& op : global_options) {
			if (op.name == name) {
				op._value = std::to_string(val); 
				found = true;
			}
		}
		if (!found) global_options.push_back(Option(name + "," + std::to_string(val))); 
	}	
	 
	void parse(std::ifstream& input); 
	void run(); 
	void cleanLine(std::string& line); 
	
	void call_hf(Command& c);
	void call_rhf(Command& c);
	void call_uhf(Command& c);
	void call_mp2(Command& c);
	void call_ccsd(Command& c);
	void call_almo(Command& c);
	void call_optg(Command& c);
	
};

#endif
