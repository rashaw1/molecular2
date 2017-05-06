/*
*
*   PURPOSE: Implements basisreader.hpp, a class for reading basis files.
*
*   DATE           AUTHOR           CHANGES
*   ==================================================================
*   30/08/15       Robert Shaw      Original code.
*   04/09/15       Robert Shaw      Now indexes primitives.
*/
 
#include "basisreader.hpp"
#include "ioutil.hpp"
#include "bf.hpp"
#include "ecp.hpp"
#include "error.hpp"
#include "basis.hpp"
#include <algorithm>
#include <iostream>
 
// Implement class BasisReader

// Open the basis file that contains the correct basis functions
void BasisReader::openFile(int q)
{
	// Find which row of elements q is in
	// Note that the first row is taken here to be H-Ne
	std::string row = "first";
	int Q = q < 0 ? -q : q; 
	if (Q > 10 && Q < 19){ row = "second"; } // Na - Ar
	else if (Q > 18 && Q < 37) { row = "third"; } // K - Kr
	else if (Q > 36 && Q < 55) { row = "fourth"; } // Rb - Xe
	else if (Q > 54 && Q < 87) { row = "fifth"; } // Cs - Rn
	else if (Q > 86) { row = "sixth"; }
	std::string filename = "basissets/";
	auto it = names.find(q);
	std::string bn;
	if (it != names.end()) bn = it->second; 
	else {
		it = names.find(0);
		if (it != names.end()) bn = it->second;
		else bn = "sto-3g";
	}
	std::transform(bn.begin(), bn.end(), bn.begin(), ::toupper);
	filename += bn; filename += row; filename += ".basis";
	// Open file, read only
	input.open(filename, std::ifstream::in);
}

void BasisReader::closeFile()
{
	if(input.is_open()){
		input.close();
	}
}

// Read in the number of contracted gaussian basis functions
// associated with an atom of atomic number q in the basis set
int BasisReader::readNbfs(int q)
{
	int nbfs = 0;
	iVector shells = readShells(q);
	for (int i = 0; i < shells.size(); i++){
		nbfs += shells(i);
	}
	return nbfs;
}

// Read in the ith basis function for atom q from the basis file
BF BasisReader::readBF(int q, int i)
{
	// Open the file
	openFile(q);

	std::string delim = ",";
	int l1 = 0, l2 = 0, l3 = 0; // Angular momenta
	Vector c; Vector e; iVector ids; // Coeffs, exps, and prim ids

	// Check it's open
	if (input.is_open()){
		// Parse
		std::string aname = getAtomName(q);
		aname += " "; aname += delim;
		std::size_t position;
		int bfcount = -1;
		std::string line, shell;
		std::getline(input, line);
    
		// Main loop
		while(!input.eof() && bfcount != i){
			position = line.find(aname);
			if (position != std::string::npos){ // Atom type found
				// Get the shell type
				position = line.find(delim);
				shell = line.substr(0, position);
	
				// Calculate total no. of exponents
				long nexps = std::count(line.begin(), line.end(), ',') - 1;
				// Copy in
				Vector tempexps((int) nexps);
				// Get rid of shell declaration
				line.erase(0, position+delim.length());
				position = line.find(delim);
				std::string temp2;
				int counter = 0;
	
				// Now front bits have been removed, copy in exps 1 by 1
				while (position != std::string::npos) {
					line.erase(0, position+delim.length());
					position = line.find(delim);
					temp2 = line.substr(0, position); // Get the exponent
					// Get rid of extraneous spaces
					temp2.erase(std::remove(temp2.begin(), temp2.end(), ' '), temp2.end()); 
					tempexps[counter] = std::stod(temp2);
					counter++;
				}
	
				// How many functions does this shell have?
				// N.b. cartesian not spherical gaussians
				int lmult = 1;
				if(shell == "p") { lmult = 3; }
				else if (shell == "sp") { lmult = 4; }
				else if (shell == "d") { lmult = 6; }
				else if (shell == "f") { lmult = 10; }
				else if (shell == "g") { lmult = 15; }
				else if (shell == "h") { lmult = 21; }
				else if (shell == "i") { lmult = 28; }
				else if (shell == "k") { lmult = 36; }
				else if (shell == "l") { lmult = 45; }
 	
				// Iterate through bfs to find right one
				int sublmult;
				std::getline(input, line);
				while (bfcount!=i && line.at(0) == 'c'){
					sublmult = 1;
					while (bfcount != i && sublmult < lmult+1){
						bfcount++;
						sublmult++;
					}
					if (bfcount!=i){
						std::getline(input, line);
					}
				}
	
				if (bfcount == i) { // Found it 
					// Find out which exponents to use
					position = line.find(delim);
	  
					// Get rid of exponent number declaration
					line.erase(0, position+delim.length());
					position = line.find(delim);
					temp2 = line.substr(0, position);
					std::size_t p = temp2.find('.');
					int start = std::stoi(temp2.substr(0, p));
					int end = std::stoi(temp2.substr(p+1, temp2.length()));

					// Now resize e and ids, and copy in
					e.resize(end - start + 1);
					ids.resize(end - start + 1);
					for (int j = 0; j < end - start + 1; j++){
						e[j] = tempexps(start+j-1); // Take account of zero-indexing
						ids[j] = start + j - 1;
					}

					// Get the coeffs
					c.resize(end - start + 1);
					for (int j = 0; j < end - start + 1; j++){
						line.erase(0, position+delim.length());
						position = line.find(delim);
						temp2 = line.substr(0, position);
						temp2.erase(std::remove(temp2.begin(), temp2.end(), ' '), temp2.end());
						c[j] = std::stod(temp2);
					}

					// Work out the lnums
					// sublmult is 1 too high
					sublmult--;
					switch(lmult){
						case 1: { // s type
							l1 = l2 = l3 = 0;
							break;
						}
						case 3: { // p type
							switch(sublmult){
								case 1: { //px 
									l1 = 1; l2 = l3 = 0;
									break;
								}
								case 2: { //py
									l1 = l3 = 0; l2 = 1;
									for (int index = 0; index < ids.size(); index++) { ids[index] += tempexps.size(); }
									break;
								}
								case 3: { //pz
									l1 = l2 = 0; l3 = 1;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 2*tempexps.size(); }
									break;
								}
							}
							break;
						}
						case 4: { // sp type
							switch(sublmult){
								case 1:{ // s
									l1 = l2 = l3 = 0;
									break;
								}
								case 2:{ // px
									l1 = 1; l2 = l3 = 0;
									for (int index = 0; index< ids.size(); index++) { ids[index] += tempexps.size(); }
									break;
								}
								case 3:{ // py
									l1 = l3 = 0; l2 = 1;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 2*tempexps.size(); }
									break;
								}
								case 4:{ //pz
									l1 = l2 = 0; l3 = 1;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 3*tempexps.size(); }
									break;
								}
							}
							break;
						}
						case 6: { // d type
							switch(sublmult){
								case 1:{ // dzz
									l3 = 2; l1 = l2 = 0;
									break;
								}
								case 2:{ // dyz
									l3 = l2 = 1; l1 = 0;
									for (int index = 0; index< ids.size(); index++) { ids[index] += tempexps.size(); }
									break;
								}
								case 3:{ // dxz
									l1 = l3 = 1; l2 = 0;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 2*tempexps.size(); }
									break;
								}
								case 4:{ // dyy
									l2 = 2; l3 = l1 = 0;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 3*tempexps.size(); }
									break;
								}
								case 5:{ // dxy 
									l3 = 0; l2 = l1 = 1;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 4*tempexps.size(); }
									break;
								}
								case 6:{ // dxx
									l1 = 2; l2 = l3 = 0;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 5*tempexps.size(); }
									break;
								}
							}
							break;
						}
						case 10:{ // f type
							switch(sublmult){
								case 1:{ //fzzz
									l1 = l2 = 0; l3 = 3;
									break;
								}
								case 2:{ //fyzz
									l1 = 0; l2 = 1; l3 = 2;
									for (int index = 0; index< ids.size(); index++) { ids[index] += tempexps.size(); }
									break;
								}
								case 3:{ //fxzz
									l1 = 1; l2 = 0; l3 = 2;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 2*tempexps.size(); }
									break;
								}
								case 4:{ //fyyz
									l1 = 0; l2 = 2; l3 = 1;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 3*tempexps.size(); }
									break;
								}
								case 5:{ //fxyz
									l1 = 1; l2 = 1; l3 = 1;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 4*tempexps.size(); }
									break;
								}
								case 6:{ //fxxz
									l1 = 2; l2 = 0; l3 = 1;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 5*tempexps.size(); }
									break;
								}
								case 7:{ //fyyy
									l1 = l3 = 0; l2 = 3;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 6*tempexps.size(); }
									break;
								}
								case 8:{ // fxyy
									l1 = 1; l2 = 2; l3 = 0;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 7*tempexps.size(); }
									break;
								}
								case 9:{ // fxxy
									l1 = 2; l2 = 1; l3 = 0;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 8*tempexps.size(); }
									break;
								}
								case 10:{ // fxxx
									l1 = 3; l2 = 0; l3 = 0;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 9*tempexps.size(); }
									break;
								}
							}
							break;
						}
						case 15: { // g-type
							switch(sublmult) { 
								case 1: { 
									l3 = 4; l2 = l1 = 0;
									break;
								}
								case 2: { 
									l3 = 3; l2 = 1;  l1 = 0;
									for (int index = 0; index< ids.size(); index++) { ids[index] += tempexps.size(); }
									break;
								}
								case 3: {
									l3 = 3; l2 = 0;  l1 = 1;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 2*tempexps.size(); }
									break;
								}
								case 4: {
									l3 = 2; l2 = 2;  l1 = 0;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 3*tempexps.size(); }
									break;
								}
								case 5: {
									l3 = 2; l2 = 1;  l1 = 1;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 4*tempexps.size(); }
									break;
								}
								case 6: {
									l3 = 2; l2 = 0;  l1 = 2;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 5*tempexps.size(); }
									break;
								}
								case 7: {
									l3 = 1; l2 = 3;  l1 = 0;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 6*tempexps.size(); }
									break;
								}
								case 8: {
									l3 = 1; l2 = 2;  l1 = 1;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 7*tempexps.size(); }
									break;
								}
								case 9: {
									l3 = 1; l2 = 1;  l1 = 2;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 8*tempexps.size(); }
									break;
								}
								case 10: {
									l3 = 1; l2 = 0;  l1 = 3;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 9*tempexps.size(); }
									break;
								}
								case 11: {
									l3 = 0; l2 = 4;  l1 = 0;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 10*tempexps.size(); }
									break;
								}
								case 12: {
									l3 = 0; l2 = 3;  l1 = 1;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 11*tempexps.size(); }
									break;
								}
								case 13: {
									l3 = 0; l2 = 2;  l1 = 2;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 12*tempexps.size(); }
									break;
								}
								case 14: {
									l3 = 0; l2 = 1;  l1 = 3;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 13*tempexps.size(); }
									break;
								}
								case 15: {
									l3 = 0; l2 = 0;  l1 = 4;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 14*tempexps.size(); }
									break;
								}
								default: { l1 = l2 = l3 = 0; }
							}
							break;
						}
						case 21: { // h-type
							switch(sublmult) { 
								case 1: { 
									l3 = 5; l2 = l1 = 0;
									break;
								}
								case 2: { 
									l3 = 4; l2 = 1;  l1 = 0;
									for (int index = 0; index< ids.size(); index++) { ids[index] += tempexps.size(); }
									break;
								}
								case 3: {
									l3 = 4; l2 = 0;  l1 = 1;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 2*tempexps.size(); }
									break;
								}
								case 4: {
									l3 = 3; l2 = 2;  l1 = 0;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 3*tempexps.size(); }
									break;
								}
								case 5: {
									l3 = 3; l2 = 1;  l1 = 1;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 4*tempexps.size(); }
									break;
								}
								case 6: {
									l3 = 3; l2 = 0;  l1 = 2;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 5*tempexps.size(); }
									break;
								}
								case 7: {
									l3 = 2; l2 = 3;  l1 = 0;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 6*tempexps.size(); }
									break;
								}
								case 8: {
									l3 = 2; l2 = 2;  l1 = 1;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 7*tempexps.size(); }
									break;
								}
								case 9: {
									l3 = 2; l2 = 1;  l1 = 2;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 8*tempexps.size(); }
									break;
								}
								case 10: {
									l3 = 2; l2 = 0;  l1 = 3;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 9*tempexps.size(); }
									break;
								}
								case 11: {
									l3 = 1; l2 = 4;  l1 = 0;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 10*tempexps.size(); }
									break;
								}
								case 12: {
									l3 = 1; l2 = 3;  l1 = 1;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 11*tempexps.size(); }
									break;
								}
								case 13: {
									l3 = 1; l2 = 2;  l1 = 2;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 12*tempexps.size(); }
									break;
								}
								case 14: {
									l3 = 1; l2 = 1;  l1 = 3;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 13*tempexps.size(); }
									break;
								}
								case 15: {
									l3 = 1; l2 = 0;  l1 = 4;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 14*tempexps.size(); }
									break;
								}
								case 16: {
									l3 = 0; l2 = 5;  l1 = 0;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 15*tempexps.size(); }
									break;
								}
								case 17: {
									l3 = 0; l2 = 4;  l1 = 1;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 16*tempexps.size(); }
									break;
								}
								case 18: {
									l3 = 0; l2 = 3;  l1 = 2;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 17*tempexps.size(); }
									break;
								}
								case 19: {
									l3 = 0; l2 = 2;  l1 = 3;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 18*tempexps.size(); }
									break;
								}
								case 20: {
									l3 = 0; l2 = 1;  l1 = 4;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 19*tempexps.size(); }
									break;
								}
								case 21: {
									l3 = 0; l2 = 0;  l1 = 5;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 20*tempexps.size(); }
									break;
								}
								default: { l1 = l2 = l3 = 0; }
							}
							break;
						}
						case 28: { // i-type
							switch(sublmult) { 
								case 1: { 
									l3 = 6; l2 = l1 = 0;
									break;
								}
								case 2: { 
									l3 = 5; l2 = 1;  l1 = 0;
									for (int index = 0; index< ids.size(); index++) { ids[index] += tempexps.size(); }
									break;
								}
								case 3: {
									l3 = 5; l2 = 0;  l1 = 1;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 2*tempexps.size(); }
									break;
								}
								case 4: {
									l3 = 4; l2 = 2;  l1 = 0;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 3*tempexps.size(); }
									break;
								}
								case 5: {
									l3 = 4; l2 = 1;  l1 = 1;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 4*tempexps.size(); }
									break;
								}
								case 6: {
									l3 = 4; l2 = 0;  l1 = 2;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 5*tempexps.size(); }
									break;
								}
								case 7: {
									l3 = 3; l2 = 3;  l1 = 0;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 6*tempexps.size(); }
									break;
								}
								case 8: {
									l3 = 3; l2 = 2;  l1 = 1;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 7*tempexps.size(); }
									break;
								}
								case 9: {
									l3 = 3; l2 = 1;  l1 = 2;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 8*tempexps.size(); }
									break;
								}
								case 10: {
									l3 = 3; l2 = 0;  l1 = 3;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 9*tempexps.size(); }
									break;
								}
								case 11: {
									l3 = 2; l2 = 4;  l1 = 0;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 10*tempexps.size(); }
									break;
								}
								case 12: {
									l3 = 2; l2 = 3;  l1 = 1;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 11*tempexps.size(); }
									break;
								}
								case 13: {
									l3 = 2; l2 = 2;  l1 = 2;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 12*tempexps.size(); }
									break;
								}
								case 14: {
									l3 = 2; l2 = 1;  l1 = 3;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 13*tempexps.size(); }
									break;
								}
								case 15: {
									l3 = 2; l2 = 0;  l1 = 4;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 14*tempexps.size(); }
									break;
								}
								case 16: {
									l3 = 1; l2 = 5;  l1 = 0;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 15*tempexps.size(); }
									break;
								}
								case 17: {
									l3 = 1; l2 = 4;  l1 = 1;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 16*tempexps.size(); }
									break;
								}
								case 18: {
									l3 = 1; l2 = 3;  l1 = 2;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 17*tempexps.size(); }
									break;
								}
								case 19: {
									l3 = 1; l2 = 2;  l1 = 3;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 18*tempexps.size(); }
									break;
								}
								case 20: {
									l3 = 1; l2 = 1;  l1 = 4;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 19*tempexps.size(); }
									break;
								}
								case 21: {
									l3 = 1; l2 = 0;  l1 = 5;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 20*tempexps.size(); }
									break;
								}
								case 22: {
									l3 = 0; l2 = 6;  l1 = 0;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 21*tempexps.size(); }
									break;
								}
								case 23: {
									l3 = 0; l2 = 5;  l1 = 1;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 22*tempexps.size(); }
									break;
								}
								case 24: {
									l3 = 0; l2 = 4;  l1 = 2;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 23*tempexps.size(); }
									break;
								}
								case 25: {
									l3 = 0; l2 = 3;  l1 = 3;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 24*tempexps.size(); }
									break;
								}
								case 26: {
									l3 = 0; l2 = 2;  l1 = 4;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 25*tempexps.size(); }
									break;
								}
								case 27: {
									l3 = 0; l2 = 1;  l1 = 5;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 26*tempexps.size(); }
									break;
								}
								case 28: {
									l3 = 0; l2 = 0;  l1 = 6;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 27*tempexps.size(); }
									break;
								}
								default: { l1 = l2 = l3 = 0; }
							}
							break;
						}
						case 36: { // k-type
							switch(sublmult) { 
								case 1: { 
									l3 = 7; l2 = l1 = 0;
									break;
								}
								case 2: { 
									l3 = 6; l2 = 1;  l1 = 0;
									for (int index = 0; index< ids.size(); index++) { ids[index] += tempexps.size(); }
									break;
								}
								case 3: {
									l3 = 6; l2 = 0;  l1 = 1;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 2*tempexps.size(); }
									break;
								}
								case 4: {
									l3 = 5; l2 = 2;  l1 = 0;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 3*tempexps.size(); }
									break;
								}
								case 5: {
									l3 = 5; l2 = 1;  l1 = 1;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 4*tempexps.size(); }
									break;
								}
								case 6: {
									l3 = 5; l2 = 0;  l1 = 2;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 5*tempexps.size(); }
									break;
								}
								case 7: {
									l3 = 4; l2 = 3;  l1 = 0;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 6*tempexps.size(); }
									break;
								}
								case 8: {
									l3 = 4; l2 = 2;  l1 = 1;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 7*tempexps.size(); }
									break;
								}
								case 9: {
									l3 = 4; l2 = 1;  l1 = 2;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 8*tempexps.size(); }
									break;
								}
								case 10: {
									l3 = 4; l2 = 0;  l1 = 3;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 9*tempexps.size(); }
									break;
								}
								case 11: {
									l3 = 3; l2 = 4;  l1 = 0;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 10*tempexps.size(); }
									break;
								}
								case 12: {
									l3 = 3; l2 = 3;  l1 = 1;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 11*tempexps.size(); }
									break;
								}
								case 13: {
									l3 = 3; l2 = 2;  l1 = 2;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 12*tempexps.size(); }
									break;
								}
								case 14: {
									l3 = 3; l2 = 1;  l1 = 3;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 13*tempexps.size(); }
									break;
								}
								case 15: {
									l3 = 3; l2 = 0;  l1 = 4;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 14*tempexps.size(); }
									break;
								}
								case 16: {
									l3 = 2; l2 = 5;  l1 = 0;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 15*tempexps.size(); }
									break;
								}
								case 17: {
									l3 = 2; l2 = 4;  l1 = 1;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 16*tempexps.size(); }
									break;
								}
								case 18: {
									l3 = 2; l2 = 3;  l1 = 2;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 17*tempexps.size(); }
									break;
								}
								case 19: {
									l3 = 2; l2 = 2;  l1 = 3;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 18*tempexps.size(); }
									break;
								}
								case 20: {
									l3 = 2; l2 = 1;  l1 = 4;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 19*tempexps.size(); }
									break;
								}
								case 21: {
									l3 = 2; l2 = 0;  l1 = 5;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 20*tempexps.size(); }
									break;
								}
								case 22: {
									l3 = 1; l2 = 6;  l1 = 0;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 21*tempexps.size(); }
									break;
								}
								case 23: {
									l3 = 1; l2 = 5;  l1 = 1;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 22*tempexps.size(); }
									break;
								}
								case 24: {
									l3 = 1; l2 = 4;  l1 = 2;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 23*tempexps.size(); }
									break;
								}
								case 25: {
									l3 = 1; l2 = 3;  l1 = 3;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 24*tempexps.size(); }
									break;
								}
								case 26: {
									l3 = 1; l2 = 2;  l1 = 4;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 25*tempexps.size(); }
									break;
								}
								case 27: {
									l3 = 1; l2 = 1;  l1 = 5;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 26*tempexps.size(); }
									break;
								}
								case 28: {
									l3 = 1; l2 = 0;  l1 = 6;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 27*tempexps.size(); }
									break;
								}
								case 29: {
									l3 = 0; l2 = 7;  l1 = 0;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 28*tempexps.size(); }
									break;
								}
								case 30: {
									l3 = 0; l2 = 6;  l1 = 1;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 29*tempexps.size(); }
									break;
								}
								case 31: {
									l3 = 0; l2 = 5;  l1 = 2;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 30*tempexps.size(); }
									break;
								}
								case 32: {
									l3 = 0; l2 = 4;  l1 = 3;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 31*tempexps.size(); }
									break;
								}
								case 33: {
									l3 = 0; l2 = 3;  l1 = 4;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 32*tempexps.size(); }
									break;
								}
								case 34: {
									l3 = 0; l2 = 2;  l1 = 5;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 33*tempexps.size(); }
									break;
								}
								case 35: {
									l3 = 0; l2 = 1;  l1 = 6;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 34*tempexps.size(); }
									break;
								}
								case 36: {
									l3 = 0; l2 = 0;  l1 = 7;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 35*tempexps.size(); }
									break;
								}
								default: { l1 = l2 = l3 = 0; }
							}
							break;
						}
						case 45: { // l-type
							switch(sublmult) { 
								case 1: { 
									l3 = 8; l2 = l1 = 0;
									break;
								}
								case 2: { 
									l3 = 7; l2 = 1;  l1 = 0;
									for (int index = 0; index< ids.size(); index++) { ids[index] += tempexps.size(); }
									break;
								}
								case 3: {
									l3 = 7; l2 = 0;  l1 = 1;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 2*tempexps.size(); }
									break;
								}
								case 4: {
									l3 = 6; l2 = 2;  l1 = 0;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 3*tempexps.size(); }
									break;
								}
								case 5: {
									l3 = 6; l2 = 1;  l1 = 1;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 4*tempexps.size(); }
									break;
								}
								case 6: {
									l3 = 6; l2 = 0;  l1 = 2;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 5*tempexps.size(); }
									break;
								}
								case 7: {
									l3 = 5; l2 = 3;  l1 = 0;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 6*tempexps.size(); }
									break;
								}
								case 8: {
									l3 = 5; l2 = 2;  l1 = 1;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 7*tempexps.size(); }
									break;
								}
								case 9: {
									l3 = 5; l2 = 1;  l1 = 2;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 8*tempexps.size(); }
									break;
								}
								case 10: {
									l3 = 5; l2 = 0;  l1 = 3;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 9*tempexps.size(); }
									break;
								}
								case 11: {
									l3 = 4; l2 = 4;  l1 = 0;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 10*tempexps.size(); }
									break;
								}
								case 12: {
									l3 = 4; l2 = 3;  l1 = 1;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 11*tempexps.size(); }
									break;
								}
								case 13: {
									l3 = 4; l2 = 2;  l1 = 2;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 12*tempexps.size(); }
									break;
								}
								case 14: {
									l3 = 4; l2 = 1;  l1 = 3;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 13*tempexps.size(); }
									break;
								}
								case 15: {
									l3 = 4; l2 = 0;  l1 = 4;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 14*tempexps.size(); }
									break;
								}
								case 16: {
									l3 = 3; l2 = 5;  l1 = 0;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 15*tempexps.size(); }
									break;
								}
								case 17: {
									l3 = 3; l2 = 4;  l1 = 1;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 16*tempexps.size(); }
									break;
								}
								case 18: {
									l3 = 3; l2 = 3;  l1 = 2;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 17*tempexps.size(); }
									break;
								}
								case 19: {
									l3 = 3; l2 = 2;  l1 = 3;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 18*tempexps.size(); }
									break;
								}
								case 20: {
									l3 = 3; l2 = 1;  l1 = 4;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 19*tempexps.size(); }
									break;
								}
								case 21: {
									l3 = 3; l2 = 0;  l1 = 5;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 20*tempexps.size(); }
									break;
								}
								case 22: {
									l3 = 2; l2 = 6;  l1 = 0;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 21*tempexps.size(); }
									break;
								}
								case 23: {
									l3 = 2; l2 = 5;  l1 = 1;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 22*tempexps.size(); }
									break;
								}
								case 24: {
									l3 = 2; l2 = 4;  l1 = 2;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 23*tempexps.size(); }
									break;
								}
								case 25: {
									l3 = 2; l2 = 3;  l1 = 3;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 24*tempexps.size(); }
									break;
								}
								case 26: {
									l3 = 2; l2 = 2;  l1 = 4;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 25*tempexps.size(); }
									break;
								}
								case 27: {
									l3 = 2; l2 = 1;  l1 = 5;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 26*tempexps.size(); }
									break;
								}
								case 28: {
									l3 = 2; l2 = 0;  l1 = 6;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 27*tempexps.size(); }
									break;
								}
								case 29: {
									l3 = 1; l2 = 7;  l1 = 0;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 28*tempexps.size(); }
									break;
								}
								case 30: {
									l3 = 1; l2 = 6;  l1 = 1;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 29*tempexps.size(); }
									break;
								}
								case 31: {
									l3 = 1; l2 = 5;  l1 = 2;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 30*tempexps.size(); }
									break;
								}
								case 32: {
									l3 = 1; l2 = 4;  l1 = 3;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 31*tempexps.size(); }
									break;
								}
								case 33: {
									l3 = 1; l2 = 3;  l1 = 4;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 32*tempexps.size(); }
									break;
								}
								case 34: {
									l3 = 1; l2 = 2;  l1 = 5;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 33*tempexps.size(); }
									break;
								}
								case 35: {
									l3 = 1; l2 = 1;  l1 = 6;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 34*tempexps.size(); }
									break;
								}
								case 36: {
									l3 = 1; l2 = 0;  l1 = 7;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 35*tempexps.size(); }
									break;
								}
								case 37: {
									l3 = 0; l2 = 8;  l1 = 0;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 36*tempexps.size(); }
									break;
								}
								case 38: {
									l3 = 0; l2 = 7;  l1 = 1;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 37*tempexps.size(); }
									break;
								}
								case 39: {
									l3 = 0; l2 = 6;  l1 = 2;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 38*tempexps.size(); }
									break;
								}
								case 40: {
									l3 = 0; l2 = 5;  l1 = 3;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 39*tempexps.size(); }
									break;
								}
								case 41: {
									l3 = 0; l2 = 4;  l1 = 4;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 40*tempexps.size(); }
									break;
								}
								case 42: {
									l3 = 0; l2 = 3;  l1 = 5;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 41*tempexps.size(); }
									break;
								}
								case 43: {
									l3 = 0; l2 = 2;  l1 = 6;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 42*tempexps.size(); }
									break;
								}
								case 44: {
									l3 = 0; l2 = 1;  l1 = 7;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 43*tempexps.size(); }
									break;
								}
								case 45: {
									l3 = 0; l2 = 0;  l1 = 8;
									for (int index = 0; index< ids.size(); index++) { ids[index] += 44*tempexps.size(); }
									break;
								}
								default: { l1 = l2 = l3 = 0; }
							}
							break;
						}
						default:{
							l1 = l2 = l3 = 0;
						}
					}
					position = std::string::npos; // Exit loop
				}

			} else { // Move on to the next line
				std::getline(input, line);
			}
		}
	} else {
		throw(Error("READBF", "Unable to open basis file."));
	}
	closeFile();
	BF b(c, l1, l2, l3, e, ids);
	return b;
}

// Return a vector of the number of basis functions in each shell
// where the length of the vector is the number of shells
iVector BasisReader::readShells(int q)
{
	iVector shells(10); // Can't cope with higher than l functions, so 10 is max!	  
	// Open the file
	openFile(q);
  
	std::string delim = " ,";
	int counter = 0;
	// Check it's open
	if (input.is_open()){
		// Parse
		std::string aname = getAtomName(q);
		aname += delim;
		std::string line, temp;
		int lmult;
		std::getline(input, line);

		while(!input.eof()){
			if (line.find(aname) != std::string::npos){
				// Count the number of c lines that follow
				temp = line.at(0);
				int nbfs = 0;

				//Work out what shell type it is -
				// note that this calculates no. of cartesian gaussians
				if (temp == "s") { lmult = 1; }
				else if (temp == "sp") { lmult = 4; }
				else if (temp == "p") { lmult = 3; }
				else if (temp == "d") { lmult = 6; }
				else if (temp == "f") { lmult = 10; }
				else if (temp == "g") { lmult = 15; }
				else if (temp == "h") { lmult = 21; }
				else if (temp == "i") { lmult = 28; }
				else if (temp == "k") { lmult = 36; }
				else if (temp == "l") { lmult = 45; }

				std::getline(input, line);
				bool section = true; 
				while (section){
					nbfs += lmult;
					std::getline(input, line);
					section = line.length() > 0;
					if (section) section = line.at(0) == 'c'; 
				}
				shells[counter] = nbfs;
				counter++;
			} else { // Get next line
				std::getline(input, line);
			}
		}
	} else {
		throw(Error("IOERR", "Could not open basis file."));
	}

	// resize the return vector
	shells.conservativeResize(counter);
  
	// Close the file
	closeFile();
  
	return shells;
}

// Return a vector with the l-angular momentum quantum numbers of each
// shell belonging to atom q.
iVector BasisReader::readLnums(int q)
{
	iVector lnums(10); // See readShells
	// Open the file
	openFile(q);
  
	std::string delim = " ,";
	int counter = 0;
	// Check it's open
	if (input.is_open()){
		// Parse
		std::string aname = getAtomName(q);
		aname += delim;
		std::string line, temp;
		int lmult;
		std::getline(input, line);

		while(!input.eof()){
			if (line.find(aname) != std::string::npos){
				// Count the number of c lines that follow
				temp = line.at(0);

				//Work out what shell type it is 
				if (temp == "s") { lmult = 0; }
				else if (temp == "sp") { lmult = 1; }
				else if (temp == "p") { lmult = 1; }
				else if (temp == "d") { lmult = 2; }
				else if (temp == "f") { lmult = 3; }
				else if (temp == "g") { lmult = 4; }
				else if (temp == "h") { lmult = 5; }
				else if (temp == "i") { lmult = 6; }
				else if (temp == "k") { lmult = 7; }
				else if (temp == "l") { lmult = 8; }

				lnums[counter] = lmult;
				counter++;
				std::getline(input, line);
			} else { // Get next line
				std::getline(input, line);
			}
		}
	} else {
		throw(Error("IOERR", "Could not open basis file."));
	}

	// resize the return vector
	lnums.conservativeResize(counter);
  
	// Close the file
	closeFile();
  
	return lnums;
}

// Now do these for a complete set of qs
iVector BasisReader::readShells(iVector& qs)
{
	iVector shells;
	iVector temp;
	shells = readShells(qs(0));
	for (int i = 1; i < qs.size(); i++){
		int s = shells.size();
		temp = readShells(qs(i));
		int t = temp.size();
		shells.conservativeResize(s+t);
		for (int j = 0; j < t; j++) {
			shells[j+s] = temp[j];
		}
	}
	return shells;
}

iVector BasisReader::readLnums(iVector& qs)
{
	iVector lnums;
	iVector temp;
	lnums = readLnums(qs(0));
	for (int i = 1; i < qs.size(); i++){
		int l = lnums.size();
		temp = readLnums(qs(i));
		int t = temp.size();
		lnums.conservativeResize(l+t);
		for (int j = 0; j < t; j++){
			lnums[j+l] = temp[j];
		}
	}
	return lnums;
}

void BasisReader::readShellBasis(Basis& b, int q, double *pos, int atom, bool df) {
	
	openFile(q);
	std::string delim = ",";
	
	if (input.is_open()) {
		
		// Parse
		std::string aname = getAtomName(q);
		aname += " "; aname += delim;
		std::string line, temp;
		int lmult;
		std::getline(input, line);

		while(!input.eof()){
			
			std::size_t position = line.find(aname);
			if (position != std::string::npos){
				temp = line.at(0);

				//Work out what shell type it is 
				if (temp == "s") { lmult = 0; }
				else if (temp == "sp") { lmult = 1; }
				else if (temp == "p") { lmult = 1; }
				else if (temp == "d") { lmult = 2; }
				else if (temp == "f") { lmult = 3; }
				else if (temp == "g") { lmult = 4; }
				else if (temp == "h") { lmult = 5; }
				else if (temp == "i") { lmult = 6; }
				else if (temp == "k") { lmult = 7; }
				else if (temp == "l") { lmult = 8; }
			  
				std::vector<libint2::real_t> exps;
				line.erase(0, position+delim.length());
				position = line.find(delim);
				std::string temp2;
	
				// Now front bits have been removed, copy in exps 1 by 1
				while (position != std::string::npos) {
					line.erase(0, position+delim.length());
					position = line.find(delim);
					temp2 = line.substr(0, position); // Get the exponent
					// Get rid of extraneous spaces
					temp2.erase(std::remove(temp2.begin(), temp2.end(), ' '), temp2.end()); 
					exps.push_back(std::stod(temp2));
				}
			  
				std::getline(input, line);
				while (line.at(0) == 'c'){
					position = line.find(delim);
					line.erase(0, position+delim.length());
					position = line.find(delim);
					temp2 = line.substr(0, position);
					std::size_t p = temp2.find('.');
					int start = std::stoi(temp2.substr(0, p));
					int end = std::stoi(temp2.substr(p+1, temp2.length()));
					
					std::vector<libint2::real_t> e; 
					for (int i = start-1; i < end; i++) e.push_back(exps[i]);
					
					std::vector<std::vector <libint2::real_t>> coeffs;
					std::vector<libint2::real_t> cvals;
					while (position != std::string::npos) {
						line.erase(0, position+delim.length());
						position = line.find(delim);
						temp2 = line.substr(0, position);
						temp2.erase(std::remove(temp2.begin(), temp2.end(), ' '), temp2.end());
						cvals.push_back(std::stod(temp2));
					}
					coeffs.push_back(cvals);
					b.addShell(lmult, e, coeffs, pos, atom, df);
					std::getline(input, line);
				}
				
				
			} else { // Get next line
				std::getline(input, line);
			}
		}
		
		closeFile();
		
	} else {
		throw(Error("IOERR", "Could not open basis file."));
	}
}

ECP BasisReader::readECP(int q, ECPBasis& ecpset, double* center) {
	openFile(-q);
	
	std::string delim = ",";

	ECP newECP(center); 
	// Check it's open
	if (input.is_open()){
		// Parse
		std::string aname = getAtomName(q);
		aname += delim;
		std::size_t position;
		std::string line, token, word;
		std::getline(input, line);
		int ncore, maxl;
    
		// Main loop
		while(!input.eof()){
			position = line.find(aname);
			if (position != std::string::npos){ // Atom type found
				// Get the shell type
				line.erase(0, position+aname.length());

				position = line.find(delim);
				if (position != std::string::npos) {
					token = line.substr(0, position);
					ncore = std::stoi(token);
						
					auto it = ecpset.core_electrons.find(q);
					if (it == ecpset.core_electrons.end()) ecpset.core_electrons[q] = ncore; 
					
					line.erase(0, position+delim.length());
					position = line.find(delim);
					
					if (position != std::string::npos) {
						token = line.substr(0, position);
						maxl = std::stoi(token);
					}
				}
				
				for (int i = 0; i<= maxl; i++) {
					std::getline(input, line);
					position = line.find(';');
					if (position != std::string::npos) {
						token = line.substr(0, position); 
						int nprims = std::stoi(token);
						line.erase(0, position+1); 
						
						int j = 0;
						position = line.find(';'); 
						while (j < nprims && position != std::string::npos) {
							token = line.substr(0, position);
							line.erase(0, position+1);
							position = token.find(',');
							int n = 0; double d = 0.0, a = 1.0; 
							if (position != std::string::npos) {
								word = token.substr(0, position); 
								n = std::stoi(word);
								token.erase(0, position+1);
								position = token.find(',');
								word = token.substr(0, position);
								a = std::stod(word);
								token.erase(0, position+1);
								d = std::stod(token);
							}
							position = line.find(';');
							int l = i == 0 ? maxl : i-1; 
							newECP.addPrimitive(n, l, a, d, true);
							j++; 
						}
						
					}
				}
			}
			std::getline(input, line);
		}
	} else {
		throw(Error("READECP", "Unable to open ECP basis file."));
	}
	closeFile();			
	return newECP; 
}
