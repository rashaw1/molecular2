#ifndef INTERNAL_COORD_HEAD
#define INTERNAL_COORD_HEAD

#include <vector>

enum Axis { 
	X = 0,
	Y = 1,
	Z = 2
}

struct CartesianXYZ {
	
	int atom;
	double weight;
	bool isAngular, isPeriodic;
	
	CartesianXYZ(int _a, double _w = 1.0) : atom(_a), weight(_w), isAngular(false), isPeriodic(false) { }
	
	bool operator==(const CartesianXYZ& other) { return atom == other.atom; }
	bool operator!=(const CartesianXYZ& other) { return !(*this == other); }
	
	double value(const Matrix& xyz, int c) { return xyz[atom][c] * weight; }
	Matrix derivative(const Matrix& xyz, int c) {
		Matrix derivatives = Matrx::Zero(xyz.rows(), xyz.cols()); 
		derivatives(atom, c) = weight; 
		return derivatives; 
	}
	
};

struct TranslationXYZ {
	
	std::vector<int> atoms;
	std::vector<double> weights; 
	bool isAngular, isPeriodic; 
	
	TranslationXYZ(std::vector<int>& _a, std::vector<double>& _w) : atoms(_a), weights(_w), isAngular(false), isPeriodic(false) { }
	
	bool operator==(const TranslationXYZ& other) { return atoms == other.atoms; }
	bool operator!=(const TranslationXYZ& other) { return !(*this == other); }
	
	double value(const Matrix& xyz, int c) { 
		double sum = 0.0;
		for (int i = 0; i < atoms.size(); i++) sum += xyz(atoms[i], c) * weights[i]; 
		return sum;
	}
	
	Matrix derivative(const Matrix& xyz, int c) {
		Matrix derivatives = Matrx::Zero(xyz.rows(), xyz.cols());
		for (int i = 0; i < atoms.size(); i++) derivatives(atoms[i], c) = weights[i]; 
		return derivatives;
	}
	
}; 

struct Rotator {
	
	std::vector<int> atoms; 
	Matrix x0, valXYZ, derXYZ;
	Vector stored_value, deriv, axis; 
	double norm, dot2; 
	bool linear; 
	
	Rotator(std::vector<int>& _a, const Matrix& xyz) : atoms(_a) {
		init(xyz);
	}
	
	void init(const Matrix& xyz) {
		x0 = xyz; 
		valXYZ = Matrix::Zero(xyz.rows(), xyz.cols()); 
		derXYZ = valXYZ;
		norm = dot2 = 0.0;
		linear = false; 
	}
	
	bool operator==(const Rotator& other) { return atoms == other.atoms; }
	bool operator!=(const Rotator& other) { return !(*this == other); }
	
	void calcAxis(); 
	Vector value(const Matrix& xyz); 
	Matrix derivative(const Matrix& xyz); 
	
};


class InternalCoordinates {
	
};

#endif