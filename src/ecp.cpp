/* Implements ecp.hpp 

   Robert A. Shaw 2016
*/

#include "ecp.hpp"

#include <cmath>
#include <iostream>
#include <algorithm>

// GaussianECP constructor and copy constructor
GaussianECP::GaussianECP() : n(0), l(0), a(0), d(0) {}
GaussianECP::GaussianECP(int _n, int _l, double _a, double _d) : n(_n-2), l(_l), a(_a), d(_d) {}
GaussianECP::GaussianECP(const GaussianECP& other) : n(other.n), l(other.l), a(other.a), d(other.d) {}


// class ECP

ECP::ECP() : N(0), L(-1) {
	center_.resize(3);	
}
ECP::ECP(const double *_center) : N(0), L(-1) {
	center_.resize(3);
	center_[0] = _center[0];
	center_[1] = _center[1];
	center_[2] = _center[2];
}

ECP::ECP(const ECP &other) {
	gaussians = other.gaussians;
	N = other.N;
	L = other.L;
	center_ = other.center_;
}

void ECP::addPrimitive(int n, int l, double a, double d, bool needSort) {
	GaussianECP newEcp(n, l, a, d);
	gaussians.push_back(newEcp);
	N++;
	L = l > L ? l : L;
	if (needSort) sort();
}

void ECP::sort() {
	std::sort(gaussians.begin(), gaussians.end(),
	 [&] (const GaussianECP& g1, const GaussianECP& g2) {return (g1.l < g2.l);});
}

// Evaluate U_l(r), assuming that gaussians sorted by angular momentum
double ECP::evaluate(double r, int l) {
	double value = 0.0;
	int am = 0;
	double r2 = r*r;
	for (int i = 0; i < N; i++) {
		if (gaussians[i].l == l) 
		  value += pow(r, gaussians[i].n) * gaussians[i].d * exp(-gaussians[i].a * r2);
	} 
	return value; 
}

ECPBasis::ECPBasis() : N(0), maxL(-1) {}

void ECPBasis::addECP(ECP &U) {
	basis.push_back(U);
	N++;
	maxL = U.getL() > maxL ? U.getL() : maxL;
}

ECP& ECPBasis::getECP(int i) { return basis[i]; }

int ECPBasis::getECPCore(int q) {
	int core = 0;
	auto it = core_electrons.find(q);
	if (it != core_electrons.end()) core = it->second;
	return core;
}

