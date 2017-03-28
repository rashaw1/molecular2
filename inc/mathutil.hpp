/*
 *
 *   PURPOSE: To declare a number of mathematical utility functions needed by 
 *            various other libraries in the molecular suite.
 * 
 *   LIST OF FUNCTIONS:
 *           fact(int i) : calculates i factorial
 *           fact2(int i) : calculates i double factorial
 *           boys(x, mmax, mmin) : calculates the boys function F_m(x)
 *                                 for m in the range mmax to mmin
 *                                 mmin = 0 by default. Returns values as
 *                                 a vector.
 *           binom(int m, int n) : calculates the binomial coefficient (n m)
 *           rmultiply(mat, vec), lmultiply(vec, mat) - right and left mult.
 *                               of a matrix by a vector.
 *
 *   DATE          AUTHOR            CHANGES
 *   =======================================================================
 *   27/08/15      Robert Shaw       Original code.
 *   28/08/15      Robert Shaw       Added boys function calculator.
 *   02/09/15      Robert Shaw       Added binomial coeff. calculator.
 *   05/09/15      Robert Shaw       Added clebsch and sphernorm.
 */

#ifndef MATHUTILHEADERDEF
#define MATHUTILHEADERDEF

#include "eigen_wrapper.hpp"
#include "tensor4.hpp"
#include <vector>

// Functions to calculate the factorial and double factorial of an integer i
unsigned long int fact(int i);
unsigned long int fact2(int i);
void factArray(int i, double *values);
void fact2Array(int i, double *values);

// Calculates the binomial coefficient (n m)T
unsigned long int binom(int n, int m); 

// Calculate the Clebsch-Gordon coefficient C^{l, m}_{t, u, v} 
// needed for transformation to spherical basis
double clebsch(int l, int m, int t, int u, double v);

// Calculate the normalisation constant for a spherical gaussian
// with quantum numbers l,m
double sphernorm(int l, int m);

// Form a section of the cartesian-to-spherical basis transformation
// matrix, from a pre-calculated list of coefficients. 
// See Schlegel and Frisch, Int. J. Q. Chem., 54, 83-87 (1995)
void formTransMat(Matrix& mat, int row, int col, int l, int m);

Vector rmultiply(const Matrix& mat, const Vector& v);
Vector lmultiply(const Vector& v, const Matrix& mat);
Vector cross(const Vector& v1, const Vector& v2);
Matrix d_cross(const Vector& v1, const Vector& v2);
Matrix d_cross_ab(const Vector& v1, const Vector& v2, const Matrix& da, const Matrix& db); 
double ncross(const Vector& v1, const Vector& v2);
Vector d_ncross(const Vector& v1, const Vector& v2);

Matrix pseudo_inverse(Matrix& mat, double threshold = 1e-8); 

Vector get_quat(const Matrix& x, const Matrix& y); 
std::vector<Matrix> get_q_der(const Matrix& x, const Matrix& y); 
Tensor4 get_R_der(const Matrix& x, const Matrix& y); 
Tensor4 get_F_der(const Matrix& x, const Matrix& y); 
Matrix build_F(const Matrix& x, const Matrix& y); 
bool is_linear(const Matrix& xyz, const Matrix& x0); 
Vector get_exp_map(const Matrix& xyz, const Matrix& x0); 
std::vector<Matrix> get_exp_map_der(const Matrix& xyz, const Matrix& x0); 

const double GAMMA[30] = {
	1.7724538509055,
	1.0,
	0.88622692545275,
	1.0,
	1.3293403881791,
	2.0,
	3.3233509704478,
	6.0,
	11.631728396567,
	24.0,
	52.342777784553,
	120.0,
	287.88527781504,
	720.0,
	1871.2543057978,
	5040.0,
	14034.407293483,
	40320.0,
	1.1929246199461e5,
	3.62880e5,
	1.1332783889488e6,
	3.628800e6,
	1.1899423083962e7,
	3.9916800e7,
	1.3684336546556e8, 
	4.79001600e8, 
	1.7105420683196e9, 
	6.227020800e9, 
	2.3092317922314e10,
	8.7178291200e10
};


#endif
