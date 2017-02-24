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

// Boys function calculator for the evaluation of molecular integrals over
// gaussian type basis functions. Returns a vector of F_m(x) for m in the range
// mmax to mmin (it uses downwards recursion).
Vector boys(double x, int mmax, int mmin = 0, double PRECISION=1e-14);

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

Matrix pseudo_inverse(Matrix& mat, double threshold = 1e-8); 

Vector get_quat(const Matrix& x, const Matrix& y); 
std::vector<Matrix> get_q_der(const Matrix& x, const Matrix& y); 
Tensor4 get_R_der(const Matrix& x, const Matrix& y); 
Tensor4 get_F_der(const Matrix& x, const Matrix& y); 
Matrix build_F(const Matrix& x, const Matrix& y); 
bool is_linear(const Matrix& xyz, const Matrix& x0); 
Vector get_exp_map(const Matrix& xyz, const Matrix& x0); 
std::vector<Matrix> get_exp_map_der(const Matrix& xyz, const Matrix& x0); 




#endif
