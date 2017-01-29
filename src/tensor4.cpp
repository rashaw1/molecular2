/* 
 *
 *   PURPOSE: To implement class Tensor4
 *
 */

#include "tensor4.hpp"
#include <iostream>
#include <cmath>
#include <algorithm>

Tensor4::Tensor4(int a, int b, int c, int d) : w(a), x(b), y(c), z(d)
{
  // Resize the data matrix
  data.resize(a*b*c*d);
}

Tensor4::Tensor4(int a, int b, int c, int d, double val) : w(a), x(b), y(c), z(d)
{
  // Resize and assign the data matrix
  data.resize(a*b*c*d);
  std::fill(data.begin(), data.end(), val);
}

Tensor4::Tensor4(const Tensor4& other)
{
  w = other.w; x = other.x; y = other.y; z = other.z;
  data.resize(w*x*y*z);
  for (int i = 0; i < w; i++){
    for (int j = 0; j < x; j++){
      for (int k = 0; k < y; k++){
        for (int l = 0; l < z; l++){
          data[i*x*y*z + j*y*z + k*z +l] = other(i, j, k, l);
        }
      } 
    }
  }
}  

void Tensor4::resize(int a, int b, int c, int d)
{
  w = a; x = b; y = c; z = d;
  data.resize(a*b*c*d);
}

void Tensor4::assign(int a, int b, int c, int d, double val)
{
  w = a; x = b; y = c; z = d;
  data.resize(a*b*c*d);
  std::fill(data.begin(), data.end(), val);
}

void Tensor4::print() const 
{
  for (int i = 0; i < w; i++){
    for (int j = 0; j < x; j++){
      for (int k = 0; k < y; k++){
	for (int l = 0; l < z; l++){
	  std::cout << i << " " << j << " " << k << " " << l << "   " << data[i*x*y*z + j*y*z + k*z + l] << "\n";
	}
      }
    }
  }
  std::cout << "\n\n";
}

double& Tensor4::operator()(int i, int j, int k, int l)
{
  return data[i*x*y*z + j*y*z + k*z + l];
}

double Tensor4::operator()(int i, int j, int k, int l) const
{
  return data[i*x*y*z + j*y*z + k*z + l];
}

Tensor4& Tensor4::operator=(const Tensor4& other)
{
  resize(other.w, other.x, other.y, other.z);
  for (int i = 0; i < w; i++){
    for (int j = 0; j < x; j++){
      for (int k = 0; k < y; k++){
	for (int l = 0; l < z; l++){
	  data[i*x*y*z + j*y*z + k*z + l] = other(i, j, k, l);
	}
      }
    }
  }
  return *this;
}

Tensor4 Tensor4::operator+(const Tensor4& other) const
{
	Tensor4 retVal(w, x, y, z);
	for (int i = 0; i < w; i++){
		for (int j = 0; j < x; j++){
			for (int k = 0; k < y; k++){
				for (int l = 0; l < z; l++){
					retVal(i, j, k, l) = data[i*x*y*z + j*y*z + k*z + l] + other(i, j, k, l);
				}
			}
		}
	}
	return retVal;
}

Tensor4 Tensor4::operator-(const Tensor4& other) const
{
	Tensor4 retVal(w, x, y, z);
	for (int i = 0; i < w; i++){
		for (int j = 0; j < x; j++){
			for (int k = 0; k < y; k++){
				for (int l = 0; l < z; l++){
					retVal(i, j, k, l) = data[i*x*y*z + j*y*z + k*z +l] - other(i, j, k, l);
				}
			}
		}
	}
	return retVal;
}

// Frobenius norm
double fnorm(const Tensor4& t)
{
	double retval = 0.0;
	for (int i = 0; i < t.w; i++)
		for (int j = 0; j < t.x; j++)
			for (int k = 0; k < t.y; k++)
				for (int l = 0; l < t.z; l++)
					retval += t(i, j, k, l) * t(i, j, k, l);
	
	return std::sqrt(retval);
}


S8EvenTensor4::S8EvenTensor4(int _w)
{
  // Resize the data matrix
  w = _w;
  resize(w);
}

S8EvenTensor4::S8EvenTensor4(int _w, double val)
{
  // Resize and assign the data matrix
	w = _w;
	assign(w, val);
}

S8EvenTensor4::S8EvenTensor4(const S8EvenTensor4& other)
{
  w = other.w;
  int size = w*(w+1)*(w+2)*(3*w+1) / 24;
  data.resize(size);
  updateMults();
  for (int i = 0; i < w; i++){
    for (int j = 0; j <= i; j++){
      for (int k = 0; k <= i; k++){
        for (int l = 0; l <= k; l++){
          data[imults[i] + j*jkmults[i+1] + jkmults[k] + l] = other(i, j, k, l);
        }
      } 
    }
  }
}  

void S8EvenTensor4::resize(int _w)
{
	w = _w;
 	int size = w*(w+1)*(w+2)*(3*w+1) / 24;
  	data.resize(size);
	updateMults();
}

void S8EvenTensor4::assign(int _w, double val)
{
	w = _w;
  	int size = w*(w+1)*(w+2)*(3*w+1) / 24;
  	data.resize(size);
  	std::fill(data.begin(), data.end(), val);
	updateMults();
}

void S8EvenTensor4::print() const 
{
  for (int i = 0; i < w; i++){
    for (int j = 0; j <= i; j++){
      for (int k = 0; k <= i; k++){
	for (int l = 0; l <= k; l++){
	  std::cout << i << " " << j << " " << k << " " << l << "   " << data[imults[i] + j*jkmults[i+1] + jkmults[k] + l] << "\n";
	}
      }
    }
  }
  std::cout << "\n\n";
}

double& S8EvenTensor4::operator()(int i, int j, int k, int l)
{
	int I = i > j ? i : j;
	int J = i > j ? j : i;
	int K = k > l ? k : l;
	int L = k > l ? l : k;
	if (K > I) {
		int tmp = I;
		I = K;
		K = tmp;
		tmp = J;
		J = L;
		L = tmp;
	}
  	return data[imults[I] + J*jkmults[I+1] + jkmults[K] + L];
}

double S8EvenTensor4::operator()(int i, int j, int k, int l) const
{
	int I = i > j ? i : j;
	int J = i > j ? j : i;
	int K = k > l ? k : l;
	int L = k > l ? l : k;
	if (K > I) {
		int tmp = I;
		I = K;
		K = tmp;
		tmp = J;
		J = L;
		L = tmp;
	}
  	return data[imults[I] + J*jkmults[I+1] + jkmults[K] + L];
}

S8EvenTensor4& S8EvenTensor4::operator=(const S8EvenTensor4& other)
{
  resize(other.getW());
  for (int i = 0; i < w; i++){
    for (int j = 0; j <= i; j++){
      for (int k = 0; k <= i; k++){
	for (int l = 0; l <= k; l++){
	  data[imults[i] + j*jkmults[i+1] + jkmults[k] + l] = other(i, j, k, l);
	}
      }
    }
  }
  return *this;
}

S8EvenTensor4 S8EvenTensor4::operator+(const S8EvenTensor4& other) const
{
	S8EvenTensor4 retVal(w);
	for (int i = 0; i < w; i++){
		for (int j = 0; j <= i; j++){
			for (int k = 0; k <= i; k++){
				for (int l = 0; l <= k; l++){
					retVal(i, j, k, l) = data[imults[i] + j*jkmults[i+1] + jkmults[k] + l] + other(i, j, k, l);
				}
			}
		}
	}
	return retVal;
}

S8EvenTensor4 S8EvenTensor4::operator-(const S8EvenTensor4& other) const
{
	S8EvenTensor4 retVal(w);
	for (int i = 0; i < w; i++){
		for (int j = 0; j <= i; j++){
			for (int k = 0; k <= i; k++){
				for (int l = 0; l <= k; l++){
					retVal(i, j, k, l) = data[imults[i] + j*jkmults[i+1] + jkmults[k] + l] - other(i, j, k, l);
				}
			}
		}
	}
	return retVal;
}

void S8EvenTensor4::updateMults() {
	
	imults.clear(); jkmults.clear();
	for (int i = 0; i < w; i++) {
		int val = ( i * (i+1) ) / 2;
		jkmults.push_back(val);
		val *= (i+2) * (3*i + 1);
		imults.push_back(val/12);
	}
	jkmults.push_back( (w * (w+1))/2 );
	
}

S8OddTensor4::S8OddTensor4(int N) : S8EvenTensor4(N) { }
S8OddTensor4::S8OddTensor4(int N, double val) : S8EvenTensor4(N, val) { }
S8OddTensor4::S8OddTensor4(const S8OddTensor4& other) : S8EvenTensor4(other) { } 

double& S8OddTensor4::operator()(int i, int j, int k, int l) {
	int I = i > j ? i : j;
	int J = i > j ? j : i;
	int K = k > l ? k : l;
	int L = k > l ? l : k;
	if (K > I) {
		int tmp = I;
		I = K;
		K = tmp;
		tmp = J;
		J = L;
		L = tmp;
	}
  	return data[imults[I] + J*jkmults[I+1] + jkmults[K] + L];
}

double S8OddTensor4::operator()(int i, int j, int k, int l) const {
	int parity = 0;
	int I = i;
	int J = j;
	int K = k;
	int L = l;
	if (j > i) { I = j; J = i; parity++; }
	if (l > k) { K = l; L = k; parity++; }
	if (K > I) {
		int tmp = I;
		I = K;
		K = tmp;
		tmp = J;
		J = L;
		L = tmp;
	}
  	return (1 - 2*(parity%2)) * data[imults[I] + J*jkmults[I+1] + jkmults[K] + L];
}