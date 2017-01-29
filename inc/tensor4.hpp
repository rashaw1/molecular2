/*
 *
 *   PURPOSE: To declare a class for a tensor of rank 4, mostly for storage purposes.
 *            It's a class instead of a struct, as I'm currently undecided about
 *            how superficial/deep it will be.
 *
 */

#ifndef TENSOR4HEADERDEF
#define TENSOR4HEADERDEF

#include <vector>

class Tensor4
{
protected:
  std::vector<double> data;
  int w, x, y, z;
public:
  Tensor4() { } 
  Tensor4(int a, int b, int c, int d);
  Tensor4(int a, int b, int c, int d, double val);
  Tensor4(const Tensor4& other);
	int getW() const { return w; }
	int getX() const { return x; }
	int getY() const { return y; }
	int getZ() const { return z; }
	void resize(int a, int b, int c, int d);
  void assign(int a, int b, int c, int d, double val);
  virtual void print() const;
  virtual double& operator()(int i, int j, int k, int l);
  virtual double operator()(int i, int j, int k, int l) const;
  Tensor4& operator=(const Tensor4& other);
  Tensor4 operator+(const Tensor4& other) const;
  Tensor4 operator-(const Tensor4& other) const;
  friend double fnorm(const Tensor4& t);
};

class S8EvenTensor4 : public Tensor4 
{
protected:
	std::vector<int> imults;
	std::vector<int> jkmults;
	void updateMults();
public:
    S8EvenTensor4() { } 
    S8EvenTensor4(int N);
    S8EvenTensor4(int N, double val);
   	S8EvenTensor4(const S8EvenTensor4& other);	
  	void resize(int N);
    void assign(int N, double val);
    virtual void print() const;
    virtual double& operator()(int i, int j, int k, int l);
    virtual double operator()(int i, int j, int k, int l) const;
    S8EvenTensor4& operator=(const S8EvenTensor4& other);
    S8EvenTensor4 operator+(const S8EvenTensor4& other) const;
    S8EvenTensor4 operator-(const S8EvenTensor4& other) const;
};

class S8OddTensor4 : public S8EvenTensor4 
{
public:
    S8OddTensor4() { } 
    S8OddTensor4(int N);
    S8OddTensor4(int N, double val);
   	S8OddTensor4(const S8OddTensor4& other);	
    virtual double& operator()(int i, int j, int k, int l);
    virtual double operator()(int i, int j, int k, int l) const;
};

#endif
