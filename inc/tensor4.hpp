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
#include <string>

class Tensor4
{
protected:
  std::vector<double> data;
  int w, x, y, z;
public:
	
	enum TYPE {
		DEFAULT,
		ODD_8,
		EVEN_8,
		ODD_4,
		EVEN_4
	};

  Tensor4() { w = x = y = z = 0; } 
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
  Tensor4 operator+() const;
  Tensor4 operator-() const;
  Tensor4 operator+(const Tensor4& other) const;
  Tensor4 operator-(const Tensor4& other) const;
  Tensor4 operator*=(const double& scalar) const; 
  
  size_t data_size; 
  std::string file_name; 
  
  virtual int writeToFile(std::string filename); 
  virtual int readFromFile();
  
  friend double fnorm(const Tensor4& t);
};

inline Tensor4 operator*(const Tensor4& t, const double& scalar)
{
	int w = t.getW(); int x = t.getX(); int y = t.getY(); int z = t.getZ();
	Tensor4 retVal(w, x, y, z);
	for (int i = 0; i < w; i++)
		for (int j = 0; j < x; j++)
			for (int k = 0; k < y; k++)
				for (int l = 0; l < z; l++)
					retVal(i, j, k, l) = t(i, j, k, l) * scalar;
	return retVal;
}
inline Tensor4 operator*(const double& scalar, const Tensor4& t)
{
	int w = t.getW(); int x = t.getX(); int y = t.getY(); int z = t.getZ();
	Tensor4 retVal(w, x, y, z);
	for (int i = 0; i < w; i++)
		for (int j = 0; j < x; j++)
			for (int k = 0; k < y; k++)
				for (int l = 0; l < z; l++)
					retVal(i, j, k, l) = t(i, j, k, l) * scalar;
	return retVal;
}

class S8EvenTensor4 : public Tensor4 
{
protected:
	std::vector<int> imults;
	std::vector<int> jkmults;
	void updateMults();
public:
    S8EvenTensor4() { w = 0; } 
    S8EvenTensor4(int N);
    S8EvenTensor4(int N, double val);
   	S8EvenTensor4(const S8EvenTensor4& other);	
  	void resize(int N);
    void assign(int N, double val);
    virtual void print() const;
    virtual double& operator()(int i, int j, int k, int l);
    virtual double operator()(int i, int j, int k, int l) const;
    S8EvenTensor4& operator=(const S8EvenTensor4& other);
    S8EvenTensor4 operator+() const;
    S8EvenTensor4 operator-() const;
    S8EvenTensor4 operator+(const S8EvenTensor4& other) const;
    S8EvenTensor4 operator-(const S8EvenTensor4& other) const;
	S8EvenTensor4 operator*=(const double& scalar) const;
	friend double fnorm(const S8EvenTensor4& t);
};

inline S8EvenTensor4 operator*(const S8EvenTensor4& t, const double& scalar)
{
	int w = t.getW();
	S8EvenTensor4 retVal(w);
	for (int i = 0; i < w; i++)
		for (int j = 0; j <= i; j++)
			for (int k = 0; k <= i; k++)
				for (int l = 0; l <= k; l++)
					retVal(i, j, k, l) = t(i, j, k, l) * scalar;
	return retVal;
}
inline S8EvenTensor4 operator*(const double& scalar, const S8EvenTensor4& t)
{
	int w = t.getW();
	S8EvenTensor4 retVal(w);
	for (int i = 0; i < w; i++)
		for (int j = 0; j <= i; j++)
			for (int k = 0; k <= i; k++)
				for (int l = 0; l <= k; l++)
					retVal(i, j, k, l) = t(i, j, k, l) * scalar;
	return retVal;
}

class S8OddTensor4 : public S8EvenTensor4 
{
public:
    S8OddTensor4() { w = 0; } 
    S8OddTensor4(int N);
    S8OddTensor4(int N, double val);
   	S8OddTensor4(const S8OddTensor4& other);	
	virtual void set(int i, int j, int k, int l, double value);
    virtual double operator()(int i, int j, int k, int l) const;
};

class S4OddTensor4 : public Tensor4 
{
protected:
	std::vector<int> mults;
	void updateMults();
public:
    S4OddTensor4() { w = 0; x=0; } 
    S4OddTensor4(int n, int m);
    S4OddTensor4(int n, int m, double val);
   	S4OddTensor4(const S4OddTensor4& other);	
  	void resize(int n, int m);
    void assign(int n, int m, double val);
    virtual void print() const;
    virtual void set(int i, int j, int k, int l, double value);
    virtual double operator()(int i, int j, int k, int l) const;
    S4OddTensor4& operator=(const S4OddTensor4& other);
    S4OddTensor4 operator+() const;
    S4OddTensor4 operator-() const;
    S4OddTensor4 operator+(const S4OddTensor4& other) const;
    S4OddTensor4 operator-(const S4OddTensor4& other) const;
	S4OddTensor4 operator*=(const double& scalar) const;
	friend double fnorm(const S4OddTensor4& t);
};

inline S4OddTensor4 operator*(const S4OddTensor4& t, const double& scalar)
{
	int w = t.getW(); int x = t.getX();
	S4OddTensor4 retVal(w, x);
	for (int i = 0; i < w; i++)
		for (int j = 0; j <= i; j++)
			for (int k = 0; k < x; k++)
				for (int l = 0; l <= k; l++)
					retVal.set(i, j, k, l, t(i, j, k, l) * scalar);
	return retVal;
}

inline S4OddTensor4 operator*(const double& scalar, const S4OddTensor4& t)
{
	int w = t.getW(); int x = t.getX();
	S4OddTensor4 retVal(w, x);
	for (int i = 0; i < w; i++)
		for (int j = 0; j <= i; j++)
			for (int k = 0; k < x; k++)
				for (int l = 0; l <= k; l++)
					retVal.set(i, j, k, l, t(i, j, k, l) * scalar);
	return retVal;
}

class S4EvenTensor4 : public S4OddTensor4 
{
public:
    S4EvenTensor4() { w = 0; x=0; } 
    S4EvenTensor4(int n, int m);
    S4EvenTensor4(int n, int m, double val);
   	S4EvenTensor4(const S4EvenTensor4& other);	
	virtual void set(int i, int j, int k, int l, double value);
    virtual double operator()(int i, int j, int k, int l) const;
};

#endif
