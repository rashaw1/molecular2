#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/Cholesky>
#include <Eigen/Eigenvalues>
#include "tensor4.hpp"

using Matrix = Eigen::MatrixXd; 
using SparseMatrix = Eigen::SparseMatrix<double>; 
using Vector = Eigen::VectorXd; 
using iVector = Eigen::VectorXi; 
using EigenSolver = Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd>; 
using GeneralizedEigenSolver = Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd>;
using LLT = Eigen::LLT<Eigen::MatrixXd>; 
using FullLU = Eigen::FullPivLU<Eigen::MatrixXd>;

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> EMatrix;
typedef Eigen::DiagonalMatrix<double, Eigen::Dynamic, Eigen::Dynamic> DiagonalMatrix;

Vector tensorToVector(const Tensor4& t, Tensor4::TYPE symm);
double angle(const Vector& u, const Vector& w);
Vector cross(const Vector& u, const Vector& w);