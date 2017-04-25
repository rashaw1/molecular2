#include "eigen_wrapper.hpp"
#include <vector> 

Vector cphf_single_solver(Vector& Q, Matrix& V, double tol = 1e-6, int maxiter = 10); 
std::vector<Vector> cphf_group_solver(std::vector<Vector>& Q, Matrix& V, double tol = 1e-6, int maxiter = 10); 

