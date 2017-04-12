#include "optimiser.hpp"
#include "logger.hpp"
#include "scf.hpp"
#include "atom.hpp"
#include "integrals.hpp"
#include <vector>
#include <Eigen/Core>

double rfo(Vector& dx, Vector& grad, Matrix& hessian, double alpha, double stepsize) {
	
	Vector step; 
	std::pair<double, double> results = rfo_prime(step, grad, hessian, alpha); 
	Vector m_dy; 
	double m_expect; 
	double step_norm = step.norm(); 
	if (step_norm > stepsize) {
		double v = alpha; 
		int niter = 0; 
		double ndy_last = 0.0; 
		
		double m_ndy = step_norm; 
		m_dy = step; 
		m_expect = results.first; 
		
		while(niter < 100) {
			v += (1-step_norm / stepsize) * step_norm / results.second; 
			results = rfo_prime(step, grad, hessian, v); 
			step_norm = step.norm(); 
			
			if (fabs(step_norm - stepsize)/stepsize < 0.001) {
				m_dy = step;
				m_expect = results.first;
				break; 
			}
			else if (niter > 10 && fabs(ndy_last - step_norm)/step_norm < 0.001) {
				m_dy = step;
				m_expect = results.first;
				break; 	
			}
			
			niter++; 
			ndy_last = step_norm; 
			if (step_norm < m_ndy) {
				m_ndy = step_norm;
				m_dy = step;
				m_expect = results.first; 
			}	
		}
	} else {
		m_dy = step; 
		m_expect = results.first; 
	}
	
	dx = m_dy; 
	return m_expect; 
}

double trust_newton(Vector& dx, Vector& grad, Matrix& hessian, double stepsize) {
	EigenSolver hsolve(hessian); 
	Vector gtilde = hsolve.eigenvectors().transpose() * grad; 
	Vector step = Vector::Zero(hessian.rows());
	for (int i = 0; i < hessian.rows(); i++) {
		double eival = hsolve.eigenvalues()[i]; 
		if (fabs(eival) > 1e-12)
			step[i] = - gtilde[i] / eival; 
	}
	
	dx = hsolve.eigenvectors() * step; 
	
	Matrix ident = Matrix::Identity(hessian.rows(), hessian.cols()); 
	if (dx.norm() > stepsize) {
		double eimin = hsolve.eigenvalues()[0]; 
		double factor = eimin < 0 ? 1.0 : 0.0; 
		double mu = (factor + 0.9)*eimin;
		
		Matrix hshift = hsolve.eigenvalues().asDiagonal();
		for (int i = 0; i < hshift.rows(); i++) {
			double hii = hshift(i, i) - mu;
			hshift(i, i) = 1.0 / (hii * hii); 
		}
		
		double delta = gtilde.transpose() * hshift * gtilde;
		delta = sqrt(delta) - stepsize; 
		
		for (int i = 1; i < 9; i++) {
			double new_mu = (factor + 0.1) * i * eimin; 

			for (int i = 0; i < hshift.rows(); i++) {
				double hii = hsolve.eigenvalues()[i] - new_mu;
				hshift(i, i) = 1.0 / (hii * hii); 
			}
		
			double new_delta = gtilde.transpose() * hshift * gtilde;
			new_delta = sqrt(new_delta) - stepsize; 
			
			if (fabs(new_delta) < fabs(delta)) {
				mu = new_mu;
				delta = new_delta;
			}
		}
		
		for (int i = 0; i < hessian.rows(); i++) {
			double eival = hsolve.eigenvalues()[i]; 
			step[i] = - gtilde[i] / (eival - mu); 
		}
		
	}
	
	dx = hsolve.eigenvectors() * step; 
	double expect = -grad.dot(dx) + dx.transpose() * hessian * dx; 
	
	return expect; 
}

std::pair<double, double> rfo_prime(Vector& dy, Vector &grad, Matrix& hessian, double alpha) {
	
	int ndim = hessian.rows();
	int eival = 2; 
	
	Matrix augH = Matrix::Zero(ndim+1, ndim+1); 
	augH.block(0, 0, ndim, ndim) = hessian; 
	augH.block(ndim, 0, 1, ndim) = -grad.transpose();
	augH.block(0, ndim, ndim, 1) = -grad; 
	
	Matrix B = alpha * Matrix::Identity(ndim+1, ndim+1);
	B(ndim, ndim) = 1.0; 

	GeneralizedEigenSolver solver(augH, B); 
	double lmin = solver.eigenvalues()[eival]; 
	Vector step = solver.eigenvectors().block(0, eival, ndim, 1);
	std::cout << solver.eigenvectors() << std::endl << std::endl;
	double vmin = solver.eigenvectors()(ndim, eival); 
	if (fabs(vmin) > 1e-12)
		step /= vmin;
	
	double nu = alpha * lmin; 
	EigenSolver hsolve(hessian); 
	const Matrix& hvec = hsolve.eigenvectors();
	const Vector& heig = hsolve.eigenvalues(); 
	double dyprime2 = 0.0;
	double dy2 = 0.0; 
	double dotp, denom; 
	for (int i = 0; i < hessian.cols(); i++) {
		dotp = hvec.col(i).dot(grad); 
		dotp *= dotp; 
		denom = heig[i] - nu; 
		dyprime2 += dotp / (denom * denom * denom);
		dy2 += dotp / (denom * denom); 
	}
	double expect = 1 + alpha * step.dot(step); 
	dyprime2 *= 2.0 * lmin / expect; 
	expect *= 0.5 * lmin; 
	double dyprime = 0.5 * dyprime2 / sqrt(dy2); 
	
	dy = step; 
	std::pair<double, double> results = {expect, dyprime}; 
	return results; 
	
}



void RHFOptimiser::optimise() {
	Logger& log = mol->control->log;
	log.title("GEOMETRY OPTIMIZATION"); 
	int MAXITER = cmd.get_option<int>("maxiter");
	double CONVERGE = cmd.get_option<double>("converge");  
	
	Vector dx = Vector::Zero(3*mol->getNActiveAtoms()); 
	
	bool converged = false;
	Vector grad;
	Matrix hessian; 
	calc.iter = 0;
	double trust = cmd.get_option<double>("trust"); 
	double expect, trust_ratio; 
	while (!converged && calc.iter < MAXITER) {
		calc(dx, grad, hessian); 
		if (calc.iter > 1) {
		
			trust_ratio = calc.delta_e / expect; 
			
			double stepnorm = dx.norm();
			if (trust_ratio > 0.75 && 1.25*stepnorm > trust) trust *= 2.0;
			else if (trust_ratio < 0.25) trust = 0.25*stepnorm; 
		}
		
		converged = calc.grad_norm < CONVERGE; 
		if (!converged)
			expect = trust_newton(dx, grad, hessian, trust); 
	}
	
	if (calc.iter < MAXITER) {
		log.print("Converged geometry:\n");
		for (int i = 0; i <mol->getNAtoms(); i++) log.print(mol->getAtom(i));
		log.result("HF Energy", calc.energy, "Hartree");
		
		if (cmd.get_option<bool>("freq")) frequencies(hessian); 
		
	} else {
		log.result("Geometry optimisation failed to converge.");
	}
}

void RHFOptimiser::frequencies(Matrix& hessian) {
	
	int natom = hessian.rows(); 
	
	// form mass-weighted hessian
	int row = 0; 
	double mi, mj; 
	Matrix mass_hessian = hessian; 
	
	std::vector<int> activeAtoms = mol->getActiveList(); 
	for (int i : activeAtoms) {
		mi = mol->getAtom(i).getMass(); 
		mi = sqrt(mi); 
		
		for (int xyz = 0; xyz < 3; xyz++) {
			
			int col = 0; 
			for (int j : activeAtoms) {
				mj = mol->getAtom(j).getMass();
				mj = sqrt(mj); 
				
				for (int abc = 0; abc < 3; abc++) {
					mass_hessian(row, col) /= mi * mj; 
					col++;	
				}
			}
			
			row++;
		}
	}
	
	EigenSolver hsolve(mass_hessian); 
	Vector freqs = hsolve.eigenvalues(); 
	Matrix modes = hsolve.eigenvectors(); 
	
	int nimaginary = 0; 
	for (int i = 0; i < hessian.rows(); i++) {
		if (freqs[i] < -1e-4) {
			freqs[i] = -sqrt(-freqs[i]); 
			nimaginary++; 
		} else if (freqs[i] < 1e-4) {
			freqs[i] = 0.0; 
		} else freqs[i] = sqrt(freqs[i]); 
		freqs[i] /= (2.0 * M_PI);  
		freqs[i] *= Logger::TOWAVENUMBERS; 
	}
	
	mol->control->log.frequencies(freqs, modes, cmd.get_option<bool>("modes"));  
	
}

double RHFCalculator::operator()(const Vector& dx, Vector& grad, Matrix& hessian) {
		
	int offset = 0;
	int natoms = mol->getNAtoms();
	
	std::vector<int> activeAtoms = mol->getActiveList(); 
	for (int i : activeAtoms) {
		Atom& a = mol->getAtom(i); 
		a.translate(dx[offset], dx[offset+1], dx[offset+2]); 
		offset += 3; 
	}
	mol->updateBasisPositions(); 
	mol->calcEnuc(); 
		
	IntegralEngine ints(mol, false);
	Fock f(cmd, ints, mol);
	SCF hf(cmd, mol, f); 
	hf.rhf(false); 
	energy = hf.getEnergy();
		
	std::vector<Atom> atomlist; 
	for (int i = 0; i < natoms; i++) atomlist.push_back(mol->getAtom(i));
	
	f.compute_forces(atomlist, mol->getNel()/2);
	Matrix g = f.getForces().transpose(); 
	Vector full_grad = Eigen::Map<Vector>(g.data(), g.cols()*g.rows());
	
	f.compute_hessian_numeric(atomlist, mol->getNel()/2, cmd); 
	if (natoms == activeAtoms.size()) {
		hessian = f.getHessian(); 
		grad = full_grad; 
	} else {
		Matrix& full_hessian = f.getHessian(); 
	
		// Restrict to only active atoms
		int nactive = activeAtoms.size(); 
		grad = Vector::Zero(3*nactive);
		Matrix temphessian = Matrix::Zero(3*nactive, 3*natoms);
		int row = 0; 
		for (int i : activeAtoms) {
			grad.segment(row, 3) = full_grad.segment(3*i, 3); 
			temphessian.block(row, 0, 3, 3*natoms) = full_hessian.block(3*i, 0, 3, 3*natoms); 
			row += 3; 
		}
	
		row = 0; 
		hessian = Matrix::Zero(3*nactive, 3*nactive); 
		for (int j : activeAtoms) {
			hessian.block(0, row, 3*nactive, 3) = temphessian.block(0, 3*j, 3*nactive, 3); 
			row += 3; 
		}
	
	}
	
	delta_e = energy - old_e; 
	grad_norm = grad.norm(); 
	mol->control->log.iteration(iter++, energy, delta_e, grad_norm);
	old_e = energy; 
	
	return energy; 
	
}

