#include "optimiser.hpp"
#include "logger.hpp"
#include "scf.hpp"
#include "atom.hpp"
#include "integrals.hpp"
#include <vector>

Vector quadratic(Matrix& hessian, Vector& gradient) {
	Matrix H = hessian.selfadjointView<Eigen::Upper>();
	int n = H.rows() + 1;
	Matrix AH = Matrix::Zero(n, n);
	AH.block(0, 0, n-1, n-1) = H;
	AH.block(n-1, 0, 1, n-1) = -gradient.transpose();
	AH.block(0, n-1, n-1, 1) = -gradient;
	
	// Rescale augmented Hessian
	double alpha = 1.0;
	double tval = 0.0; 
	bool converged = false;
	int iter = -1;
	Matrix SAH(n, n); 
	EigenSolver hes(H); 
	Matrix hesvects = hes.eigenvectors();
	Vector h = hes.eigenvalues();
	Vector dq; 
	while (!converged && iter < 10) {
		++iter;
		if (iter == 10) alpha = 1.0; 
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n-1; j++) SAH(j, i) = AH(j, i) / alpha;
			SAH(n-1, i) = AH(n-1, i);
		}
		
		EigenSolver es(SAH);
		SAH = es.eigenvectors();
		for (int i = 0; i < n; i++) {
			tval = fabs( SAH.block(0, i, n-1, 1).lpNorm<Eigen::Infinity>() / SAH(i, n-1) );
			if ( tval < 1e5 ) {
				for (int j = 0; j < n; j++) SAH(i, j) /= SAH(i, n-1); 
			}
		}
		
		tval = 1e6; 
		int root = 0; 
		while (tval > 1e5) {
			tval = fabs( SAH.block(0, root, n-1, 1).lpNorm<Eigen::Infinity>() / SAH(root, n-1));
			if (tval > 1e5) root++; 
			if (root == n) root = 0; 
		}
		
		dq = SAH.block(0, root, n-1, 1); 
		double dqnorm = dq.norm(); 
		
		if (dqnorm < 0.4 && fabs(alpha) < 1e8) converged = true; 
		
		double lambda = -1.0 * gradient.dot(dq); 
		double sum = 0.0;
		for (int i = 0; i < n-1; i++) {
			double sq1 = hesvects.col(i).dot(gradient);
			double sq2 =  h(i) - lambda * alpha; 
			sum += sq1 * sq1 / (sq2 * sq2 * sq2);
		}
		
		tval = 2 * lambda * sum / (1.0 + alpha*dqnorm*dqnorm); 
		alpha += 2 * (0.4 * dqnorm - dqnorm*dqnorm) / tval; 
		
			std::cout << tval << " " << dqnorm << " " << alpha << "\n"; 
	}
	if (!converged) std::cout << "Not converged\n";
	return dq; 
}

void quadratic_scf(Command& cmd, Molecule& mol, Fock& blah)
{
	Logger& log = mol.control.log;
	log.title("QUADRATIC GEOMETRY OPTIMIZATION"); 
	
	bool unrestricted = (mol.getMultiplicity() > 1) || (mol.getNel()%2 != 0);
	
	double grad_norm = 1.0;
	double old_e = 0.0;
	double energy = 0.0; 
	int iter = 0;
	bool converged = false;
	int MAXITER = cmd.get_option<int>("maxiter");
	double CONVERGE = cmd.get_option<double>("converge");
	
	while (!converged && iter < MAXITER) {
		
		IntegralEngine ints(mol, false);
		Fock f(cmd, ints, mol);
		SCF hf(cmd, mol, f); 
		if (unrestricted) hf.uhf(false);
		else hf.rhf(false);
		energy = hf.getEnergy();
		
		std::vector<Atom> atomlist; 
		for (int i = 0; i < mol.getNAtoms(); i++) atomlist.push_back(mol.getAtom(i));
		f.getDens() = 0.5 * f.getDens();
		f.compute_forces(atomlist, mol.getNel()/2);
		Matrix g = f.getForces().transpose(); 
		Vector gradient(Eigen::Map<Vector>(g.data(), g.cols()*g.rows()));
		grad_norm = gradient.norm();
		converged = grad_norm < CONVERGE;
		log.iteration(iter++, energy, fabs(energy-old_e), grad_norm);
		old_e = energy; 
		
		if (!converged) {	
			f.compute_hessian(atomlist, mol.getNel()/2);
			Vector delta_geom = quadratic(f.getHessian(), gradient);
			int offset = 0; 
			for (int i = 0; i < mol.getNAtoms(); i++) {
				mol.getAtom(i).translate(delta_geom[offset], delta_geom[offset+1], delta_geom[offset+2]);
				mol.updateBasisPositions();
				offset+=3; 
			} 
		}
	}
	
	if (converged) {
		log.print("Converged geometry:\n");
		for (int i = 0; i <mol.getNAtoms(); i++) log.print(mol.getAtom(i));
		log.result("HF Energy", energy, "Hartree");
	} else {
		log.result("Geometry optimisation failed to converge.");
	}
}