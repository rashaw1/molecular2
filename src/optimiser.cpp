#include "optimiser.hpp"
#include "logger.hpp"
#include "scf.hpp"
#include "atom.hpp"
#include "integrals.hpp"
#include <vector>

void optimise(Command& cmd, SharedMolecule mol) {
	Logger& log = mol->control->log;
	log.title("CONJUGATE GEOMETRY OPTIMIZATION"); 
	
	bool unrestricted = (mol->getMultiplicity() > 1) || (mol->getNel()%2 != 0);
	
	double grad_norm = 1.0;
	double old_e = 0.0;
	double energy = 0.0; 
	int iter = 0;
	bool converged = false;
	int MAXITER = cmd.get_option<int>("maxiter");
	double CONVERGE = cmd.get_option<double>("converge");
	
	// Get initial Hessian
	
	while (!converged && iter < MAXITER) {
		
		IntegralEngine ints(mol, false);
		Fock f(cmd, ints, mol);
		SCF hf(cmd, mol, f); 
		if (unrestricted) hf.uhf(false);
		else hf.rhf(false);
		energy = hf.getEnergy();
		
		std::vector<Atom> atomlist; 
		for (int i = 0; i < mol->getNAtoms(); i++) atomlist.push_back(mol->getAtom(i));
		f.getDens() = 0.5 * f.getDens();
		f.compute_forces(atomlist, mol->getNel()/2);
		Matrix g = f.getForces().transpose(); 
		Vector gradient(Eigen::Map<Vector>(g.data(), g.cols()*g.rows()));
		
		// Calculate internal coordinates 
		// Diagonalise hessian
		// Set v0 = 1 for RFO
		// Take an NR step (get_delta_prime)
		// inorm = norm of step
		// cnorm = getCartesianNorm
		// if cnorm > 1.1*trust: (reduce step length) {
		// Froot
		// iopt = brent_wiki
		// if IC.bork, ForceRebuild
		// if ForceRebuild {
		// if LastForce {
		// Switch to Cartesian coordinates }
		// recover}
		// dy,expect = trust_step
		// cnorm = getCartesianNorm }
		// Dot = dot(dy/norm(dy), gradient/norm(gradient))
		// bump = cnorm > 0.8*trust
		// Copy current variables to previous before updating
		// Update internal coordinates
		// Calculate energy and gradient
		// Add new coords and gradients to history
		// Check convergence
		// Adjust trust radius or reject step
		// Check coordinate system evern (N) steps
		// Update the hessian
		
/*		double delta_e = energy - old_e; 
		grad_norm = cg_solver.grad_norms[cg_solver.grad_norms.size()-1]; 
		log.iteration(iter++, energy, delta_e, grad_norm);
		old_e = energy; 
		
		converged = fabs(delta_e) < CONVERGE; 
		
		if (!converged) {	
			int offset = 0; 
			for (int i = 0; i < mol->getNAtoms(); i++) {
				mol->getAtom(i).translate(step[offset], step[offset+1], step[offset+2]);
				offset+=3; 
			} 
			mol->updateBasisPositions();
		}*/
	}
	
	if (converged) {
		log.print("Converged geometry:\n");
		for (int i = 0; i <mol->getNAtoms(); i++) log.print(mol->getAtom(i));
		log.result("HF Energy", energy, "Hartree");
	} else {
		log.result("Geometry optimisation failed to converge.");
	}
}

