#include "InternalCoords.hpp"
#include "mathutil.hpp"
#include <cmath>

void Rotator::calcAxis() {
	// Assumes molecule is linear 
	int natoms = atoms.size();
	if (natoms > 1) {
		Vector vy = x0.row(atoms[natoms-1]) - x0.row(atoms[0]).normalized(); 
		
		double theta_x = acos(vy[0]);
		double theta_y = acos(vy[1]);
		double theta_z = acos(vy[2]); 
		
		Vector e0 = Vector::Zero(3); 
		if (theta_x < theta_y && theta_x < theta_z) e0[0] = 1.0;
		else if (theta_y < theta_x && theta_y < theta_z) e0[1] = 1.0;
		else e0[2] = 1.0; 
		
		axis = vy.cross(e0).normalized(); 
	}
}

Vector Rotator::value(const Matrix& xyz) {
	Vector rVal;
	if ((xyz - valXYZ).norm() < 1e-12) rVal = stored_value; 
	else {
		int natoms = atoms.size();
		
		Matrix subXYZ = Matrix::Zero(natoms, 3);
		Matrix subX0 = subXYZ;
		for (int i = 0; i < natoms; i++) {
			subXYZ.row(i) = xyz.row(atoms[i]); 
			subX0.row(i) = x0.row(atoms[i]); 
		}
		
		Vector xmean(3), ymean(3);
		for (int i = 0; i < 3; i++) {
			xmean[i] = subXYZ.row(i).sum() / ((double) natoms);
			ymean[i] = subX0.row(i).sum() / ((double) natoms);
		}
		 
		if (!linear && is_linear(subXYZ, subX0)) linear = true; 
		
		if (linear) {
			
			Vector vx = subXYZ.row(natoms-1) - subXYZ.row(0);
			Vector vy = subX0.row(natoms-1) - subX0.row(0); 
			
			calcAxis();
			
			Vector ev = vx.normalized();
			dot2 = ev.transpose() * axis; 
			dot2 *= dot2; 
			
			// Add dummy atom one Bohr from molecular center
			Vector xdum = vx.cross(axis).normalized();
			Vector ydum = vy.cross(axis).normalized();
			
			subXYZ.conservativeResize(subXYZ.rows() + 1, Eigen::NoChange); 
			subX0.conservativeResize(subX0.rows() + 1, Eigen::NoChange); 
			subXYZ.row(natoms) = xdum+xmean;
			subX0.row(natoms) = ydum+ymean; 
		}
		
		rVal = get_exp_map(subXYZ, subX0);
		norm = rVal.norm(); 
		valXYZ = xyz; 
		stored_value = rVal;  
	}
	return rVal; 
}

Matrix Rotator::derivative(const Matrix& xyz) {
	
}

