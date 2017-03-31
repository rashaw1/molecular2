#include "ecpint.hpp"
#include "Faddeeva.hpp"
#include "gshell.hpp"
#include <cmath>
#include <iostream>

const double ROOT_PI = 1.772453850905516;

void compute_base_integrals(int N_min, int N_max, double p, double o_root_p, double P1, double P2, double P1_2, double P2_2, double X1, double X2, double oP1, double oP2, double* values) {
	
	int hmax = N_max / 2;
	int hmin = (N_min + 1) / 2;
	int gmax = (N_max - 1) / 2;
	int gmin = N_min / 2;
	
	double P1_2k = 1.0;
	double P2_2k = 1.0;
	
	for (int k = 2; k < hmin; k++) {
		P1_2k *= P1_2;
		P2_2k *= P2_2;
	}
	
	double ck, dk, ek, val; 
	double C0 = o_root_p * ROOT_PI;
	for (int n = hmin; n <= hmax; n++) {
		ck = C0; 
		dk = P1_2k * X1;
		ek = P2_2k * X2; 
		val = ck * (dk - ek); 	
		
		for (int k = n - 1; k > 0; k--) {
			ck *= 2*k*(2*k - 1)*(n-k-0.5) / ((2*n - 2*k) * (2*n - 2*k - 1) * p);
			dk *= oP1;
			ek *= oP2; 
			val += ck * (dk - ek);
		}
		
		values[2*n - N_min] = val;
		
		P1_2k *= P1_2;
		P2_2k *= P2_2;
	}
	
	P1_2k = P1;
	P2_2k = P2;
	for (int k = 1; k < gmin; k++) {
		P1_2k *= P1_2;
		P2_2k *= P2_2;
	} 
	
	
	for (int n = gmin; n <= gmax; n++) {
		ck = C0; 
		dk = P1_2k * X1;
		ek = P2_2k * X2; 
		val = ck * (dk - ek);
		
		for (int k = n-1; k >0; k--) {
			ck *= 2*k*(2*k+1)*(n-k-0.5) / ((2*n-2*k) * (2*n - 1 - 2*k) * p);
			dk *= oP1; 
			ek *= oP2; 
			val += ck * (dk - ek);
		}
		
		values[2*n + 1 - N_min] = val;
		
		P1_2k *= P1_2;
		P2_2k *= P2_2; 
	}
	
}

void G000(ECP& U, GaussianShell& shellA, GaussianShell& shellB, ShellPairData& data, ThreeIndex<double> &values) {
	double value = 0.0;
	
	int npA = shellA.nprimitive();
	int npB = shellB.nprimitive();
	int npC = U.getN();

	double A = data.Am; 
	double B = data.Bm;
	double A2 = data.A2;
	double B2 = data.B2; 

	double zA, zB, zC, dA, dB, dC; 
	double p, P1_2, P2_2, aAbB; 
	for (int c = 0; c < npC; c++) {
		GaussianECP& g = U.getGaussian(c);
		
		if (g.l == 0) {
			zC = g.a;
			dC = g.d;
			
			for (int a = 0; a < npA; a++) {
				zA = shellA.exp(a);
				dA = shellA.coef(a);
				
				for (int b = 0; b < npB; b++) {
					zB = shellB.exp(b);
					dB = shellB.coef(b); 
					
					p = zA + zB + zC;
					P1_2 = (zA*A + zB*B) / p;
					P2_2 = (zB*B - zA*A) / p;
					P1_2 *= P1_2;
					P2_2 *= P2_2;
					aAbB = zA*A2 + zB*B2; 
					
					value += dA * dB * dC * (exp(p * P1_2 - aAbB) - exp(p*P2_2 - aAbB))/(zA*zB*A*B*sqrt(p)); 
				}
			}
		}
	}
	
	values(0, 0, 0) = 0.25 * M_PI * ROOT_PI * value;
}

void G001(ECP& U, GaussianShell& shellA, GaussianShell& shellB, ShellPairData& data, ThreeIndex<double>& values) {
	int npA = shellA.nprimitive();
	int npB = shellB.nprimitive();
	int npC = U.getN();

	double A = data.Am; 
	double B = data.Bm;
	double A2 = data.A2;
	double B2 = data.B2; 

	double zA, zB, zC, dA, dB, dC; 
	double p, P1, P2, P1_2, P2_2, X1, X2, aAbB; 
	double radial_value = 0.0;
	for (int c = 0; c < npC; c++) {
		GaussianECP& g = U.getGaussian(c);
		
		if (g.l == 1) {
			zC = g.a;
			dC = g.d;
			
			for (int a = 0; a < npA; a++) {
				zA = shellA.exp(a);
				dA = shellA.coef(a);
				
				for (int b = 0; b < npB; b++) {
					zB = shellB.exp(b);
					dB = shellB.coef(b); 
					
					p = zA + zB + zC;
					P1 = (zA*A + zB*B) / p;
					P2 = (zB*B - zA*A) / p;
					P1_2 = P1 * P1;
					P2_2 = P2 * P2;
					aAbB = zA*A2 + zB*B2; 
					X1 = exp(p*P1_2 - aAbB);
					X2 = exp(p*P2_2 - aAbB); 
					
					radial_value += dA * dB * dC * (p*(P1*X1 - P2*X2) - (zB*B + 0.5*p/(zB*B))*(X1 - X2))/(zA*zA*zB*A*A*B*sqrt(p)); 
				}
			}
		}
	}
	
	double xA = A > 0 ? data.A[2] / A : 0.0;
	double xB = B > 0 ? data.B[2] / B : 0.0;
	double phiA = atan2(data.A[1], data.A[0]);
	double phiB = atan2(data.B[1], data.B[0]);
	
	std::vector<double> fac = facArray(3);
	std::vector<double> dfac = dfacArray(3); 
	TwoIndex<double> SA = realSphericalHarmonics(1, xA, phiA, fac, dfac);
	TwoIndex<double> SB = realSphericalHarmonics(1, xB, phiB, fac, dfac);
	
	radial_value *= ROOT_PI * M_PI * M_PI;
	values(0, 0, 0) = radial_value * SA(1, 0) * SB(1, 0);
	values(0, 0, 1) = radial_value * SA(1, 1) * SB(1, 1);
	values(0, 0, 2) = radial_value * SA(1, 2) * SB(1, 2); 
}

void G002(ECP& U, GaussianShell& shellA, GaussianShell& shellB, ShellPairData& data, ThreeIndex<double>& values) {
	int npA = shellA.nprimitive();
	int npB = shellB.nprimitive();
	int npC = U.getN();

	double A = data.Am; 
	double B = data.Bm;
	double A2 = data.A2;
	double B2 = data.B2; 

	double zA, zB, zC, dA, dB, dC; 
	double p, P1, P2, P1_2, P2_2, X1, X2, aAbB, root_p, Kab, o2bB; 
	double radial_value = 0.0;
	for (int c = 0; c < npC; c++) {
		GaussianECP& g = U.getGaussian(c);
		
		if (g.l == 2) {
			zC = g.a;
			dC = g.d;
			
			for (int a = 0; a < npA; a++) {
				zA = shellA.exp(a);
				dA = shellA.coef(a);
				
				for (int b = 0; b < npB; b++) {
					zB = shellB.exp(b);
					dB = shellB.coef(b); 
					
					p = zA + zB + zC;
					root_p = sqrt(p);
					P1 = (zA*A + zB*B) / p;
					P2 = (zB*B - zA*A) / p;
					P1_2 = P1 * P1;
					P2_2 = P2 * P2;
					aAbB = zA*A2 + zB*B2; 
					Kab = 1.0 / (zA * A * zB * B);
					X1 = exp(p*P1_2 - aAbB) * Kab;
					X2 = exp(p*P2_2 - aAbB) * Kab; 
					
					double H4 = 0.5 * (X1 - X2) / p + P1_2 * X1 - P2_2 * X2;
					double G3 = P1 * X1 - P2 * X2;
					double H2 = X1 - X2;
					
					o2bB = 0.5 / (zB * B); 
					double temp = p*p*H4 - (2.0*p*zB*B + 3.0*p*p*o2bB)*G3 + (zB*B*zB*B + p + 3.0*p*p*o2bB*o2bB)*H2;
					temp /= (zA * A * zA * A * root_p); 
					
					radial_value += dA * dB * dC * temp; 
				}
			}
		}
	}
	
	double xA = A > 0 ? data.A[2] / A : 0.0;
	double xB = B > 0 ? data.B[2] / B : 0.0;
	double phiA = atan2(data.A[1], data.A[0]);
	double phiB = atan2(data.B[1], data.B[0]);
	
	std::vector<double> fac = facArray(5);
	std::vector<double> dfac = dfacArray(5); 
	TwoIndex<double> SA = realSphericalHarmonics(2, xA, phiA, fac, dfac);
	TwoIndex<double> SB = realSphericalHarmonics(2, xB, phiB, fac, dfac);
	
	radial_value *= ROOT_PI * M_PI * M_PI;
	values(0, 0, 0) = radial_value * SA(2, 0) * SB(2, 0);
	values(0, 0, 1) = radial_value * SA(2, 1) * SB(2, 1);
	values(0, 0, 2) = radial_value * SA(2, 2) * SB(2, 2); 
	values(0, 0, 3) = radial_value * SA(2, 3) * SB(2, 3);
	values(0, 0, 4) = radial_value * SA(2, 4) * SB(2, 4);
}

void G003(ECP& U, GaussianShell& shellA, GaussianShell& shellB, ShellPairData& data, ThreeIndex<double>& values) {
	int npA = shellA.nprimitive();
	int npB = shellB.nprimitive();
	int npC = U.getN();

	double A = data.Am; 
	double B = data.Bm;
	double A2 = data.A2;
	double B2 = data.B2; 

	double zA, zB, zC, dA, dB, dC; 
	double p, P1, P2, P1_2, P2_2, X1, X2, aAbB, root_p, Kab, J; 
	double radial_value = 0.0;
	for (int c = 0; c < npC; c++) {
		GaussianECP& g = U.getGaussian(c);
		
		if (g.l == 2) {
			zC = g.a;
			dC = g.d;
			
			for (int a = 0; a < npA; a++) {
				zA = shellA.exp(a);
				dA = shellA.coef(a);
				
				for (int b = 0; b < npB; b++) {
					zB = shellB.exp(b);
					dB = shellB.coef(b); 
					
					p = zA + zB + zC;
					root_p = sqrt(p);
					P1 = (zA*A + zB*B) / p;
					P2 = (zB*B - zA*A) / p;
					P1_2 = P1 * P1;
					P2_2 = P2 * P2;
					aAbB = zA*A2 + zB*B2; 
					Kab = ROOT_PI / (16.0*zA * A * zB * B*root_p);
					X1 = exp(p*P1_2 - aAbB) * Kab;
					X2 = exp(p*P2_2 - aAbB) * Kab; 
					
					double G5 = 1.5 * (P1 * X1 - P2 * X2) / p + (P1_2 * P1 * X1 - P2_2 * P2 * X2);
					double H4 = 0.5 * (X1 - X2) / p + P1_2 * X1 - P2_2 * X2;
					double G3 = P1 * X1 - P2 * X2;
					double H2 = X1 - X2;
					
					J = 0.5 / (zB * B); 
					double bB = zB * B; 
					double aA = zA * A;
					double temp = p*p*p*G5 - 3.0*p*p*(bB + 2.0*p*J)*H4;
					temp += 3.0*p*(bB*bB + 1.5*p + 5.0*p*p*J*J)*G3;
					temp -= (bB*bB*bB + 15.0*p*p*p*J*J*J + 1.5*bB*p + 4.5*p*p*J)*H2;
					temp /= (aA * aA * aA);
					
					temp = temp < 1e-12 ? 0.0 : temp;
					
					radial_value += dA * dB * dC * temp; 
				}
			}
		}
	}
	
	double xA = A > 0 ? data.A[2] / A : 0.0;
	double xB = B > 0 ? data.B[2] / B : 0.0;
	double phiA = atan2(data.A[1], data.A[0]);
	double phiB = atan2(data.B[1], data.B[0]);
	
	std::vector<double> fac = facArray(7);
	std::vector<double> dfac = dfacArray(7); 
	TwoIndex<double> SA = realSphericalHarmonics(3, xA, phiA, fac, dfac);
	TwoIndex<double> SB = realSphericalHarmonics(3, xB, phiB, fac, dfac);
	
	radial_value *= 16.0* M_PI * M_PI;
	values(0, 0, 0) = radial_value * SA(3, 0) * SB(3, 0);
	values(0, 0, 1) = radial_value * SA(3, 1) * SB(3, 1);
	values(0, 0, 2) = radial_value * SA(3, 2) * SB(3, 2); 
	values(0, 0, 3) = radial_value * SA(3, 3) * SB(3, 3);
	values(0, 0, 4) = radial_value * SA(3, 4) * SB(3, 4);
	values(0, 0, 5) = radial_value * SA(3, 5) * SB(3, 5);
	values(0, 0, 6) = radial_value * SA(3, 6) * SB(3, 6);
}