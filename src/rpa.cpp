
#include "rpa.hpp"
#include "eigen_wrapper.hpp"
#include <ctf.hpp>
#include <libint2.hpp>
#include "logger.hpp"

// Constructor
RPA::RPA(Command& c, Fock& _focker) : cmd(c), focker(_focker)
{
	N = focker.getCP().rows();
	nocc = focker.getMolecule()->getNel()/2;
	nvirt = N - nocc; 
	energy = 0.0;
}

void RPA::compute(bool print) {
	
	Logger& log = focker.getMolecule()->control->log; 
	
	auto &shells = focker.getMolecule()->getBasis().getIntShells();
	std::vector<CTF::Tensor<> > V;
	
	if (print) log.print("Transforming integrals"); 
	if (cmd.get_option<bool>("longrange")) {
		if(print) log.print("Using long-range potential, with mu = " + std::to_string(cmd.get_option<double>("mu"))); 
		longrange_eris(shells, V);
	} else { 
		full_eris(V); 
	}
	if (print) log.localTime();  
	
	if (print) log.print("\nForming excitation matrices"); 
	int64_t sz1, sz2, *i1, *i2;
	double *viajb, *vijab;
	V[0].read_local(&sz1, &i1, &viajb);
	V[1].read_local(&sz2, &i2, &vijab);  
	
	int dim = nvirt*nocc; 
	int nvnono = nocc*dim; 
	Matrix A = Matrix::Zero(dim, dim); 
	Matrix B = Matrix::Zero(dim, dim);
	Matrix K = Matrix::Zero(dim, dim); 
	
	Vector& eps = focker.getEps(); 
	bool sosex = cmd.get_option<bool>("sosex"); 
	
	int axis = nocc*(nocc+1);
	axis /= 2; 
	int ix1, ix2;  
	for (int i = 0; i < nocc; i++) {
		for (int a = 0; a < nvirt; a++) {
			ix1 = i*nvirt+a;
			
			for (int j = 0; j < nocc; j++) {
				for (int b = 0; b < nvirt; b++) {
					ix2 = j*nvirt + b;
					
					K(ix1, ix2) = 2.0 * viajb[i + a*nocc + j*dim + b*nvnono]; 
					A(ix1, ix2) = K(ix1, ix2);
					B(ix1, ix2) = K(ix1, ix2);
					if (sosex) {
						A(ix1, ix2) -= vijab[i + (j*(j+1))/2 + a*axis + axis*(b*(b+1))/2]; 
						B(ix1, ix2) -= viajb[j + a*nocc + i*dim + b*nvnono]; 
					}
					if (i==j && a==b) A (ix1, ix2) += eps[a+nocc] - eps[i]; 
				}
			}
		}
	}
	delete i1;
	delete i2;
	delete vijab;
	delete viajb; 
	if (print) log.localTime(); 
	
	if (print) log.print("\nSolving Riccatti equations");
	if (cmd.get_option<bool>("longrange")) {
		
		Matrix Ap = A; 
		for (int i = 0; i < dim ;i++) Ap(i, i) = 0.0; 
		
		Matrix T = -B; 
		
		double delta = 1.0;
		int iter = 0;
		Matrix newT; 
		while (delta > 1e-4 && iter < 15) {
			newT = -(B + T*Ap);  
			newT -= (T*B + Ap)*T; 
			
			for (int p = 0; p < dim; p++)
				for (int q = 0; q < dim; q++)
					newT(p, q) /= A(p, p) + A(q, q); 
			
			delta = (newT - T).norm(); 
			iter++; 
			T = newT; 
		}
		
		energy = 0.5 * (K*T).trace(); 
		
	} else {
		Matrix M = (A+B)*(A-B); 
		Eigen::EigenSolver<Matrix> solver(M);
		const Eigen::VectorXcd& evals = solver.eigenvalues(); 
		Vector omega = Vector::Zero(evals.size()); 
	
		energy = 0.0; 
		for (int i = 0; i < omega.size(); i++) { 
			std::complex<double> ei = evals[i]; 
			omega[i] = std::sqrt(ei.real()); 
			energy += omega[i] - A(i, i); 
		}	
		if(print) log.print("\nRPA excitation energies: ");
		if(print) log.print(omega); 
	
		energy /= 2.0;
		if (sosex) energy /= 2.0; 
	}
	if (print) log.localTime(); 
}

void RPA::longrange_eris(const std::vector<libint2::Shell>& shells, std::vector<CTF::Tensor<> >& moInts) {
	using libint2::Shell;
	using libint2::Engine;
	using libint2::Operator;

	IntegralEngine& integrals = focker.getIntegrals(); 

	const auto n = integrals.nbasis(shells);
	int l2 = (n*(n+1))/2;
	
	int aoshape[4] = {N, N, N, N};  
	int iajbshape[4] = {nocc, nvirt, nocc, nvirt};
	int ijabshape[4] = {nocc, nocc, nvirt, nvirt}; 
	int sysy[4] = {SY, NS, SY, NS}; 
	int nsns[4] = {NS, NS, NS, NS};
	int syns[4] = {SY, NS, NS, NS}; 
	int nssy[4] = {NS, NS, SY, NS};
	
	CTF::Tensor<> aoInts(4, aoshape, sysy, dw); 
	
	int64_t sz, *ix; 
	double *values; 
	aoInts.read_local(&sz, &ix, &values); 
	
	double mu = cmd.get_option<double>("mu"); 
	Engine engine(Operator::erf_coulomb, integrals.max_nprim(shells), integrals.max_l(shells), 0);
	engine.set_params(mu); 
	
	auto shell2bf = integrals.map_shell_to_basis_function(shells);

	const auto& buf = engine.results();

	// loop over shell quartets
	for (auto s1=0; s1 != shells.size(); ++s1) {
    
		auto bf1_first = shell2bf[s1];
		auto n1 = shells[s1].size();

		for (auto s2=0; s2 <= s1; ++s2) {
      
			auto bf2_first = shell2bf[s2];
			auto n2 = shells[s2].size();
      
			for (auto s3=0; s3 <= s1; ++s3) {
	
				auto bf3_first = shell2bf[s3];
				auto n3 = shells[s3].size();

				const auto s4_max = (s1 == s3) ? s2 : s3;
				for(auto s4=0; s4<=s4_max; ++s4) {
	  
					auto bf4_first = shell2bf[s4];
					auto n4 = shells[s4].size();

					engine.compute(shells[s1], shells[s2], shells[s3], shells[s4]);

					const auto* buf_1234 = buf[0];
					if (buf_1234 == nullptr)
						continue;

					for(auto f1=0, f1234=0; f1!=n1; ++f1) {
            
						const auto bf1 = f1 + bf1_first;
            
						for(auto f2=0; f2!=n2; ++f2) {
              
							const auto bf2 = f2 + bf2_first;
             
							for(auto f3=0; f3!=n3; ++f3) {
                
								const auto bf3 = f3 + bf3_first;
             
								for(auto f4=0; f4!=n4; ++f4, ++f1234) {
                  
									const auto bf4 = f4 + bf4_first;
									int mu = std::min(bf1, bf2);
									int nu = std::max(bf1, bf2);
									int sigma = std::min(bf3, bf4);
									int lambda = std::max(bf3, bf4); 
									int l1 = (nu*(nu+1))/2;
									int l3 = (lambda*(lambda+1))/2; 
									values[mu + l1 + sigma*l2 + l3*l2] = buf_1234[f1234]; 
									values[sigma + l3 + mu*l2 + l1*l2] = buf_1234[f1234]; 
								}
							}
						}
					}
				}
			}
		}
	}
	aoInts.write(sz, ix, values);
	delete ix; 
	delete values; 
	
	int64_t sz1, sz2, *ix1, *ix2;
	double *v1, *v2; 
	
	CTF::Matrix<> cp_occ(N, nocc, NS, dw); 
	CTF::Matrix<> cp_virt(N, nvirt, NS, dw); 
	Matrix& CP = focker.getCP(); 

	int ctr = 0; 
	cp_occ.read_local(&sz1, &ix1, &v1);
	for (int i = 0; i < nocc; i++)
		for (int mu = 0; mu < N; mu++) 
			v1[ctr++] = CP(mu, i); 
	cp_occ.write(sz1, ix1, v1); 
	
	delete ix1;
	delete v1; 
	
	ctr = 0; 
	cp_virt.read_local(&sz2, &ix2, &v2);
	for (int a = nocc; a < N; a++)
		for (int mu = 0; mu < N; mu++)
			v2[ctr++] = CP(mu, a); 
	cp_virt.write(sz2, ix2, v2); 
	
	delete ix2; 
	delete v2;
	
	CTF::Tensor<> Viajb(4, iajbshape, nsns, dw); 
	CTF::Tensor<> Vijab(4, ijabshape, sysy, dw);
	 
	{
		int t1_lens[4] = {N, N, N, nvirt}; 
		int t2_lens[4] = {N, N, nvirt, nvirt};
		int t3_lens[4] = {N, nocc, nvirt, nvirt};
		CTF::Tensor<> temp1(4, t1_lens, syns, dw);
		CTF::Tensor<> temp2(4, t2_lens, sysy, dw); 
		CTF::Tensor<> temp3(4, t3_lens, nssy, dw);	
		
		temp1["uvwb"] = cp_virt["xb"] * aoInts["uvwx"]; 
		temp2["uvab"] = cp_virt["wa"] * temp1["uvwb"]; 
		temp3["ujab"] = cp_occ["vj"] * temp2["uvab"];
		Vijab["ijab"] = 0.25 * cp_occ["ui"] * temp3["ujab"]; 
	} 
	
	
	{
		int t1_lens[4] = {N, N, N, nvirt};
		int t2_lens[4] = {N, N, nocc, nvirt};
		int t3_lens[4] = {N, nvirt, nocc, nvirt}; 
		CTF::Tensor<> temp1(4, t1_lens, syns, dw);
		CTF::Tensor<> temp2(4, t2_lens, syns, dw); 
		CTF::Tensor<> temp3(4, t3_lens, nsns, dw);	
		
		temp1["uvwb"] = cp_virt["xb"] * aoInts["uvwx"]; 
		temp2["uvjb"] = cp_occ["wj"] * temp1["uvwb"]; 
		temp3["uajb"] = cp_virt["va"] * temp2["uvjb"];
		Viajb["iajb"] = cp_occ["ui"] * temp3["uajb"]; 
	}  
	
	moInts.push_back(Viajb);
	moInts.push_back(Vijab);
}

void RPA::full_eris(std::vector<CTF::Tensor<> > &moInts) {
	
	int ao_lens[4] = {N, N, N, N};
	int iajbshape[4] = {nocc, nvirt, nocc, nvirt};
	int ijabshape[4] = {nocc, nocc, nvirt, nvirt}; 
	int sysy[4] = {SY, NS, SY, NS}; 
	int nssy[4] = {NS, NS, SY, NS};
	int nsns[4] = {NS, NS, NS, NS};
	int syns[4] = {SY, NS, NS, NS}; 
	CTF::Tensor<> ao_integrals(4, ao_lens, sysy, dw); 
	int64_t sz, *i1, *i2, *i3;
	double *v1, *v2, *v3;
	int ctr = 0;
			
	ao_integrals.read_local(&sz, &i1, &v1);
	
	IntegralEngine& intengine = focker.getIntegrals();
	if(focker.getMolecule()->control->get_option<bool>("direct")) {
		S8EvenTensor4 aoInts = intengine.compute_eris(focker.getMolecule()->getBasis().getIntShells());
		for (int a = 0; a < N; a++) 
			for (int b = 0; b <= a; b++)
				for (int c = 0; c < N; c++)
					for (int d = 0; d <= c; d++) 
						v1[ctr++] = aoInts(d, c, b, a); 
	} else {
		for (int a = 0; a < N; a++) 
			for (int b = 0; b <= a; b++)
				for (int c = 0; c < N; c++)
					for (int d = 0; d <= c; d++) 
						v1[ctr++] = intengine.getERI(d, c, b, a); 
	}

	ao_integrals.write(sz, i1, v1);
	
	delete v1; 
	delete i1; 
	
	CTF::Matrix<> cp_occ(N, nocc, NS, dw); 
	CTF::Matrix<> cp_virt(N, nvirt, NS, dw); 
	Matrix& CP = focker.getCP(); 
	
	ctr = 0; 
	cp_occ.read_local(&sz, &i2, &v2);
	for (int i = 0; i < nocc; i++)
		for (int mu = 0; mu < N; mu++)
			v2[ctr++] = CP(mu, i); 
	cp_occ.write(sz, i2, v2); 
	
	delete i2;
	delete v2; 
	
	ctr = 0; 
	cp_virt.read_local(&sz, &i3, &v3);
	for (int a = nocc; a < N; a++)
		for (int mu = 0; mu < N; mu++)
			v3[ctr++] = CP(mu, a); 
	cp_virt.write(sz, i3, v3); 
	
	delete i3; 
	delete v3; 
	
	CTF::Tensor<> Viajb(4, iajbshape, nsns, dw); 
	CTF::Tensor<> Vijab(4, ijabshape, sysy, dw);
	 
	{
		int t1_lens[4] = {N, N, N, nvirt}; 
		int t2_lens[4] = {N, N, nvirt, nvirt};
		int t3_lens[4] = {N, nocc, nvirt, nvirt};
		CTF::Tensor<> temp1(4, t1_lens, syns, dw);
		CTF::Tensor<> temp2(4, t2_lens, sysy, dw); 
		CTF::Tensor<> temp3(4, t3_lens, nssy, dw);	
		
		temp1["uvwb"] = cp_virt["xb"] * ao_integrals["uvwx"]; 
		temp2["uvab"] = cp_virt["wa"] * temp1["uvwb"]; 
		temp3["ujab"] = cp_occ["vj"] * temp2["uvab"];
		Vijab["ijab"] = 0.25 * cp_occ["ui"] * temp3["ujab"]; 
	} 
	
	{
		int t1_lens[4] = {N, N, N, nvirt};
		int t2_lens[4] = {N, N, nocc, nvirt};
		int t3_lens[4] = {N, nvirt, nocc, nvirt}; 
		CTF::Tensor<> temp1(4, t1_lens, syns, dw);
		CTF::Tensor<> temp2(4, t2_lens, syns, dw); 
		CTF::Tensor<> temp3(4, t3_lens, nsns, dw);	
		
		temp1["uvwb"] = cp_virt["xb"] * ao_integrals["uvwx"]; 
		temp2["uvjb"] = cp_occ["wj"] * temp1["uvwb"]; 
		temp3["uajb"] = cp_virt["va"] * temp2["uvjb"];
		Viajb["iajb"] = cp_occ["ui"] * temp3["uajb"]; 
	}  
	
	moInts.push_back(Viajb);
	moInts.push_back(Vijab);
	 
}

