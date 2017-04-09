#ifndef MP2HEADERDEF
#define MP2HEADERDEF

#include "tensor4.hpp"
#include "fock.hpp"
#include <ctf.hpp>

class IntegralEngine;

class Integrals {
  public:
  CTF::World * dw;
  CTF::Tensor<> * aa;
  CTF::Tensor<> * ii;
  CTF::Tensor<> * ab;
  CTF::Tensor<> * ai;
  CTF::Tensor<> * ia;
  CTF::Tensor<> * ij;
  CTF::Tensor<> * abcd;
  CTF::Tensor<> * abci;
  CTF::Tensor<> * aibc;
  CTF::Tensor<> * aibj;
  CTF::Tensor<> * abij;
  CTF::Tensor<> * ijab;
  CTF::Tensor<> * aijk;
  CTF::Tensor<> * ijak;
  CTF::Tensor<> * ijkl;

  Integrals(int no, int nv, CTF::World &dw_){
    int shapeASAS[] = {AS,NS,AS,NS};
    int shapeASNS[] = {AS,NS,NS,NS};
    int shapeNSNS[] = {NS,NS,NS,NS};
    int shapeNSAS[] = {NS,NS,AS,NS};
    int vvvv[]      = {nv,nv,nv,nv};
    int vvvo[]      = {nv,nv,nv,no};
    int vovv[]      = {nv,no,nv,nv};
    int vovo[]      = {nv,no,nv,no};
    int vvoo[]      = {nv,nv,no,no};
    int oovv[]      = {no,no,nv,nv};
    int vooo[]      = {nv,no,no,no};
    int oovo[]      = {no,no,nv,no};
    int oooo[]      = {no,no,no,no};
    
    dw = &dw_;
    
    aa = new CTF_Vector(nv,dw_);
    ii = new CTF_Vector(no,dw_);
    
    ab = new CTF_Matrix(nv,nv,AS,dw_,"Vab",1);
    ai = new CTF_Matrix(nv,no,NS,dw_,"Vai",1);
    ia = new CTF_Matrix(no,nv,NS,dw_,"Via",1);
    ij = new CTF_Matrix(no,no,AS,dw_,"Vij",1);

    abcd = new CTF::Tensor<>(4,vvvv,shapeASAS,dw_,"Vabcd",1);
    abci = new CTF::Tensor<>(4,vvvo,shapeASNS,dw_,"Vabci",1);
    aibc = new CTF::Tensor<>(4,vovv,shapeNSAS,dw_,"Vaibc",1);
    aibj = new CTF::Tensor<>(4,vovo,shapeNSNS,dw_,"Vaibj",1);
    abij = new CTF::Tensor<>(4,vvoo,shapeASAS,dw_,"Vabij",1);
    ijab = new CTF::Tensor<>(4,oovv,shapeASAS,dw_,"Vijab",1);
    aijk = new CTF::Tensor<>(4,vooo,shapeNSAS,dw_,"Vaijk",1);
    ijak = new CTF::Tensor<>(4,oovo,shapeASNS,dw_,"Vijak",1);
    ijkl = new CTF::Tensor<>(4,oooo,shapeASAS,dw_,"Vijkl",1);
  }

  ~Integrals(){
    delete aa;
    delete ii;
    
    delete ab;
    delete ai;
    delete ia;
    delete ij;
    
    delete abcd;
    delete abci;
    delete aibc;
    delete aibj;
    delete abij;
    delete ijab;
    delete aijk;
    delete ijak;
    delete ijkl;
  }
  
  CTF::Idx_Tensor operator[](char const * idx_map_){
    int i, lenm, no, nv;
    lenm = strlen(idx_map_);
    char new_idx_map[lenm+1];
    new_idx_map[lenm]='\0';
    no = 0;
    nv = 0;
    for (i=0; i<lenm; i++){
      if (idx_map_[i] >= 'a' && idx_map_[i] <= 'h'){
        new_idx_map[i] = 'a'+nv;
        nv++;
      } else if (idx_map_[i] >= 'i' && idx_map_[i] <= 'n'){
        new_idx_map[i] = 'i'+no;
        no++;
      }
    }
//    printf("indices %s are %s\n",idx_map_,new_idx_map);
    if (0 == strcmp("a",new_idx_map)) return (*aa)[idx_map_];
    if (0 == strcmp("i",new_idx_map)) return (*ii)[idx_map_];
    if (0 == strcmp("ab",new_idx_map)) return (*ab)[idx_map_];
    if (0 == strcmp("ai",new_idx_map)) return (*ai)[idx_map_];
    if (0 == strcmp("ia",new_idx_map)) return (*ia)[idx_map_];
    if (0 == strcmp("ij",new_idx_map)) return (*ij)[idx_map_];
    if (0 == strcmp("abcd",new_idx_map)) return (*abcd)[idx_map_];
    if (0 == strcmp("abci",new_idx_map)) return (*abci)[idx_map_];
    if (0 == strcmp("aibc",new_idx_map)) return (*aibc)[idx_map_];
    if (0 == strcmp("aibj",new_idx_map)) return (*aibj)[idx_map_];
    if (0 == strcmp("abij",new_idx_map)) return (*abij)[idx_map_];
    if (0 == strcmp("ijab",new_idx_map)) return (*ijab)[idx_map_];
    if (0 == strcmp("aijk",new_idx_map)) return (*aijk)[idx_map_];
    if (0 == strcmp("ijak",new_idx_map)) return (*ijak)[idx_map_];
    if (0 == strcmp("ijkl",new_idx_map)) return (*ijkl)[idx_map_];
    printf("Invalid integral indices\n");
    assert(0);
//shut up compiler
    return (*aa)[idx_map_];
  }
};

class Amplitudes {
  public:
  CTF::Tensor<> * ai;
  CTF::Tensor<> * abij;
  CTF::World * dw;

  Amplitudes(int no, int nv, CTF::World &dw_){
    dw = &dw_;
    int shapeASAS[] = {AS,NS,AS,NS};
    int vvoo[]      = {nv,nv,no,no};

    ai = new CTF_Matrix(nv,no,NS,dw_,"Tai",1);
    abij = new CTF::Tensor<>(4,vvoo,shapeASAS,dw_,"Tabij",1);
  }

  ~Amplitudes(){
    delete ai;
    delete abij;
  }

  CTF::Idx_Tensor operator[](char const * idx_map_){
    if (strlen(idx_map_) == 4) return (*abij)[idx_map_];
    else return (*ai)[idx_map_];
  }
};

class MP2
{
private:
	int N, nocc;
	double energy;
	S8EvenTensor4 moInts;
	std::shared_ptr<Integrals> spinInts;
	std::shared_ptr<Amplitudes> amplitudes;  
	bool spinBasis;
	Fock& focker;
public:
	CTF::World dw; 
	
	MP2(Fock& _focker);
	void transformIntegrals();
	void transformThread(int start, int end, Tensor4& moTemp);
	void calculateEnergy();
	double getEnergy() const { return energy; }
	std::shared_ptr<Integrals>& getSpinInts() { 
		return spinInts;
	}
	std::shared_ptr<Amplitudes>& getAmplitudes() {
		return amplitudes;
	}
	S8EvenTensor4& getMOInts() {
		return moInts;
	}
	int getN() const { return N; }
	int getNocc() const { return nocc; }
	Fock& getFock() { return focker; }
};

#endif 
