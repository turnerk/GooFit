#include "DalitzVetoPdf.hh"
#include "DalitzPlotHelpers.hh" 

EXEC_TARGET fptype device_DalitzVeto (fptype* evt, fptype* p, unsigned int* indices) {
  int idx[6];
  idx[0] = indices[0];
  idx[1] = indices[1];
  idx[2] = indices[2];
  idx[3] = indices[3];
  idx[4] = indices[4];
  idx[5] = indices[5];

  fptype x         = evt[indices[2 + idx[0] + 0]]; 
  fptype y         = evt[indices[2 + idx[0] + 1]]; 

  fptype motherM   = p[idx[1]];
  fptype d1m       = p[idx[2]];
  fptype d2m       = p[idx[3]];
  fptype d3m       = p[idx[4]];

  fptype motherM2  = motherM*motherM;
  fptype d1m2      = d1m*d1m;
  fptype d2m2      = d2m*d2m;
  fptype d3m2      = d3m*d3m;

  fptype massSum   = motherM2 + d1m2 + d2m2 + d3m2;

  fptype ret = inDalitz(x, y, motherM, d1m, d2m, d3m) ? 1.0 : 0.0; 
  unsigned int numVetos = idx[5];

  fptype z         = massSum - x - y;

  //unroll?
#pragma unroll
  for (int i = 0; i < numVetos; ++i) {
    int i3 = i*3;

    unsigned int varIndex =   indices[6 + i3 + 0];
    fptype minimum        = p[indices[6 + i3 + 1]];
    fptype maximum        = p[indices[6 + i3 + 2]];
    fptype currDalitzVar = (PAIR_12 == varIndex ? x : PAIR_13 == varIndex ? y : z);

    ret *= ((currDalitzVar < maximum) && (currDalitzVar > minimum)) ? 0.0 : 1.0;
  }

  return ret; 
}

MEM_DEVICE device_function_ptr ptr_to_DalitzVeto = device_DalitzVeto;

__host__ DalitzVetoPdf::DalitzVetoPdf (std::string n, Variable* _x, Variable* _y, Variable* motherM, Variable* d1m, Variable* d2m, Variable* d3m, vector<VetoInfo*> vetos) 
  : GooPdf(0, n) 
{
  registerObservable(_x);
  registerObservable(_y);

  std::vector<unsigned int> pindices;
  pindices.push_back(registerParameter(motherM));
  pindices.push_back(registerParameter(d1m));
  pindices.push_back(registerParameter(d2m));
  pindices.push_back(registerParameter(d3m));

  pindices.push_back(vetos.size()); 
  for (vector<VetoInfo*>::iterator v = vetos.begin(); v != vetos.end(); ++v) {
    pindices.push_back((*v)->cyclic_index);
    pindices.push_back(registerParameter((*v)->minimum));
    pindices.push_back(registerParameter((*v)->maximum));
  }

  GET_FUNCTION_ADDR(ptr_to_DalitzVeto);
  initialise(pindices); 
}
