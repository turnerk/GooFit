#include "ExpGausPdf.hh"

EXEC_TARGET fptype device_ExpGaus (fptype* evt, fptype* p, unsigned int* indices) {
  float x     = evt[indices[2 + indices[0]]]; 
  float mean  = p[indices[1]];
  float sigma = p[indices[2]];
  float alpha = p[indices[3]];

  float ret = 0.5*alpha; 
  float exparg = ret * (2*mean + alpha*sigma*sigma - 2*x);
  float erfarg = (mean + alpha*sigma*sigma - x) / (sigma * 1.4142135623);

  ret *= EXP(exparg); 
  ret *= ERFC(erfarg); 

  return ret; 
}

MEM_DEVICE device_function_ptr ptr_to_ExpGaus = device_ExpGaus; 

ExpGausPdf::ExpGausPdf (std::string n, Variable* _x, Variable* mean, Variable* sigma, Variable* tau) 
  : GooPdf(_x, n)
{
  std::vector<unsigned int> pindices;
  pindices.push_back(registerParameter(mean));
  pindices.push_back(registerParameter(sigma));
  pindices.push_back(registerParameter(tau));
  GET_FUNCTION_ADDR(ptr_to_ExpGaus);
  initialise(pindices); 
}


