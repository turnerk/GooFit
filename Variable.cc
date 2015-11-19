#include "Variable.hh"
#include <cmath> 

Variable::Variable (std::string n) 
  : Indexable(n) 
  , numbins(100)
  , fixed(false)
  , blind(0)
{
} 

Variable::Variable (std::string n, fptype v) 
  : Indexable(n, v)
  , error(0.002) 
  , lowerlimit(v - 0.01)
  , upperlimit(v + 0.01)
  , numbins(100)
  , fixed(true)
  , blind(0)
{
}

Variable::Variable (std::string n, fptype dn, fptype up) 
  : Indexable(n)
  , upperlimit(up)
  , lowerlimit(dn)
  , numbins(100)
  , fixed(false)
  , blind(0)
{
}

Variable::Variable (std::string n, fptype v, fptype dn, fptype up) 
  : Indexable(n, v)
  , error(0.1*(up-dn))
  , upperlimit(up)
  , lowerlimit(dn)
  , numbins(100)
  , fixed(false)
  , blind(0)
{
}

Variable::Variable (std::string n, fptype v, fptype e, fptype dn, fptype up) 
  : Indexable(n, v)
  , error(e)
  , upperlimit(up)
  , lowerlimit(dn)
  , numbins(100)
  , fixed(false)
  , blind(0)
{
}

CountVariable::CountVariable (std::string n) : Variable(n)
{
}

CountVariable::CountVariable (std::string n, fptype v) : Variable (n, v)
{
}

CountVariable::CountVariable (std::string n, fptype dn, fptype up) : Variable (n, dn, up)
{
}

CountVariable::CountVariable (std::string n, fptype v, fptype dn, fptype up) : Variable (n, v, dn, up)
{
}

CountVariable::CountVariable (std::string n, fptype v, fptype e, fptype dn, fptype up) : Variable (n, v, dn, up)
{
}

Variable::~Variable () {
}

