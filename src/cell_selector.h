#pragma once

#include "cysift.h"

// forward declar
class CellFlag;

class CellSelector {

public:

 CellSelector(cy_uint po, cy_uint pa, cy_uint pn,
	      cy_uint co, cy_uint ca, cy_uint cn) : 
  plogor(po), plogand(pa), pnot(pn),
    clogor(co), clogand(ca), cnot(cn) {}

 CellSelector() : plogor(0), plogand(0), pnot(0),
    clogor(0), clogand(0), cnot(0) {}
  
  bool TestFlags(cy_uint pflag, cy_uint cflag) const;
  
  bool TestFlags(CellFlag pflag, CellFlag cflag) const;
  
private:
  
  cy_uint plogor = 0;
  cy_uint plogand = 0;
  bool pnot = 0;

  cy_uint clogor = 0;
  cy_uint clogand = 0;
  bool cnot = 0;
  
};
