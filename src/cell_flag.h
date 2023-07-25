#pragma once

#include "cysift.h"

#include <bitset>
#include <stdexcept>
#include <string>
#include <iostream>
#include <sstream>
#include <vector>
#include <cmath>

class CellFlag {
public:
  
  CellFlag() { bitmap = 0; }
  
  explicit CellFlag(cy_uint num) {
    bitmap = num;
  }
  
  std::string toBitString() const;
  
  void setFlagOn(int n); 
  
  void setFlagOff(int n);
  
  std::string toString() const;
  
  friend std::ostream& operator<<(std::ostream& os, const CellFlag& cellFlag);
  
  cy_uint toBase10() const;
  
  bool testAndOr(cy_uint logor, cy_uint logand) const;
  
private:
  
  // data
  cy_uint bitmap = 0;
  
  void check_bounds(int n) const;
};

	    
