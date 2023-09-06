#pragma once

#include "cysift.h"

// forward declar
class CellFlag;

struct SelectionUnit {

  cy_uint pand = 0;
  cy_uint pnot = 0;
  cy_uint cand = 0;
  cy_uint cnot = 0;

  bool TestFlags(cy_uint pflag, cy_uint cflag) const;

  friend std::ostream& operator<<(std::ostream& os, const SelectionUnit& unit);
};

class CellSelector {

public:

  CellSelector() {}

  void AddSelectionUnit(SelectionUnit u) { m_ors.push_back(u); }
  
  bool TestFlags(cy_uint pflag, cy_uint cflag) const;
  
  friend std::ostream& operator<<(std::ostream& os, const CellSelector& cs);
  
private:

  std::vector<SelectionUnit> m_ors;
  
};
