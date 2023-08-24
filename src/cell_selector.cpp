#include "cell_selector.h"
#include "cell_flag.h"

bool CellSelector::TestFlags(cy_uint pflag, cy_uint cflag) const {
  return TestFlags(CellFlag(pflag), CellFlag(cflag));
}

bool CellSelector::TestFlags(CellFlag pflag, CellFlag cflag) const {
  
  if (plogor == 0 &&
      plogand == 0 &&
      clogand == 0 &&
      clogor == 0)
    return true;
      
  bool pflags_met = pflag.testAndOr(plogor, plogand);
  bool cflags_met = cflag.testAndOr(clogor, clogand);
  if ( (pflags_met != pnot) && (cflags_met != cnot))
    return true;
  return false;
}

