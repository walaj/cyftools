#include "cell_selector.h"
#include "cell_flag.h"

std::ostream& operator<<(std::ostream& os, const SelectionUnit& unit) {
    os << "P-AND: " << unit.pand << ", ";
    os << "P-NOT: " << unit.pnot << ", ";
    os << "C-AND: " << unit.cand << ", ";
    os << "C-NOT: " << unit.cnot;
    return os;
}

std::ostream& operator<<(std::ostream& os, const CellSelector& cs) {

  if (!cs.m_ors.size()) {
    os << "NO CONDITIONS SPECIFIED" << std::endl;
    return os;
  }
  
  for (size_t i = 0; i < cs.m_ors.size(); i++) {
    os << cs.m_ors.at(i) << std::endl;
    if (i != cs.m_ors.size() - 1)
      os << "---- OR ----" << std::endl;
  }
  return os;
}

bool SelectionUnit::TestFlags(cy_uint pflag, cy_uint cflag) const {

  // Check the AND conditions
  if ((pflag & pand) != pand) return false; // All bits in pand should be set in pflag
  if ((cflag & cand) != cand) return false; // All bits in cand should be set in cflag
  
  // Check the NOT conditions
  if (pflag & pnot) return false; // Any bit in pnot should NOT be set in pflag
  if (cflag & cnot) return false; // Any bit in cnot should NOT be set in cflag
  
  return true; // If passed all checks
  
}

bool CellSelector::TestFlags(cy_uint pflag, cy_uint cflag) const {

  for (const auto& s : m_ors)
    if (s.TestFlags(pflag, cflag))
      return true;

  return false;
}

/*bool CellSelector::TestFlags(CellFlag pflag, CellFlag cflag) const {
  
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
*/
