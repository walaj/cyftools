#include "cell_selector.h"
#include "cell_flag.h"

#include <iomanip> // For std::setw and std::setfill

size_t CellSelector::size() const {
  return m_ors.size();
}

bool SelectionUnit::isEmpty() const {
  return (pand + pnot + cand + cnot) == 0;
}

std::ostream& operator<<(std::ostream& os, const SelectionUnit& unit) {
  os << "P-AND: " << std::left << std::setw(15) << std::setfill(' ') << unit.pand << ", ";
  os << "P-NOT: " << std::left << std::setw(15) << std::setfill(' ') << unit.pnot << ", ";
  os << "C-AND: " << std::left << std::setw(15) << std::setfill(' ') << unit.cand << ", ";
  os << "C-NOT: " << std::left << std::setw(15) << std::setfill(' ') << unit.cnot << ", ";
  os << "P-OR: "  << std::left << std::setw(15) << std::setfill(' ') << unit.por << ", ";
  os << "C-OR: "  << std::left << std::setw(15) << std::setfill(' ') << unit.cor;
  return os;
}

std::ostream& operator<<(std::ostream& os, const CellSelector& cs) {

  if (!cs.m_ors.size()) {
    //os << "NO CONDITIONS SPECIFIED" << std::endl;
    return os;
  }

  os << "---- OR CONDITIONS BELOW ----" << std::endl;
  for (size_t i = 0; i < cs.m_ors.size(); i++) {
    os << cs.m_ors.at(i) << std::endl;
  }
  return os;
}

bool SelectionUnit::TestFlags(cy_uint pflag, cy_uint cflag) const {

  // bool pass = true;
  // std::cerr << "por " << por << " cor " << cor << " pand " << pand << " cand " << cand << " flag " << pflag << " cflag " << cflag << std::endl;
  // std::cerr << "pflag | por " << (pflag | por) << std::endl;
  
  // check the or conditions. If an OR, turn off if doesn't pass
  // this should pass if ANY of the flag bits are on (check if on with &)
  // remember that you check if bits are on with &, and if you want a logical AND
  // then the bits that are ON should leave you with same full number as flag combo to test
  // otherwise, if you want to OR the bits that are on, just check that ANY of them are on (> 0)
  if (por > 0 && (pflag & por) == 0)
    return false;
  if (cor > 0 && (cflag & cor) == 0)
    return false;
  
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
