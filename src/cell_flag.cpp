#include "cell_flag.h"

std::ostream& operator<<(std::ostream& os, const CellFlag& cellFlag) {
  os << cellFlag.toString(); 
  return os;
}

bool CellFlag::testAndOr(cy_uint logor, cy_uint logand) const {

  // no conditions specified, so return true
  if (static_cast<cy_uint>(-1) == logor &&
      static_cast<cy_uint>(-1) == logand)
    return true;
  
  // First, test the OR condition
  bool orCondition = logor == 0 ||
    static_cast<cy_uint>(-1) || (bitmap & logor) != 0; 

  // If OR is true and nothing to test for and, return true
  /*  bool andNotSpecified = logand == static_cast<cy_uint>(-1) ||
    (logand == 0)
  if (orCondition && ( || logand == 0)
  */
  // if the OR condition is true, then test the AND condition
  if (orCondition) {
    bool andCondition = (bitmap & logand) == logand; // True if all AND flags are turned on
    return andCondition || logand == 0; // return true if all AND flags are turned on or if logand not specified
  }

  return false;
}

std::string CellFlag::toBitString() const {

  std::bitset<64> bs(bitmap);
  //return std::to_string(bitmap);
  std::string bitmap_str = bs.to_string();
  return bitmap_str;
  
}

cy_uint CellFlag::toBase10() const {
  return bitmap; //bitmap.to_ullong();
}

void CellFlag::setFlagOn(int n) {
  check_bounds(n);
  bitmap |= (1ULL << n);
  //bitmap.set(n, true);
}

void CellFlag::setFlagOff(int n) {
  check_bounds(n);
  bitmap &= ~(1ULL << n);
  //Xbitmap.set(n, false);
}

std::string CellFlag::toString() const {
  return std::to_string(bitmap);
  //  return toBase10();
}

void CellFlag::check_bounds(int n) const {
  
#ifdef USE_64_BIT
  const int BITMAP_SIZE = 64;
#else
  const int BITMAP_SIZE = 32;
#endif
  
  if (n < 0 || n >= BITMAP_SIZE) {
    throw std::out_of_range("Flag index out of bounds");
  }
}
