#include "cell_flag.h"

std::ostream& operator<<(std::ostream& os, const CellFlag& cellFlag) {
  os << cellFlag.toString(); 
  return os;
}

bool CellFlag::testAndOr(cy_uint logor, cy_uint logand) const {

  //bool resultN;
  //bool result;
  //std::bitset<64> bs(bitmap);
  
  // new way
  bool orCondition = logor == 0 || (bitmap & logor) != 0; // True if logor not specified or any OR flags are turned on

  if (orCondition) {
    bool andCondition = (bitmap & logand) == logand; // True if all AND flags are turned on
    return andCondition || logand == 0; // return true if all AND flags are turned on or if logand not specified
  }

  return false;
  /*
  // old way
  std::bitset<BITMAP_SIZE> orBitset(logor);
  std::bitset<BITMAP_SIZE> andBitset(logand);

  bool orCondition = (bs & orBitset).any(); // True if any OR flags are turned on
  bool andCondition = (bs & andBitset) == andBitset; // True if all AND flags are turned on

  if (logor == 0) // don't gate on OR if not specified
    orCondition = true;

  if (orCondition && andCondition) {
    return true;
  } else if (orCondition && !andBitset.any()) {
    return true;
  } else {
    return false;
  }

  //std::cerr << " OLD " << result << " NEW " << resultN << " flag " << bitmap << " logor " << logor << " logand " << logand << std::endl;
  return false;
  */
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
