#ifndef CELL_FLAG_H
#define CELL_FLAG_H

#include <bitset>
#include <stdexcept>
#include <string>
#include <iostream>
#include <sstream>
#include <vector>
#include <cmath>

class CellFlag {
public:
  static const int BITMAP_SIZE = 64;
  static const std::string BASE64_CHARS;
  
  CellFlag() {}
  
  explicit CellFlag(const std::string& base64_string) {
      //fromBase10(base64_string);
    }
  
  explicit CellFlag(uint64_t num) {
      fromBase10(num);
      //fromBase10(base64_string);
    }
  
  std::string toBitString() const;
  
  void setFlagOn(int n); 
  
  void setFlagOff(int n);
  
  std::string toString() const;
  
  friend std::ostream& operator<<(std::ostream& os, const CellFlag& cellFlag);
  
  uint64_t toBase10_uint64_t() const;
  
  long long base64_to_base10(const std::string& base64_num);
  
  bool testAndOr(uint64_t logor, uint64_t logand) const;
  
  bool test(uint64_t on, uint64_t off) const;
    
  std::string toBase10() const;
  
  void fromBase10(uint64_t base10);
  
private:
  
  // data
  std::bitset<BITMAP_SIZE> bitmap;
  
  void check_bounds(int n) const;
  
  std::string to_base64() const;
  
  void from_base64(const std::string& base64_string);
};

#endif

	    
