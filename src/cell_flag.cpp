#include "cell_flag.h"

const std::string CellFlag::BASE64_CHARS = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";

std::ostream& operator<<(std::ostream& os, const CellFlag& cellFlag) {
  os << cellFlag.toString(); 
  return os;
}

std::string CellFlag::toBitString() const {
  
  std::string bitmap_str = bitmap.to_string();
  return bitmap_str;
  
}

void CellFlag::setFlagOn(int n) {
  check_bounds(n);
  bitmap.set(n, true);
}

void CellFlag::setFlagOff(int n) {
  check_bounds(n);
  bitmap.set(n, false);
}

std::string CellFlag::toString() const {
  return toBase10();
}

void CellFlag::check_bounds(int n) const {
  if (n < 0 || n >= BITMAP_SIZE) {
    throw std::out_of_range("Flag index out of bounds");
  }
}

std::string CellFlag::to_base64() const {
  int num_bytes = static_cast<int>(std::ceil(BITMAP_SIZE / 8.0));
  std::vector<unsigned char> bytes(num_bytes, 0);
  
  for (int i = 0; i < BITMAP_SIZE; ++i) {
    if (bitmap[i]) {
      bytes[i / 8] |= (1 << (7 - i % 8));
    }
  }
  
  std::stringstream ss;
  for (int i = 0; i < num_bytes; i += 3) {
    int a = bytes[i];
    int b = i + 1 < num_bytes ? bytes[i + 1] : 0;
    int c = i + 2 < num_bytes ? bytes[i + 2] : 0;
    
    int triple = (a << 16) | (b << 8) | c;
    ss << BASE64_CHARS[(triple >> 18) & 63]
       << BASE64_CHARS[(triple >> 12) & 63]
       << BASE64_CHARS[(triple >> 6) & 63]
       << BASE64_CHARS[triple & 63];
  }
  
  std::string base64_string = ss.str();
  int padding = (4 - base64_string.size() % 4) % 4;
  base64_string.replace(base64_string.end() - padding, base64_string.end(), padding, '=');
  
  return base64_string;
}



void CellFlag::from_base64(const std::string& base64_string) {

  fromBase10(10);
  
  //bitmap.set(4, true);
  std::string bitmap_str = bitmap.to_string();
    std::cerr << " bitmap_str " << bitmap_str << std::endl;
  std::cerr <<  " -- " << std::string(128 - bitmap_str.size(), '0') << std::endl;
  std::cerr << " 10(" << toBase10()<< ")" << std::endl;
  
  return;
  
  std::cerr << " BASE64 " << base64_string << std::endl;
  std::string decoded_base64 = base64_string;
  std::replace(decoded_base64.begin(), decoded_base64.end(), '=', 'A');

  std::cerr << " DECODED BASE64 " << decoded_base64 << std::endl;
  int num_bytes = static_cast<int>(std::ceil(BITMAP_SIZE / 8.0));
  std::vector<unsigned char> bytes(num_bytes, 0);

  std::cerr << " NUM BYES " << num_bytes << std::endl;
  
  for (size_t i = 0, j = 0; i < decoded_base64.size(); i += 4, j += 3) {
    int a = BASE64_CHARS.find(decoded_base64[i]);
    int b = BASE64_CHARS.find(decoded_base64[i + 1]);
    int c = BASE64_CHARS.find(decoded_base64[i + 2]);
    int d = BASE64_CHARS.find(decoded_base64[i + 3]);
    
    int triple = (a << 18) | (b << 12) | (c << 6) | d;
    bytes[j] = (triple >> 16) & 0xFF;
    if (j + 1 < num_bytes) bytes[j + 1] = (triple >> 8) & 0xFF;
    if (j + 2 < num_bytes) bytes[j + 2] = triple & 0xFF;
  }

  std::cerr << "\t" << base64_to_base10(decoded_base64) <<std::endl;

  bitmap.to_string();
  std::cerr << " bitmap_str " << bitmap_str << std::endl;
  std::cerr <<  " -- " << std::string(128 - bitmap_str.size(), '0') << std::endl;
  
  for (int i = 0; i < BITMAP_SIZE; ++i) {
    bitmap.set(i, (bytes[i / 8] & (1 << (7 - i % 8))) != 0);
  }

  bitmap_str = bitmap.to_string();
  std::cerr << " AFTER bitmap_str " << bitmap_str << std::endl;
  std::cerr <<  " -- " << std::string(128 - bitmap_str.size(), '0') << std::endl;
  
}


long long CellFlag::base64_to_base10(const std::string& base64_num) {
  const std::string base64_chars = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
  const int base = 64;

  // Calculate the length of the base-64 number
  const int len = base64_num.length();

  // Initialize the result to 0
  long long result = 0;

  // Iterate over each character in the base-64 number
  for (int i = 0; i < len; ++i) {
    // Find the position of the current character in the base-64 character set
    const int pos = base64_chars.find(base64_num[i]);

    // Multiply the current result by the base and add the value of the current character
    result = result * base + pos;
  }

  // Return the final result
  return result;
}


std::string CellFlag::toBase10() const {

  std::string result = "0";

  for (int i = bitmap.size() - 1; i >= 0; --i) {

    //  for (size_t i = 0; i < bitmap.size(); ++i) {
        // Shift left (multiply by 2)
        bool carry = false;
        for (auto it = result.rbegin(); it != result.rend(); ++it) {
            int value = (*it - '0') * 2 + (carry ? 1 : 0);
            carry = value >= 10;
            *it = '0' + (value % 10);
        }
        if (carry) {
            result.insert(result.begin(), '1');
        }

        // If the current bit is set, add 1 to the result
        if (bitmap[i]) {
            carry = true;
            for (auto it = result.rbegin(); it != result.rend() && carry; ++it) {
                int value = *it - '0' + 1;
                carry = value >= 10;
                *it = '0' + (value % 10);
            }
            if (carry) {
                result.insert(result.begin(), '1');
            }
        }
   }
    
    return result;
}

void CellFlag::fromBase10(uint64_t base10) {
  
  bitmap = std::bitset<BITMAP_SIZE>(base10);
  
}

bool CellFlag::test(uint64_t on, uint64_t off) const {

  std::bitset<BITMAP_SIZE> onBitset(on);
  std::bitset<BITMAP_SIZE> offBitset(off);
  
  for (size_t i = 0; i < bitmap.size(); ++i) {
    if (onBitset[i] && !bitmap[i]) {
      return false;
    }
    if (offBitset[i] && bitmap[i]) {
      return false;
    }
  }
  
  return true;
  
}


