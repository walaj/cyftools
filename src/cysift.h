#pragma once
#include <unordered_map>
#include <string>
#include <vector>
#include <sstream>
#include <iostream>

#define TUMOR_FLAG 0b1
#define MARK_FLAG 0b10
#define MARGIN_FLAG 0b100

// Macro to test if a flag is set
#define IS_FLAG_SET(flags, flag) ((flags) & (flag))

// Macro to set a flag
#define SET_FLAG(flags, flag) ((flags) |= (flag))

// Macro to clear a flag
#define CLEAR_FLAG(flags, flag) ((flags) &= ~(flag))

// define phenotype map
typedef std::pair<float,float> Pheno;
typedef std::unordered_map<std::string, Pheno> PhenoMap;

#ifdef USE_64_BIT
typedef uint64_t cy_uint;
#else
typedef uint32_t cy_uint;
#endif

// object to parse the radial density scripts
struct RadialSelector {

  explicit RadialSelector(const std::string& line) {

    std::vector<std::string> tokens = split(line, ',');
    
    if (tokens.size() != 7) {
      throw std::runtime_error("There must be exactly 7 tokens: " + line);
    }

    int_data.resize(6);
    for (int i = 0; i < 6; i++) {
      try {
	int_data[i] = std::stoi(tokens[i]);
      } catch (const std::invalid_argument &e) {
	throw std::runtime_error("The first 6 tokens must be integers.");
      }
    }

    label = tokens[6];
    
  }

  // split tokens
  std::vector<std::string> split(const std::string &input, char delimiter) {
    std::vector<std::string> tokens;
    std::istringstream stream(input);
    std::string token;
    
    while (std::getline(stream, token, delimiter)) {
      tokens.push_back(token);
    }
    
    return tokens;
  }

  friend std::ostream& operator<<(std::ostream& os, const RadialSelector& rs) {
    for (size_t i = 0; i < rs.int_data.size(); ++i) {
      os << rs.int_data[i];
      if (i < rs.int_data.size() - 1) {
	os << ",";
      }
    }
    os << "," << rs.label;
    return os;
  }
  
  std::vector<cy_uint> int_data;
  std::string label;
  
};

