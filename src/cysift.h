#pragma once
#include <unordered_map>
#include <string>
#include <vector>
#include <sstream>
#include <iostream>
#include <cassert>

#define TUMOR_FLAG 1
#define MARK_FLAG 2
#define MARGIN_FLAG 4
#define TUMOR_MANUAL_FLAG 8
#define TCELL_FLAG 16
#define TLS_FLAG 32
#define MARGIN_MANUAL_FLAG 64
#define NORMAL_FLAG 128
#define NORMAL_MARGIN_FLAG 256
#define BUILD_GRAPH_FLAG 2097152 

#define GLEASON_GRADE_GROUP_1 4096
#define GLEASON_GRADE_GROUP_2 8192
#define GLEASON_GRADE_GROUP_3 16384
#define GLEASON_GRADE_GROUP_4 32768
#define GLEASON_GRADE_GROUP_5 65536
#define PERINEURAL_INVASION 131072
#define SEMINAL_VESICLES 262144

#define PROSTATE_AMCAR 4
#define PROSTATE_HMWCK 8
#define PROSTATE_CD20 64
#define PROSTATE_CD4 1024
#define PROSTATE_CD3 2048
#define PROSTATE_PD1 32768
#define PROSTATE_CD8 4096
#define PROSTATE_FOXP3 16384
#define PROSTATE_CD57 65536

#define ORION_CD4 64
#define ORION_FOXP3 128
#define ORION_CD8 256
#define ORION_CD20 1024
#define ORION_PDL1 2048
#define ORION_CD3 4096
#define ORION_CD163 8192
#define ORION_PD1 32768
#define ORION_PANCK 131072

#define ORIONJHU_CD3 2048
#define ORIONJHU_CD20 64
#define ORIONJHU_CD8 256
#define ORIONJHU_CD163 32
#define ORIONJHU_FOXP3 16385
#define ORIONJHU_PANCK 65536
#define ORIONJHU_PD1 8192

// Macro to test if a flag is set
// Data flag (cflag or pflag) is one, the flag to query is two
// eg IS_FLAG_SET(cflag, TUMOR_FLAG)
#define IS_FLAG_SET(flags, flag) (((flags) & (flag)) == (flag))

// query flag by bit position (0 is OK)
#define IS_FLAG_I_SET(flags, i) ((flags & (1ull << (i))) != 0)

// Macro to test if any of the flags are set
#define IS_FLAG_SET_OR(flags, flag) ((flags) & (flag))

// are all flags off
#define ARE_FLAGS_OFF(flags, flag) (((flags) & (flag)) == 0)

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
	int_data[i] = std::stoll(tokens[i]);
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

