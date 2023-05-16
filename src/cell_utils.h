#pragma once

#include <cmath>
#include <vector>
#include <string>
#include <sstream>
#include <iomanip>
#include <unordered_map>
#include <numeric>
#include <iostream>
#include <random>
#include <set>
#include <variant>
#include <unordered_set>

#include "cell_header.h"
#include "cell_column.h"

/// @brief Alias for a single datum from a cell, which can be int, float or string
using CellDatum = std::variant<uint64_t, float, std::string>;

/// @brief Alias for a row of cells, which can contain integers, floats, or strings.
using CellRow = std::vector<std::variant<uint64_t, float, std::string>>;

/// @brief Alias for a function pointer that handles a CellRow
using CellRowFunc = CellRow(*)(const CellRow&);

/** Format an integer to include commas
 * @param data Number to format
 * @return String with formatted number containing commas
 */
template <typename T> inline std::string AddCommas(T data) {
  std::stringstream ss; 
  ss << std::fixed << data; 
  std::string s = ss.str();
  if (s.length() > 3)
    for (int i = s.length()-3; i > 0; i -= 3)
      s.insert(i,",");
  return s;
}

void process_token_to_variant(const std::string_view& token,
			      const Tag& tag,
			      CellDatum& value);

void get_two_elements_as_floats(const std::string_view& str, size_t n, size_t m,
				float &x, float &y);

std::string round_string(const std::string& str, int precision);

std::string tokens_to_comma_string(const std::vector<std::string>& input);

std::vector<std::string> tokenize_comma_delimited(const std::string& str);

int get_nth_element_as_integer(const std::string_view& str, size_t n);

std::string exclude_elements(const std::string_view& str, const std::unordered_set<size_t>& exclude_set);

int read_one_line_to_cellrow(const std::string& line,
			     CellRow& values,
			     const CellHeader& m_header);

/*std::set<int> sampleWithoutReplacement(int n, int N, int seed) {

  assert(n <= N);
  std::set<int> s;
  std::vector<int> v(N);
  for (int i = 0; i < N; i++) {
    v[i] = i; // initialize the vector with 0-based indices
  }
  
  // use a random device to generate random numbers
  std::random_device rd;
  std::mt19937 gen(seed);
  std::uniform_int_distribution<> dis(0, N-1);
  
  // randomly sample n indices without replacement
    while (s.size() < n) {
      int index = dis(gen); // generate a random index
      if (s.count(index) == 0) { // check if the index has already been sampled
	s.insert(index); // if not, insert it into the set
      }
    }
    
    // convert the sampled indices back to 0-based indexing and return the set
    std::set<int> sampledIndices;
    for (int index : s) {
      sampledIndices.insert(v[index]);
    }
    return sampledIndices;
}
*/

enum class ColumnType;
std::string columnTypeToString(ColumnType type);

std::string variantToString(const std::variant<uint64_t, float, std::string>& value);
