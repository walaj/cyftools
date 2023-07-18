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
#include <unordered_set>
#include <string_view>

#include "cysift.h"
#include "cell_header.h"
#include "cell_column.h"
#include "cell_row.h"

#ifdef HAVE_HDF5
#include <H5Cpp.h>
#endif

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

void column_to_row_major(std::vector<float>& data, int nobs, int ndim);

PhenoMap phenoread(const std::string& filename);

void get_two_elements_as_floats(const std::string_view& str, size_t n, size_t m,
				float &x, float &y);

std::string round_string(const std::string& str, int precision);

std::string tokens_to_comma_string(const std::vector<std::string>& input);

std::vector<std::string> tokenize_comma_delimited(const std::string& str);

int get_nth_element_as_integer(const std::string_view& str, size_t n);

std::string exclude_elements(const std::string_view& str, const std::unordered_set<size_t>& exclude_set);

float euclidean_distance(float x1, float y1, float x2, float y2);

bool check_readable(const std::string& filename);

float euclidean_distance_squared(float x1, float y1, float x2, float y2);

float manhattan_distance(float x1, float y1, float x2, float y2);

#ifdef HAVE_HDF5
void write_hdf5_dataframe_attributes(H5::Group& group);
#endif
enum class ColumnType;
std::string columnTypeToString(ColumnType type);

