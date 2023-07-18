#pragma once

#include <vector>
#include <string>
#include <sstream>
#include <iomanip>
#include <numeric>
#include <iostream>

#include "cell_column.h"

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

std::vector<std::string> tokenize_comma_delimited(const std::string& str);

float euclidean_distance(float x1, float y1, float x2, float y2);

bool check_readable(const std::string& filename);

float euclidean_distance_squared(float x1, float y1, float x2, float y2);

#ifdef HAVE_HDF5
void write_hdf5_dataframe_attributes(H5::Group& group);
#endif

enum class ColumnType;
std::string columnTypeToString(ColumnType type);
