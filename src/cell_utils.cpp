#include "cell_utils.h"

#include <string_view>
#include <stdexcept>
#include <limits>
#include <cmath>
#include <fstream>

void column_to_row_major(std::vector<float>& data, int nobs, int ndim) {

  float* temp = new float[data.size()];
  
  for (int row = 0; row < nobs; ++row) {
    for (int col = 0; col < ndim; ++col) {
      int old_index = col * nobs + row;
      int new_index = row * ndim + col;
      temp[new_index] = data[old_index];
    }
  }
  
  for (size_t i = 0; i < data.size(); ++i) {
    data[i] = temp[i];
  }
  
  delete[] temp;
  
}

PhenoMap phenoread(const std::string& filename) {
  
  PhenoMap data;
  
  // open the phenotype file
  // of format: marker(string), low_bound(float), upper_bound(float)
  std::ifstream file(filename);
  if (!file.is_open()) {
    std::cerr << "Error opening file: " << filename << std::endl;
    return data;
  }

  // read and load the phenotype file
  std::string line;
  while (std::getline(file, line)) {
    std::istringstream lineStream(line);
        std::string key;
        float value1, value2;
	
        if (std::getline(lineStream, key, ',')) {
	  lineStream >> value1;
	  if (lineStream.get() == ',') {
                lineStream >> value2;
                data[key] = std::make_pair(value1, value2);
	  }
        }
  }
  
  file.close();
  return data;
  
}


std::string columnTypeToString(ColumnType type) {
  switch (type) {
    case ColumnType::INT:
      return "INT";
    case ColumnType::FLOAT:
      return "FLOAT";
    case ColumnType::STRING:
      return "STRING";
    case ColumnType::FLAG:
      return "FLAG";
    case ColumnType::GRAPH:
      return "GRAPH";
    default:
      return "UNKNOWN";
  }
}


std::string tokens_to_comma_string(const std::vector<std::string>& input) {
    std::string result;

    if (!input.empty()) {
        result += input.front();
        
        for (auto it = input.begin() + 1; it != input.end(); ++it) {
            result += ',';
            result += *it;
        }
    }

    return result;
}

std::vector<std::string> tokenize_comma_delimited(const std::string& str) {
    std::vector<std::string> tokens;
    size_t start = 0;
    size_t end = str.find(',');

    while (end != std::string::npos) {
        tokens.push_back(str.substr(start, end - start));
        start = end + 1;
        end = str.find(',', start);
    }

    // Add the last token
    tokens.push_back(str.substr(start));

    return tokens;
}

std::string exclude_elements(const std::string_view& str, const std::unordered_set<size_t>& exclude_set) {
  
    std::string result;
    size_t start = 0;
    size_t end = str.find(',');
    
    for (size_t index = 0; end != std::string_view::npos; ++index) {
      std::string_view token = str.substr(start, end - start);
      
      if (exclude_set.find(index) == exclude_set.end()) {
	if (!result.empty()) {
	  result += ',';
	}
	result += token;
      }
      
      start = end + 1;
      end = str.find(',', start);
    }
    
    // Process the last token
    std::string_view last_token = str.substr(start);
    if (exclude_set.find(str.size() - 1) == exclude_set.end()) {
      if (!result.empty()) {
	result += ',';
      }
      result += last_token;
    }
    
    return result;
}

int get_nth_element_as_integer(const std::string_view& str, size_t n) {
  size_t start = 0;
    size_t end = str.find(',');
    
    for (size_t index = 1; index < n; ++index) {
      if (end == std::string_view::npos) {
	throw std::out_of_range("The nth element is out of range.");
      }
      
      start = end + 1;
      end = str.find(',', start);
    }
    
    std::string_view token = (end != std::string_view::npos) ? str.substr(start, end - start) : str.substr(start);
    
    // Check if the token is an integer
    char* end_ptr;
    long value = std::strtol(token.data(), &end_ptr, 10);
    
    if (end_ptr != token.data() && static_cast<size_t>(end_ptr - token.data()) == token.size() && value >= std::numeric_limits<int>::min() && value <= std::numeric_limits<int>::max()) {
        return static_cast<int>(value);
    } else {
        throw std::invalid_argument("The nth element is not an integer.");
    }
}

void get_two_elements_as_floats(const std::string_view& str, size_t n, size_t m,
				float &x, float &y) {
    size_t start = 0;
    size_t end = str.find(',');
    size_t index = 1;

    while (index <= std::max(n, m) && end != std::string_view::npos) {
        std::string_view token = str.substr(start, end - start);

        if (index == n || index == m) {
            char* end_ptr;
            float value = std::strtof(token.data(), &end_ptr);

            if (end_ptr != token.data() && static_cast<size_t>(end_ptr - token.data()) == token.size()) {
                if (index == n) {
                    x = value;
                } else {
                    y = value;
                }
            } else {
                throw std::invalid_argument("Element at index " + std::to_string(index) + " is not a float.");
            }
        }

        start = end + 1;
        end = str.find(',', start);
        ++index;
    }

    if (index <= std::max(n, m)) {
        throw std::out_of_range("One or both indices are out of range.");
    }
}

void write_hdf5_dataframe_attributes(H5::Group& group) {

  H5::DataSpace scalar_dataspace(H5S_SCALAR);
  H5::StrType strdatatype(H5::PredType::C_S1, H5T_VARIABLE);
    
  // write the word "_index"
  H5::Attribute attr = group.createAttribute("_index", strdatatype, scalar_dataspace);
  const char* temp_cstr1 = "_index";
  attr.write(strdatatype, &temp_cstr1);

  // write the encoding-type
  attr = group.createAttribute("encoding-type", strdatatype, scalar_dataspace);
  const char* temp_cstr2 = "dataframe";
  attr.write(strdatatype, &temp_cstr2);

  // write the encoding-version
  attr = group.createAttribute("encoding-version", strdatatype, scalar_dataspace);
  const char* temp_cstr3 = "0.1.0";
  attr.write(strdatatype, &temp_cstr3);
  
}

float euclidean_distance(float x1, float y1, float x2, float y2) {
    float dx = x2 - x1;
    float dy = y2 - y1;
    return std::sqrt(dx*dx + dy*dy);
}

float euclidean_distance_squared(float x1, float y1, float x2, float y2) {
    float dx = x2 - x1;
    float dy = y2 - y1;
    return std::sqrt(dx*dx + dy*dy);
}

bool check_readable(const std::string& filename) {

  std::ifstream file(filename);

  bool answer = file.is_open();
  
  // Close the file if it was opened
  if (file) {
    file.close();
  }

  return answer;
}
