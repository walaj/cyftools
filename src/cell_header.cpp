#include "cell_header.h"


Tag::Tag(const std::string& line) {
  std::istringstream iss(line);
  std::string token;
  
  // Parse the record type
  getline(iss, token, '\t');
  record_type = token.substr(1); // Remove the '@' character
  
  // Parse the tag fields and their data
  while (getline(iss, token, '\t')) {
    std::string field = token.substr(0, 2);
    std::string value = token.substr(3);
    values[field] = value;
  }
}

void Tag::addValue(const std::string& field, const std::string& value) {
  values[field] = value;
}

std::ostream& operator<<(std::ostream& os, const Tag& tag) {
  os << "@" << tag.record_type;
  for (const auto& kv : tag.values) {
    os << "\t" << kv.first << ":" << kv.second;
  }
  return os;

}

std::ostream& operator<<(std::ostream& os, const CellHeader& h) {
  for (const auto& k : h.tags)
    os << k << std::endl;
  return os;
}
