#ifndef CELL_HEADER_H
#define CELL_HEADER_H

#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>
#include <sstream>

struct Tag {
  
  std::string record_type; // e.g., "MA", "CD"
  std::unordered_map<std::string, std::string> values; // Stores the tag fields and their data
  
  // Constructor to initialize the structure from a given line
  Tag(const std::string& line);
  
  // Constructor to initialize the structure from a given line
  // Function to add a tag field and its data to the structure
  void addValue(const std::string& field, const std::string& value);

  // Overloaded << operator
  friend std::ostream& operator<<(std::ostream& os, const Tag& tag);
  
};

class CellHeader {

public:
  // Constructor
  CellHeader() = default;
  
  // Function to add a Tag to the header
  void addTag(const Tag& tag) {
    tags.push_back(tag);
  }
  
  // Function to access the tags
  const std::vector<Tag>& getTags() const {
    return tags;
  }

  friend std::ostream& operator<<(std::ostream& os, const CellHeader& h);
  
private:
  std::vector<Tag> tags; // Vector to store the Tag objects
  
};

#endif
