#ifndef CELL_HEADER_H
#define CELL_HEADER_H

#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>
#include <set>
#include <sstream>

struct Tag {
  
  std::string record_type; // e.g., "MA", "CD"
  std::unordered_map<std::string, std::string> values; // Stores the tag fields and their data
  std::string name;

  // Default constructor
  Tag() = default;
  
  // Constructor to initialize the structure from a given line
  Tag(const std::string& line);

  // Constructor with tag name and NM
  Tag(const std::string& type, const std::string& nm);
  
  // Constructor to initialize the structure from a given line
  // Function to add a tag field and its data to the structure
  void addValue(const std::string& field, const std::string& value);

  // Overloaded << operator
  friend std::ostream& operator<<(std::ostream& os, const Tag& tag);

  bool isMarkerTag() const;

  bool isMetaTag() const;
  
  bool isXDim() const;
  bool isYDim() const;
  bool isZDim() const;
  bool isIDTag() const;    
  bool isOrderTag() const;
  bool isVersionTag() const;
  
  std::string GetName() const;

private:

  std::istream& safeGetline(std::istream& is, std::string& line, char delim) {
    std::getline(is >> std::ws, line, delim);
    if (!is && line.empty()) {
      throw std::runtime_error("Error: empty or invalid input.");
    }
    return is;
  }

};

class CellHeader {

public:
  // Constructor
  CellHeader() = default;
  
  // Function to add a Tag to the header
  void addTag(const Tag& tag);
  
  // Function to access the tags
  const std::vector<Tag>& getTags() const {
    return tags;
  }

  void Cut(const std::set<std::string>& tokens);
  
  void Remove(const std::string& token);
  
  void Print() const;

  const std::vector<std::string>& GetColOrder() const { return col_order; }
  
  friend std::ostream& operator<<(std::ostream& os, const CellHeader& h);

  bool hasMarker(const std::string& m) const;
  bool hasMeta(const std::string& m) const;  

  std::string GetX() const { return x_; }
  std::string GetY() const { return y_; }
  std::string GetZ() const { return z_; }
  std::string GetID() const { return id_; }


  
private:

  friend class CellTable;
  
  std::vector<Tag> tags; // Vector to store the Tag objects

  std::string version_;
  std::string id_;
  std::string x_;
  std::string y_;
  std::string z_;  
  std::set<std::string> markers_;
  std::set<std::string> meta_;  

  std::vector<std::string> col_order;
  
};

#endif
