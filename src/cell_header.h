#ifndef CELL_HEADER_H
#define CELL_HEADER_H

#include <iostream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <set>
#include <sstream>

struct Tag {

  static const std::unordered_set<std::string> ALLOWED_DATA_TAGS;
  
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

  bool hasName() const;

  bool isColumnTag() const {
    return ALLOWED_DATA_TAGS.count(record_type); 
  }

  // Overloaded << operator
  friend std::ostream& operator<<(std::ostream& os, const Tag& tag);

  bool isMarkerTag() const;
  bool isMetaTag() const;
  bool isGraphTag() const;  
  bool isStringTag() const;
  
  bool isXDim() const;
  bool isYDim() const;
  bool isZDim() const;
  bool isDimTag() const;
  bool isFlagTag() const;
  bool isIDTag() const;    
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
  
  const std::vector<Tag>& GetTags() const { return tags; }
  
  // Function to add a Tag to the header
  void addTag(const Tag& tag);

  // Function to access the tags
  const std::vector<Tag>& getTags() const {
    return tags;
  }

  void Cut(const std::unordered_set<std::string>& include);
  
  void Remove(const std::string& token);
  
  void Print() const;

  friend std::ostream& operator<<(std::ostream& os, const CellHeader& h);

  size_t whichMarkerColumn(const std::string& str) const;

  size_t whichColumn(const std::string& str) const;
  
  bool hasMarker(const std::string& m) const;
  bool hasMeta(const std::string& m) const;  
  bool hasFlag(const std::string& m) const;
  bool hasGraph(const std::string& m) const;  

  std::string GetX() const { return x_; }
  std::string GetY() const { return y_; }
  std::string GetZ() const { return z_; }
  std::string GetID() const { return id_; }

  const Tag& GetColumnTag(int i) const {

    size_t count = 0;
    for (const auto& t : tags)
      if (t.isColumnTag()) {
	if (count == i)
	  return t;
	count++;
      }
    assert(false);
    return tags[0]; // never will get called
  }
  
  bool validate() const {

    // collect the columns
    std::unordered_set<std::string> columns;
    for (const auto& t : tags)
      if (t.isColumnTag())
	columns.insert(t.GetName());

    std::string x = GetX();
    if (!x.empty() && !columns.count(x)) {
      std::cerr << "Header error: can't find XD column " << x << " in header " << std::endl;
      return false;
    }

    std::string y = GetY();
    if (!y.empty() && !columns.count(y)) {
      std::cerr << "Header error: can't find YD column " << y << " in header " << std::endl;
      return false;
    }
    
    std::string id = GetID();
    if (!id.empty() && !columns.count(id)) {
      std::cerr << "Header error: can't find ID column " << id << " in header " << std::endl;
      return false;
    }
    
    return true;
  }

  std::vector<std::string> GetColOrder() const {

    std::vector<std::string> col_order;
    for (const auto& t : tags)
      if (t.isColumnTag())
	col_order.push_back(t.GetName());
    return col_order;
    
  }
  
  int ColumnCount() const {
    int out = 0;
    for (const auto& t : tags)
      if (t.isColumnTag())
	out++;
    return out;
  }
  
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
  std::set<std::string> graph_;
  std::set<std::string> flag_;
  
};

#endif
