#include "cell_header.h"

#include <regex>
#include <set>

Tag::Tag(const std::string& type, const std::string& nm) {

  record_type = type;
  values["NM"] = nm;
  
}

void CellHeader::Cut(const std::set<std::string>& tokens) {
  
  std::set<std::string> col_order_set(col_order.begin(), col_order.end());

  std::set<std::string> tags_set;
  
  for (const auto& t : tags)
    tags_set.insert(t.name);
  
  // Remove elements from m_table and col_order if they are not in tokens
  for (const auto& token : tokens) {
    if (col_order_set.count(token) == 0) {
      std::cerr << "Warning: '" << token << "' not found in col_order" << std::endl;
    }
    
    if (tags_set.count(token) == 0) {
      std::cerr << "Warning: '" << token << "' not found in tags" << std::endl;
    }
  }

  // Remove elements of tags if their name is not in the tokens vector
  auto new_end = std::remove_if(tags.begin(), tags.end(), [&tokens](const Tag& tag) {
    return std::find(tokens.begin(), tokens.end(), tag.name) == tokens.end();
  });
  
  // Resize tags to remove the elements whose names are not in tokens
  tags.resize(new_end - tags.begin());

  // remove col_order elements
  col_order.erase(std::remove_if(col_order.begin(), col_order.end(), [&tokens](const std::string& item) {
    return std::find(tokens.begin(), tokens.end(), item) == tokens.end();
  }), col_order.end());
}

void CellHeader::Remove(const std::string& token) {

  bool token_found = false;
      
  if (markers_.count(token) > 0) {
    markers_.erase(token);
  } else if (meta_.count(token) > 0) {
    meta_.erase(token);
  } else {
    std::cerr << "Warning: '" << token << "' not found in the markers or meta of table" << std::endl;
  }
  
  // Remove elements of col_order if they are in the tokens vector
  auto new_end = std::remove_if(col_order.begin(), col_order.end(), [&token, &token_found](const std::string& item) {
    if (item == token) {
      token_found = true;
      return true;
    }
    return false;
  });

  // Emit a warning to std::cerr if the token is not found in col_order
  if (!token_found) {
    std::cerr << "Warning: '" << token << "' not found in col_order" << std::endl;
  } else {
    // Resize col_order to remove the elements that were found in tokens
    col_order.resize(new_end - col_order.begin());
  }

  
}



void CellHeader::Print() const {
  for (const auto& tag : tags) 
    std::cout << tag << std::endl;
}

Tag::Tag(const std::string& line) {

  std::regex ws_re("\\s+");
  std::sregex_token_iterator it(line.begin(), line.end(), ws_re, -1); 
  std::sregex_token_iterator end;
  
  if (it != end) {
    std::string token = *it;
    record_type = token.substr(1); // Remove the '@' character
    ++it;
  }

  // Parse the tag fields and their data
  for (; it != end; ++it) {
    std::string token = *it;
    std::string field = token.substr(0, 2);
    std::string value = token.substr(3);
    this->addValue(field, value);
  }
}

bool Tag::isOrderTag() const {
  bool v = record_type == "OD";

  if (!v)
    return v;

  return v;
}


bool Tag::isVersionTag() const {
  bool v = record_type == "HD";

  if (!v)
    return v;
  
  if (values.empty()) {
    throw std::runtime_error("Error: @HD tag malformed, needs VN key");
  }
  
  const auto& iter = values.find("VN");
  if (iter == values.end()) {
    throw std::runtime_error("Error: 'VN' required tag of @HD");
  }

  return v;

}

bool Tag::isMarkerTag() const {
  return record_type == "MA";
}

bool Tag::isMetaTag() const {
  return record_type == "CA";
}

bool Tag::isXDim() const {
  return record_type == "XD";
}

bool Tag::isYDim() const {
  return record_type == "YD";
}

bool Tag::isZDim() const {
  return record_type == "ZD";
}

bool Tag::isIDTag() const {
  return record_type == "ID";
}


std::string Tag::GetName() const {
  if (values.empty()) {
    
    throw std::runtime_error("Error: tag " + record_type + "is empty, needs NM tag");
  }

  const auto& iter = values.find("NM");
  if (iter == values.end()) {
    throw std::runtime_error("Error: tag " + record_type + "requires NM tag");
  }
  
  return iter->second;
}

void Tag::addValue(const std::string& field, const std::string& value) {
  values[field] = value;

  if (field == "NM")
    name = value;
}

void CellHeader::addTag(const Tag& tag) {
    tags.push_back(tag);
  
  if (tag.isMarkerTag()) {
    markers_.insert(tag.GetName());
  } else if (tag.isMetaTag()) {
    meta_.insert(tag.GetName());
  } else if (tag.isXDim()) {
    x_ = tag.GetName();
  } else if (tag.isYDim()) {
    y_ = tag.GetName();
  } else if (tag.isZDim()) {
    z_ = tag.GetName();
  } else if (tag.isIDTag()) {
    id_ = tag.GetName();
  } else if (tag.isVersionTag()) {
    ;//
  } else {
    throw std::runtime_error("Tag: " + tag.record_type + " - must be one of: HD, MA, CA, XD, YD, ZD, ID");
  }
  
}


std::ostream& operator<<(std::ostream& os, const Tag& tag) {
  os << "@" << tag.record_type;
  for (const auto& kv : tag.values) {
    os << "\t" << kv.first << ":" << kv.second;
  }
  return os;

}

bool CellHeader::hasMarker(const std::string& m) const {
  return (markers_.count(m));
}

bool CellHeader::hasMeta(const std::string& m) const {
  return (meta_.count(m));
}

std::ostream& operator<<(std::ostream& os, const CellHeader& h) {
  for (const auto& k : h.tags)
    os << k << std::endl;
  return os;
}
