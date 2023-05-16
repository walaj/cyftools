#include "cell_header.h"

#include <regex>
#include <set>

const std::unordered_set<std::string> Tag::ALLOWED_DATA_TAGS =
  {"MA","ID","GA","CA","FA","XD","YD","ZD"};

const std::unordered_set<std::string> Tag::ALLOWED_INFO_TAGS =
  {"HD","CN","FN"};


Tag::Tag(const std::string& type, const std::string& nm) {

  // tag has to be allowed
  if (!ALLOWED_DATA_TAGS.count(type) && !ALLOWED_INFO_TAGS.count(type)) {
    throw std::runtime_error("Tag with type " + type + " is not allowed");
  }
  
  record_type = type;
  name = nm;
  
}

void CellHeader::Cut(const std::unordered_set<std::string>& include) {

  // build <set> of all tags
  std::unordered_set<std::string> tags_set;
  for (const auto& t : tags)
    tags_set.insert(t.name);
  
  // warn if elements from include list are in header
  for (const auto& token : include) {
    if (tags_set.count(token) == 0) {
      std::cerr << "Warning: '" << token << "' not found in header" << std::endl;
    }
  }

  // Remove elements of tags if their name is not in the vector
  tags.erase(std::remove_if(tags.begin(), tags.end(),
                            [&include](const Tag& tag) {
                              return include.find(tag.GetName()) == include.end();
                            }),
             tags.end());
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
  
  // Remove elements of tags if they are in the tokens vector
  auto new_end = std::remove_if(tags.begin(), tags.end(), [&token, &token_found](const Tag& item) {
    if (item.name == token) {
      token_found = true;
      return true;
    }
    return false;
  });

  // Emit a warning to std::cerr if the token is not found in tags
  if (!token_found) {
    std::cerr << "Warning: '" << token << "' not found in tags" << std::endl;
  } else {
    // Resize col_order to remove the elements that were found in tokens
   tags.resize(new_end - tags.begin());
  }

  
}

void CellHeader::Print() const {
  for (const auto& tag : info_tags)
    std::cout << tag << std::endl;
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

bool Tag::isFlagTag() const {
  bool v = record_type == "FA";

  if (!v)
    return v;
  
  // enforce any requirements for further tags
  //***
  
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

bool Tag::isStringTag() const {
  // right now only accepted string tag is graph
  // but this can be modified in future
  return isGraphTag();
}

bool Tag::isGraphTag() const {
  return record_type == "GA";
}

bool Tag::isDimTag() const {
  return isXDim() || isYDim() || isZDim();
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
  
  if (!isInfoTag() && name.empty())
    throw std::runtime_error("Error: tag " + record_type + " requires NM tag");
  
  return name;
}

void Tag::addValue(const std::string& field, const std::string& value) {
  
  if (field == "NM") {
    name = value;
    return;
  }
  
  values[field] = value;

}

bool Tag::hasName() const {
  return !name.empty();
}

bool CellHeader::hasTag(const std::string& tagname) const {

  for (const auto& t : tags)
    if (tagname == t.GetName())
      return true;

  return false;
}

bool CellHeader::identicalColumns(const CellHeader& header) const {

  size_t i = 0;

  if (this->tags.size() != header.tags.size())
    return false;

  for (size_t i = 0; i < this->tags.size(); i++) {
    if (this->tags.at(i).GetName() != header.tags.at(i).GetName())
      return false;
  }

  return true;
  
  
  
}

void CellHeader::addTag(const Tag& tag) {

  /// check that tags has a name if its a column tag
  if (!tag.isInfoTag()) {
    if (!tag.hasName()) {
      throw std::runtime_error("Tag of record type " + tag.record_type + " needs an NM:<name> field");
    }
  }

  // throw error if already exists
  if (hasTag(tag.GetName())) {
    throw std::runtime_error("Trying to add tag " + tag.GetName() + " that already exists in header");
  }

  // add the tag
  if (tag.isInfoTag())
    info_tags.push_back(tag);
  else if (tag.isColumnTag())
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
  } else if (tag.isInfoTag()) {
    ;//
  } else if (tag.isGraphTag()) {
    graph_.insert(tag.GetName());
  } else if (tag.isFlagTag()) {
    flag_.insert(tag.GetName());
  } else {
    throw std::runtime_error("Tag: " + tag.record_type + " - must be one of: HD, MA, CA, XD, YD, ZD, ID");
  }
}

std::ostream& operator<<(std::ostream& os, const Tag& tag) {
  os << "@" << tag.record_type;
  if (!tag.name.empty())
    os << "\tNM:" << tag.name;
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

bool CellHeader::hasFlag(const std::string& m) const {
  return (flag_.count(m));
}

bool CellHeader::hasGraph(const std::string& m) const {
  return (graph_.count(m));
}

std::ostream& operator<<(std::ostream& os, const CellHeader& h) {
  for (const auto& k : h.tags)
    os << k << std::endl;
  return os;
}

size_t CellHeader::whichColumn(const std::string& str) const {

  size_t count = 0;
  for (const auto& t : tags) {
    if (t.isColumnTag()) {
      if (t.GetName() == str)
	return count;
      count++;
    }
  }
  
  std::cerr << "Column: " << str << " never found in header" << std::endl;
  return static_cast<size_t>(-1);

}

size_t CellHeader::whichMarkerColumn(const std::string& str) const {

  size_t count = 0;
  for (const auto& t : tags) {
    if (t.isColumnTag() && t.isMarkerTag()) {
      if (t.GetName() == str)
	return count;
      count++;
    }
  }
  
  std::cerr << "Column: " << str << " never found in header" << std::endl;
  return static_cast<size_t>(-1);
  
}
