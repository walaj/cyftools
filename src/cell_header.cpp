#include "cell_header.h"

#include <iostream>

// ADD: need to fill out
bool CellHeader::validate() const {
  return true;
}

std::ostream& operator<<(std::ostream& os, const Tag& tag) {

  // get string representation of tag type
  std::string ttype;
  switch (tag.type) {
  case Tag::MA_TAG: ttype = "MA"; break;
  case Tag::CA_TAG: ttype = "CA"; break;
  case Tag::GA_TAG: ttype = "GA"; break;
  case Tag::PG_TAG: ttype = "PG"; break;
  default:
    throw std::runtime_error("Unknown tag type with uint8_t value: " + std::to_string(tag.type));
  }
  os << "@" << ttype;

  if (ttype != "PG")
    os << "\tID:" << tag.id;
  
  os << "\t" << tag.data;
  
  return os;
}


std::vector<Tag> CellHeader::GetDataTags() const {

  std::vector<Tag> data_tags;
  for (const auto& t : tags)
    if (t.type == Tag::MA_TAG || t.type == Tag::CA_TAG)
      data_tags.push_back(t);
      
  return data_tags;
}

std::vector<Tag> CellHeader::GetInfoTags() const {

  std::vector<Tag> info_tags;
  for (const auto& t : tags)
    if (t.type == Tag::GA_TAG || t.type == Tag::PG_TAG)
      info_tags.push_back(t);
      
  return info_tags;
  
  return info_tags;
}

void CellHeader::SortTags() {

  // we are using the order provided by the tag definitions
  std::sort(tags.begin(), tags.end(), [](const Tag &a, const Tag &b) {
    return a.type < b.type;
  });

}

void CellHeader::Cut(const std::unordered_set<size_t> to_remove) {
  
  std::vector<Tag> new_tags;
  for (size_t i = 0; i < tags.size(); i++) {
    if (to_remove.count(i))
      continue;
    new_tags.push_back(tags.at(i));
  }
  tags = new_tags;
}

void CellHeader::addTag(const Tag& tag) {
  
  // ADD: check tag type is valid
  ///

  // add tags to data_tags that store
  // float values in the column
  tags.push_back(tag);

  // info tags that don't store column data
  //if (tag.type == Tag::PG_TAG ||
  //    tag.type == Tag::GA_TAG) {
  //  info_tags.push_back(tag);
  //  return;
  //}
}


std::vector<Tag> CellHeader::GetMarkerTags() const {

  // return only tags that hold marker data
  std::vector<Tag> marker_tags;
  for (const auto& t : tags)
    if (t.type == Tag::MA_TAG)
      marker_tags.push_back(t);
  
  return marker_tags;
}

/*const Tag& CellHeader::GetDataTag(int i) const {
  assert(i < data_tags.size());
  return data_tags.at(i);
  }*/


void CellHeader::Print() const {

  for (const auto& tag : tags) {
    std::cout << tag << std::endl;
  }
  //  for (const auto& tag : info_tags)
  //  std::cout << tag << std::endl;
  
}

size_t CellHeader::WhichColumn(const std::string& str, uint8_t tag_type) const {

  // ADD: error check on tag_type to make sure allowed
  ////

  // find the tag index
  size_t count = 0;
  for (const auto& t : tags) {
    if (t.type == tag_type) {
      if (t.id == str)
	return count;
      count++;
    }
  }
  
  throw std::runtime_error("Column: " + str + " not in header");
  return static_cast<size_t>(-1);
  
}

Tag::Tag(const std::string& line) {

  // setup a regex to parse on spaces
  std::regex ws_re("\\s+");
  std::sregex_token_iterator it(line.begin(), line.end(), ws_re, -1); 
  std::sregex_token_iterator end;
  
  if (it != end) {
    std::string token = *it;

    if      (token.substr(1) == "MA") { type = Tag::MA_TAG; }
    else if (token.substr(1) == "GA") { type = Tag::GA_TAG; }
    else if (token.substr(1) == "CA") { type = Tag::CA_TAG; }
    else if (token.substr(1) == "PG") { type = Tag::PG_TAG; }
    else { throw std::runtime_error("Tag of type" + token.substr(1) + "not an allowed tag"); }
    
    ++it;
  }
  
  // Parse the tag fields and their data
  for (; it != end; ++it) {
    std::string token = *it;
    std::string field = token.substr(0, 2);
    std::string value = token.substr(3);

    // special case for ID
    if (field == "ID") {
      id = value;
    } else {
      data = data + field + ":" + value;
    }
  }
}
