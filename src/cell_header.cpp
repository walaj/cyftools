#include "cell_header.h"

#include <regex>
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
  
  if (ttype == "MA" || ttype == "CA")
    os << "\tIN:" << tag.i;
  
  os << "\t" << tag.data;
  
  return os;
}

void CellHeader::ClearMeta() {
  
  // remove meta tags
  tags.erase(std::remove_if(tags.begin(), tags.end(), 
			    [](const Tag& t) { return t.type == Tag::CA_TAG; }), 
	     tags.end());
  
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
    if (a.type != b.type)
      return a.type < b.type;
    return a.i < b.i;
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

bool CellHeader::isConcatenatable(const CellHeader& header) const {

  std::vector<Tag> this_tags = GetDataTags();
  std::vector<Tag> h_tags    = GetDataTags();

  if (h_tags.size() != this_tags.size()) {
    std::cerr << "Warning: Headers not concatentable. Have MA/CA (data) sizes of: " << this_tags.size() <<
      " and " << h_tags.size() << std::endl;
    return false;
  }
					 
  
  for (size_t i = 0; i < this_tags.size(); i++) {

    const Tag& tt = this_tags.at(i);
    const Tag& ht = h_tags.at(i);
    
    if (tt.id != ht.id || ht.type != tt.type) {
      std::cerr << "Warning: Headers not concatenatable. Have different tags in the " << i <<
	" position - " << tt.id << " and " << ht.id << " with types " <<
	tt.PrintType() << " and " << ht.PrintType() << std::endl;
      return false;
    }
  }

  return true;
  
}

void CellHeader::addTag(const Tag& tag) {
  
  Tag thistag = tag;
  
  // add the index if not already added
  if ( (tag.type == Tag::CA_TAG || tag.type == Tag::MA_TAG) && tag.i == -1 )
    thistag.i = GetDataTags().size();

  // check if tag already in header
  for (const auto& t : GetAllTags()) {
    if (t.id == tag.id && tag.type == t.type && t.type != Tag::PG_TAG) {
      std::cerr << "Warning: Tag " << t.id << " already in header" << std::endl;
      return;
    }
  }

  // find what the new tag id should be
  // this is just the highest existing + 1
  int new_tag_i = 0; 
  for (const auto& t : GetAllTags()) {
    if (t.i > new_tag_i)
      new_tag_i = t.i;
  }
  
  // add tags to data_tags that store
  // float values in the column
  tags.push_back(thistag);

  tags.back().i = new_tag_i + 1;
  
  // if it's a program tag, update the counter for temporal sorting
  //if (tag.type == Tag::PG_TAG) {
  //  tags.back().i = GetProgramTags().size() - 1;
  //} else if (tag.type == Tag::MA_TAG) {
  //  tags.back().i = GetMarkerTags().size() - 1;
  //} else if (tag.type == Tag::CA_TAG) {
  //  tags.back().i = GetMetaTags().size() -1;
  //}

  // make sure tags are sorted
  SortTags();
}

std::ostream& operator<<(std::ostream& os, const CellHeader& h) {
  for (const auto& tag : h.tags) {
    os << tag << std::endl;
  }
  return os;
}

std::vector<Tag> CellHeader::GetMarkerTags() const {

  // return only tags that hold marker data
  std::vector<Tag> marker_tags;
  for (const auto& t : tags)
    if (t.type == Tag::MA_TAG)
      marker_tags.push_back(t);
  
  return marker_tags;
}

std::vector<Tag> CellHeader::GetMetaTags() const {

  // return only tags that hold meta data
  std::vector<Tag> meta_tags;
  for (const auto& t : tags)
    if (t.type == Tag::CA_TAG)
      meta_tags.push_back(t);
  
  return meta_tags;
}

std::vector<Tag> CellHeader::GetProgramTags() const {

  // return only tags that hold marker data
  std::vector<Tag> program_tags;
  for (const auto& t : tags)
    if (t.type == Tag::PG_TAG)
      program_tags.push_back(t);
  
  return program_tags;
}

void CellHeader::Print() const {

  for (const auto& tag : tags) {
    std::cout << tag << std::endl;
  }
  
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
