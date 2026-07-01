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
  case Tag::HD_TAG: ttype = "HD"; break;
  case Tag::SA_TAG: ttype = "SA"; break;
  case Tag::FL_TAG: ttype = "FL"; break;
  case Tag::IM_TAG: ttype = "IM"; break;
  case Tag::RO_TAG: ttype = "RO"; break;
  default:
    throw std::runtime_error("Unknown tag type with uint8_t value: " + std::to_string(tag.type));
  }
  os << "@" << ttype;

  if (!tag.id.empty())
    os << "\tID:" << tag.id;

  if (!tag.data.empty())
    os << "\t" << tag.data;

  return os;
}

void CellHeader::ClearMeta() {
  
  // remove meta tags
  tags.erase(std::remove_if(tags.begin(), tags.end(), 
			    [](const Tag& t) { return t.type == Tag::CA_TAG; }), 
	     tags.end());
  
}

void CellHeader::EnsureHeaderLine(const std::string& version) {
  for (const auto& t : tags)
    if (t.type == Tag::HD_TAG)
      return;                                  // already present
  Tag hd(Tag::HD_TAG, "", "VN:" + version + "\tSO:unsorted");
  tags.insert(tags.begin(), hd);               // front, so it is the first line
}

std::string CellHeader::GetHeaderField(const std::string& key) const {
  for (const auto& t : tags)
    if (t.type == Tag::HD_TAG)
      return t.GetField(key);
  return "";
}

void CellHeader::SetHeaderField(const std::string& key, const std::string& value) {
  EnsureHeaderLine();
  // @HD has an empty id, so this merges the field into the existing @HD line.
  UpsertTag(Tag(Tag::HD_TAG, "", key + ":" + value));
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
    if (t.type == Tag::GA_TAG || t.type == Tag::PG_TAG ||
	t.type == Tag::HD_TAG || t.type == Tag::SA_TAG ||
	t.type == Tag::FL_TAG || t.type == Tag::IM_TAG ||
	t.type == Tag::RO_TAG)
      info_tags.push_back(t);
      
  return info_tags;
}

void CellHeader::SortTags() {

  // we are using the order provided by the tag definitions
  std::sort(tags.begin(), tags.end(), [](const Tag &a, const Tag &b) {
    // @HD (the format/version line) is always first
    if ((a.type == Tag::HD_TAG) != (b.type == Tag::HD_TAG))
      return a.type == Tag::HD_TAG;
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

void CellHeader::UpsertTag(const Tag& tag) {

  // find an existing tag of the same class + id and merge into it
  for (auto& t : tags) {
    if (t.type == tag.type && t.id == tag.id) {

      // collect KEY:VALUE tokens: existing first (order preserved), then the new
      // ones, which override a matching KEY or append at the end.
      std::vector<std::pair<std::string, std::string>> kv;  // (key, "KEY:VALUE")
      auto absorb = [&kv](const std::string& blob) {
        size_t start = 0;
        while (start <= blob.size()) {
          size_t end = blob.find('\t', start);
          std::string tok = (end == std::string::npos) ? blob.substr(start)
                                                       : blob.substr(start, end - start);
          if (!tok.empty()) {
            std::string key = tok.substr(0, tok.find(':'));
            bool found = false;
            for (auto& p : kv) if (p.first == key) { p.second = tok; found = true; break; }
            if (!found) kv.emplace_back(key, tok);
          }
          if (end == std::string::npos) break;
          start = end + 1;
        }
      };
      absorb(t.data);     // existing fields
      absorb(tag.data);   // new fields override

      std::string merged;
      for (size_t i = 0; i < kv.size(); ++i) { if (i) merged += "\t"; merged += kv[i].second; }
      t.data = merged;
      return;
    }
  }

  // none matched: append a fresh tag (info tags sort by class with i = -1)
  tags.push_back(tag);
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

std::vector<Tag> CellHeader::GetSampleTags() const {

  // return only tags that hold marker data
  std::vector<Tag> sample_tags;
  for (const auto& t : tags)
    if (t.type == Tag::SA_TAG)
      sample_tags.push_back(t);

  return sample_tags;
}

std::vector<Tag> CellHeader::GetFlagTags() const {

  // return only the @FL flag-bit definition tags
  std::vector<Tag> flag_tags;
  for (const auto& t : tags)
    if (t.type == Tag::FL_TAG)
      flag_tags.push_back(t);

  return flag_tags;
}


void CellHeader::Print() const {

  for (const auto& tag : tags) {
    std::cout << tag << std::endl;
  }
  
}

void CellHeader::PrintMarkers() const {

  auto markers = GetMarkerTags();

  for (size_t i = 0; i < markers.size(); ++i) {
    std::cout << markers[i].id << std::endl;
    //    if (i < markers.size() - 1) {
    //  std::cout << std::endl;
    //}
  }

  return;
  
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
    else if (token.substr(1) == "SA") { type = Tag::SA_TAG; }
    else if (token.substr(1) == "HD") { type = Tag::HD_TAG; }
    else if (token.substr(1) == "FL") { type = Tag::FL_TAG; }
    else if (token.substr(1) == "IM") { type = Tag::IM_TAG; }
    else if (token.substr(1) == "RO") { type = Tag::RO_TAG; }
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

// ---------------------------------------------------------------------------
// Typed-column schema accessors. The `data` blob is the tab-delimited remainder
// of a header line ("TY:f\tCH:1\tLV:a,b,c"), so sub-fields are parsed from it on
// demand. Columns with no TY: field default to Float, keeping legacy (all-float)
// headers fully backward compatible.
// ---------------------------------------------------------------------------

std::string Tag::GetField(const std::string& key) const {
  size_t start = 0;
  while (start <= data.size()) {
    size_t end = data.find('\t', start);
    std::string tok = (end == std::string::npos) ? data.substr(start)
                                                  : data.substr(start, end - start);
    size_t colon = tok.find(':');
    if (colon != std::string::npos && tok.substr(0, colon) == key)
      return tok.substr(colon + 1);
    if (end == std::string::npos) break;
    start = end + 1;
  }
  return "";
}

CyfValueType Tag::ValueType() const {
  std::string ty = GetField("TY");
  if (ty.empty()) return CyfValueType::Float;   // legacy default: everything is float
  return cyfTypeFromCode(ty[0]);
}

std::vector<std::string> Tag::CategoryLevels() const {
  std::vector<std::string> levels;
  std::string lv = GetField("LV");
  if (lv.empty()) return levels;
  size_t start = 0;
  while (true) {
    size_t comma = lv.find(',', start);
    if (comma == std::string::npos) { levels.push_back(lv.substr(start)); break; }
    levels.push_back(lv.substr(start, comma - start));
    start = comma + 1;
  }
  return levels;
}

std::vector<CyfValueType> CellHeader::ColumnTypes() const {
  std::vector<CyfValueType> types;
  for (const auto& t : tags)
    if (t.type == Tag::MA_TAG || t.type == Tag::CA_TAG)
      types.push_back(t.ValueType());
  return types;
}

void CellHeader::CleanProgramTags() {
  // remove PG tags
  tags.erase(std::remove_if(tags.begin(), tags.end(),
			    [](const Tag& t) { return t.type == Tag::PG_TAG; }),
	     tags.end());

}

size_t CellHeader::RemoveRoiTags(const std::string& name_filter, long sample_filter,
                                 const std::string& id_filter) {
  const size_t before = tags.size();
  tags.erase(std::remove_if(tags.begin(), tags.end(),
    [&](const Tag& t) {
      if (t.type != Tag::RO_TAG) return false;
      if (!id_filter.empty() && t.id != id_filter)   // exact ID match (disambiguates dup names)
        return false;
      if (!name_filter.empty() && t.GetField("NM").find(name_filter) == std::string::npos)
        return false;
      if (sample_filter >= 0) {
        try { if (std::stol(t.GetField("SA")) != sample_filter) return false; }
        catch (...) { return false; }
      }
      return true;   // matches -> remove
    }), tags.end());
  return before - tags.size();
}
