#pragma once
#include <cstdint>
#include <vector>
#include <string>
#include <unordered_set>

#include "cyf_field.h"   // CyfValueType: typed-column schema

#include <cereal/archives/binary.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/string.hpp> 

class Tag {

 public:
  
  static const uint8_t ANY = 0;
  static const uint8_t MA_TAG = 1; // marker tag
  static const uint8_t CA_TAG = 2; // meta tag
  static const uint8_t GA_TAG = 3; // graph tag
  static const uint8_t PG_TAG = 4; // program tag
  static const uint8_t HD_TAG = 5; // header line: format version + sort order (like SAM @HD)
  static const uint8_t SA_TAG = 6; // sample tag
  static const uint8_t FL_TAG = 7; // flag-bit definition tag (self-describing cflag/pflag bits)
  static const uint8_t IM_TAG = 8; // image/acquisition tag (source TIFF, microns/pixel, microscope)
  static const uint8_t RO_TAG = 9; // region of interest (polygon / shape stored in the header)

  uint8_t type;
  std::string id;
  std::string data;  
  int i = -1; // sort order. -1 means don't factor this in for sort. if > -1, then really store the temporal order of creation for the PG tags
  
  Tag() = default;
  
  Tag(uint8_t type, const std::string& mid, const std::string& mdata) :
    type(type), id(mid), data(mdata), i(-1) {}

  std::string PrintType() const {
    std::string typen;
    switch(type) {
    case MA_TAG: typen = "MA"; break;
    case CA_TAG: typen = "CA"; break;
    case GA_TAG: typen = "GA"; break;
    case PG_TAG: typen = "PG"; break;
    case HD_TAG: typen = "HD"; break;
    case SA_TAG: typen = "SA"; break;
    case FL_TAG: typen = "FL"; break;
    case IM_TAG: typen = "IM"; break;
    case RO_TAG: typen = "RO"; break;
    default: typen = "UNKNOWN TAG TYPE";
    }
    return (typen);
  }

  Tag(uint8_t type, const std::string& mid, const std::string& mdata, int mi) :
    type(type), id(mid), data(mdata), i(mi) {}
  
  Tag(const std::string& line);

  bool isData() const { return type == Tag::MA_TAG || type == Tag::CA_TAG; }

  // ---- typed-column schema accessors (parsed on demand from `data`) ----
  // Extract a "KEY:VALUE" sub-field from the tab-delimited data blob; "" if absent.
  std::string GetField(const std::string& key) const;
  // Declared value type of this column (from TY:), defaulting to Float ('f').
  CyfValueType ValueType() const;
  // Category levels for a TY:A column (from LV:, comma-separated); empty otherwise.
  std::vector<std::string> CategoryLevels() const;

  friend std::ostream& operator<<(std::ostream& os, const Tag& tag);

  template <class Archive>
  void serialize(Archive & ar)
  {
    ar(type, id, i, data);
  }
  
};

class CellHeader {

 public:
  
  CellHeader() = default;

  std::vector<Tag> GetDataTags() const;

  std::vector<Tag> GetInfoTags() const;
  
  std::vector<Tag> GetMarkerTags() const;

  std::vector<Tag> GetProgramTags() const;

  Tag GetVersionTags() const;

  std::vector<Tag> GetSampleTags() const;

  std::vector<Tag> GetFlagTags() const;

  std::vector<Tag> GetMetaTags() const;

  std::vector<Tag> GetAllTags() const { return tags; }

  size_t size() const { return tags.size(); }

  bool isConcatenatable(const CellHeader& header) const;

  void ClearMeta();
  
  void SortTags();
  
  // Function to add a Tag to the header
  void addTag(const Tag& tag);

  // Add a tag, or if one of the same type+id already exists, merge its KEY:VALUE
  // fields in (new keys override, existing keys updated, order preserved). This is
  // what `cyftools addtag` uses so metadata can be filled in incrementally.
  void UpsertTag(const Tag& tag);

  // Append a tag preserving exact order (no sort, no dedup, no warnings).
  // Used by the CYF codec to rebuild a header from on-disk text without
  // disturbing the data-column ordering that record values are bound to.
  void appendRawTag(const Tag& tag) { tags.push_back(tag); }

  // Ensure the header begins with an @HD format/version line (the SAM @HD analog);
  // inserts "@HD  VN:<version>  SO:unsorted" at the front if none is present. Idempotent.
  void EnsureHeaderLine(const std::string& version = "1.0");

  // Read a KEY:VALUE sub-field of the @HD format line (e.g. "MP" microns-per-pixel,
  // "UN" coordinate units); returns "" when @HD or the field is absent.
  std::string GetHeaderField(const std::string& key) const;
  // Set/replace a sub-field on the @HD line, creating @HD if it is missing.
  void SetHeaderField(const std::string& key, const std::string& value);

  void Print() const;

  void PrintMarkers() const;  

  bool validate() const;

  size_t WhichColumn(const std::string& str, uint8_t tag_type) const;

  // The value type of each data column (MA+CA tags, in header order). This is the
  // schema that drives typed record encode/decode. All-float for legacy headers.
  std::vector<CyfValueType> ColumnTypes() const;

  void CleanProgramTags();

  // Remove @RO polygon tags. name_filter: only those whose NM contains it
  // ("" = any); sample_filter: only those with this SA sample id (< 0 = any).
  // Returns the number removed.
  size_t RemoveRoiTags(const std::string& name_filter, long sample_filter);

  void Cut(const std::unordered_set<size_t> to_remove);

  // overload [] operator
  Tag& operator[](std::size_t index) {
    return tags[index];
  }
  
  // overload [] operator for const objects
  const Tag& operator[](std::size_t index) const {
    return tags[index];
  }
  
  // overload at() function for const objects
  const Tag& at(std::size_t index) const {
    return tags.at(index);
  }
  
  // iterators
  std::vector<Tag>::iterator begin() {
    return tags.begin();
  }
  
  std::vector<Tag>::iterator end() {
    return tags.end();
  }
  
  // const iterators
  std::vector<Tag>::const_iterator begin() const {
    return tags.begin();
  }
  
  std::vector<Tag>::const_iterator end() const {
      return tags.end();
  }
  
  std::vector<Tag>::const_iterator cbegin() const {
    return tags.cbegin();
  }
  
  std::vector<Tag>::const_iterator cend() const {
    return tags.cend();
  }
  
  template <class Archive>
    void serialize(Archive & ar)
    {
      ar(tags);
    }

  friend std::ostream& operator<<(std::ostream& os, const CellHeader& h);

private:
  std::vector<Tag> tags;
  
};

