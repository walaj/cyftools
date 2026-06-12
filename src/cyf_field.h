#pragma once
//
// cyf_field.h — the typed value model for CYF columns (docs/CYF_FORMAT.md §3).
//
// A quantification table has a fixed schema: each data column declares one type
// in the header, and every record stores a bare value of that type. CyfField is
// that value; CyfValueType is the type. This is header-only and depends on
// nothing heavy, so it can become part of a portable libcyf.
//
// Backward compatibility: the legacy model stores everything as float (type 'f').
// A float value is just a CyfField of type Float, and asDouble()/asFloat() give
// the numeric view the existing all-float code paths expect.
//
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

// On-disk single-character type codes (BAM-flavored; see spec §3).
enum class CyfValueType : char {
  Int32    = 'i',
  UInt32   = 'I',
  Int64    = 'l',
  UInt64   = 'L',
  Float    = 'f',
  Double   = 'd',
  String   = 'Z',
  Category = 'A',   // small integer index into a header-declared level list (LV:)
  Array    = 'B',   // typed numeric array: a subtype code + values
};

inline char cyfTypeCode(CyfValueType t) { return static_cast<char>(t); }

inline bool cyfIsNumericCode(char c) {
  return c == 'i' || c == 'I' || c == 'l' || c == 'L' || c == 'f' || c == 'd';
}

inline bool cyfIsIntegerCode(char c) {
  return c == 'i' || c == 'I' || c == 'l' || c == 'L';
}

inline CyfValueType cyfTypeFromCode(char c) {
  switch (c) {
    case 'i': return CyfValueType::Int32;
    case 'I': return CyfValueType::UInt32;
    case 'l': return CyfValueType::Int64;
    case 'L': return CyfValueType::UInt64;
    case 'f': return CyfValueType::Float;
    case 'd': return CyfValueType::Double;
    case 'Z': return CyfValueType::String;
    case 'A': return CyfValueType::Category;
    case 'B': return CyfValueType::Array;
    default:
      throw std::runtime_error(std::string("CYF: unknown value type code '") + c + "'");
  }
}

// A single typed value. Kept as an explicit tagged struct (rather than
// std::variant) for clarity and easy reasoning; one payload member is used per
// type. This is a reference-quality representation; it can be made more compact
// later without changing the on-disk format.
struct CyfField {
  CyfValueType type = CyfValueType::Float;

  int64_t      i   = 0;      // integer types (Int32/UInt32/Int64/UInt64)
  double       d   = 0.0;    // Float / Double
  std::string  s;            // String
  uint32_t     cat = 0;      // Category: index into the column's LV: levels

  char                 sub = 'f';  // Array: numeric subtype code
  std::vector<int64_t> ai;         // Array values when sub is an integer code
  std::vector<double>  ad;         // Array values when sub is a float code

  CyfField() = default;

  // ---- factories -------------------------------------------------------
  static CyfField Int(int64_t v, CyfValueType t = CyfValueType::Int64) {
    CyfField f; f.type = t; f.i = v; return f;
  }
  static CyfField Real(double v, CyfValueType t = CyfValueType::Float) {
    CyfField f; f.type = t; f.d = v; return f;
  }
  static CyfField Str(std::string v) {
    CyfField f; f.type = CyfValueType::String; f.s = std::move(v); return f;
  }
  static CyfField Cat(uint32_t idx) {
    CyfField f; f.type = CyfValueType::Category; f.cat = idx; return f;
  }
  static CyfField IntArray(std::vector<int64_t> v, char subtype = 'l') {
    CyfField f; f.type = CyfValueType::Array; f.sub = subtype; f.ai = std::move(v); return f;
  }
  static CyfField RealArray(std::vector<double> v, char subtype = 'f') {
    CyfField f; f.type = CyfValueType::Array; f.sub = subtype; f.ad = std::move(v); return f;
  }

  // ---- numeric view (for the all-float compatibility path) -------------
  double asDouble() const {
    switch (type) {
      case CyfValueType::Int32: case CyfValueType::UInt32:
      case CyfValueType::Int64: case CyfValueType::UInt64:
        return static_cast<double>(i);
      case CyfValueType::Float: case CyfValueType::Double:
        return d;
      case CyfValueType::Category:
        return static_cast<double>(cat);
      default:
        throw std::runtime_error("CYF: value is not numeric (cannot asDouble)");
    }
  }
  float asFloat() const { return static_cast<float>(asDouble()); }

  int64_t asInt64() const {
    switch (type) {
      case CyfValueType::Int32: case CyfValueType::UInt32:
      case CyfValueType::Int64: case CyfValueType::UInt64:
        return i;
      case CyfValueType::Category:
        return static_cast<int64_t>(cat);
      case CyfValueType::Float: case CyfValueType::Double:
        return static_cast<int64_t>(d);
      default:
        throw std::runtime_error("CYF: value is not numeric (cannot asInt64)");
    }
  }

  bool operator==(const CyfField& o) const {
    if (type != o.type) return false;
    switch (type) {
      case CyfValueType::Int32: case CyfValueType::UInt32:
      case CyfValueType::Int64: case CyfValueType::UInt64:
        return i == o.i;
      case CyfValueType::Float: case CyfValueType::Double:
        return d == o.d;
      case CyfValueType::String:
        return s == o.s;
      case CyfValueType::Category:
        return cat == o.cat;
      case CyfValueType::Array:
        return sub == o.sub && ai == o.ai && ad == o.ad;
    }
    return false;
  }
  bool operator!=(const CyfField& o) const { return !(*this == o); }
};
