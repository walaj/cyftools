#pragma once
//
// cyf_flags.h — the built-in cflag (structural) flag vocabulary, plus helpers to
// declare it inside a CYF header as self-describing @FL tags.
//
// This is the first step off the compile-time #define flag table in cysift.h: the
// code still uses the #defines internally, but every CYF file can now declare what
// each cflag bit means (e.g. "@FL  ID:TLS  RG:cflag  BI:5"), so a reader needs no
// source to interpret a cell's flags. Phenotype (pflag) bits are panel-specific
// (the ORION_*/PROSTATE_* tables), and are expected to be declared per-dataset
// rather than hardcoded here — that is the next step.
//
#include <string>
#include <utility>
#include <vector>

#include "cell_header.h"

namespace cyf {

struct FlagDef {
  const char* reg;   // "cflag" or "pflag"
  int         bit;   // 0-based bit index
  const char* name;  // symbolic name
};

// The structural cflag bits — mirrors the *_FLAG #defines in cysift.h as bit indices
// (TUMOR_FLAG=1=bit0, MARK_FLAG=2=bit1, ... BUILD_GRAPH_FLAG=2097152=bit21).
inline const std::vector<FlagDef>& standardFlags() {
  static const std::vector<FlagDef> defs = {
    {"cflag",  0, "Tumor"},
    {"cflag",  1, "Mark"},
    {"cflag",  2, "Margin"},
    {"cflag",  3, "TumorManual"},
    {"cflag",  4, "Tcell"},
    {"cflag",  5, "TLS"},
    {"cflag",  6, "MarginManual"},
    {"cflag",  7, "Normal"},
    {"cflag",  8, "NormalMargin"},
    {"cflag",  9, "Artifact"},
    {"cflag", 21, "BuildGraph"},
  };
  return defs;
}

// One @FL tag: id = name, data = "RG:<reg>\tBI:<bit>". The sort index is the bit
// number, so SortTags orders the @FL block by ascending bit.
inline Tag makeFlagTag(const FlagDef& f) {
  Tag t(Tag::FL_TAG, f.name,
        std::string("RG:") + f.reg + "\tBI:" + std::to_string(f.bit));
  t.i = f.bit;
  return t;
}

// Append the standard cflag @FL declarations to a header. Skips any (register, bit)
// already declared, so it is safe to call more than once.
inline void addStandardFlagTags(CellHeader& h) {
  std::vector<std::pair<std::string, std::string>> existing;
  for (const auto& t : h.GetFlagTags())
    existing.emplace_back(t.GetField("RG"), t.GetField("BI"));

  for (const auto& f : standardFlags()) {
    const std::string bi = std::to_string(f.bit);
    bool present = false;
    for (const auto& e : existing)
      if (e.first == f.reg && e.second == bi) { present = true; break; }
    if (!present)
      h.appendRawTag(makeFlagTag(f));
  }
}

} // namespace cyf
