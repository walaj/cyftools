// cyf_test.cpp — dependency-light self-test for the CYF codec and the
// cereal/CYF compatibility adapter. Builds with just cereal (no mlpack/CGAL/etc):
//
//     make cyf-test          (from src/)
//   or
//     g++ -std=c++17 -DUSE_64_BIT -I. -I../external/cereal/include \
//         cyf_test.cpp cyf_io.cpp cell_header.cpp -o cyf_test && ./cyf_test
//
// Covers: binary round-trip, header text round-trip, format auto-detection,
// truncation detection, and round-tripping through both backends of OutArchive/
// InArchive. Exit code 0 on success.
#include "cyf_io.h"
#include "cell_archive.h"
#include "cyf_flags.h"
#include <cassert>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

static int failures = 0;
static void check(bool ok, const char* what) {
  std::cout << "  [" << (ok ? "PASS" : "FAIL") << "] " << what << "\n";
  if (!ok) ++failures;
}

static Tag mk(uint8_t t, std::string id, std::string d, int i) { Tag x(t, id, d); x.i = i; return x; }

static CellHeader buildHeader() {
  CellHeader h; int i = 0;
  h.appendRawTag(mk(Tag::HD_TAG, "",     "VN:1.0\tSO:unsorted",                    i++));
  h.appendRawTag(mk(Tag::SA_TAG, "S1",   "SM:tonsil_01\tPA:Orion",                 i++));
  h.appendRawTag(mk(Tag::MA_TAG, "DAPI", "TY:f\tCH:1",                             i++));
  h.appendRawTag(mk(Tag::MA_TAG, "CD3",  "TY:f\tCH:5",                             i++));
  h.appendRawTag(mk(Tag::CA_TAG, "area", "TY:f",                                   i++));
  h.appendRawTag(mk(Tag::CA_TAG, "ecc",  "TY:f",                                   i++));
  h.appendRawTag(mk(Tag::PG_TAG, "pg1",  "PN:cyftools\tVN:2.0\tCL:convert in.csv", i++));
  return h;
}
static std::vector<Cell> buildCells() {
  Cell a; a.id = 1; a.x = 100.5f; a.y = 200.25f; a.cflag = 1; a.pflag = 2; a.cols = {12.0f, 3.5f, 42.0f, 0.8f};
  Cell b; b.id = 2; b.x = 101.0f; b.y = 202.0f;  b.cflag = 0; b.pflag = 0; b.cols = {8.0f, 0.0f, 37.0f, 0.6f};
  return {a, b};
}
static bool cellEq(const Cell& a, const Cell& b) {
  return a.id == b.id && a.x == b.x && a.y == b.y && a.cflag == b.cflag &&
         a.pflag == b.pflag && a.cols == b.cols;
}

// A header exercising every column type: float marker + int/string/categorical/
// array metadata columns. Data columns are MA+CA (what the schema binds to).
static CellHeader buildTypedHeader() {
  CellHeader h; int i = 0;
  h.appendRawTag(mk(Tag::HD_TAG, "",         "VN:1.0",                              i++));
  h.appendRawTag(mk(Tag::MA_TAG, "DAPI",     "TY:f",                                i++));
  h.appendRawTag(mk(Tag::CA_TAG, "count",    "TY:i",                                i++));
  h.appendRawTag(mk(Tag::CA_TAG, "label",    "TY:Z",                                i++));
  h.appendRawTag(mk(Tag::CA_TAG, "celltype", "TY:A\tLV:tumor,stroma,immune",        i++));
  h.appendRawTag(mk(Tag::CA_TAG, "knn",      "TY:B\tKD:knn",                        i++));
  return h;
}
static cyf::CyfRow typedRow(uint64_t id, double dapi, int64_t cnt, std::string lab,
                            uint32_t ct, std::vector<int64_t> knn) {
  cyf::CyfRow r; r.id = id; r.x = 1.0f; r.y = 2.0f; r.cflag = 5; r.pflag = 6;
  r.values = { CyfField::Real(dapi, CyfValueType::Float),
               CyfField::Int(cnt, CyfValueType::Int32),
               CyfField::Str(std::move(lab)),
               CyfField::Cat(ct),
               CyfField::IntArray(std::move(knn), 'L') };
  return r;
}
static bool rowEq(const cyf::CyfRow& a, const cyf::CyfRow& b) {
  if (a.id != b.id || a.x != b.x || a.y != b.y || a.cflag != b.cflag || a.pflag != b.pflag) return false;
  return a.values == b.values;
}

int main() {
  const CellHeader h = buildHeader();
  const std::vector<Cell> cells = buildCells();

  std::cout << "CYF codec:\n";

  // round-trip through CyfWriter/CyfReader
  {
    std::ostringstream os(std::ios::binary);
    { cyf::CyfWriter w(os); w.writeHeader(h); for (auto& c : cells) w.writeCell(c); }
    std::istringstream is(os.str(), std::ios::binary);
    cyf::CyfReader r(is);
    CellHeader h2;
    check(r.readHeader(h2), "readHeader");
    check(cyf::renderHeaderText(h2) == cyf::renderHeaderText(h), "header text round-trip");
    std::vector<Cell> got; Cell c;
    cyf::CyfReader::Status st;
    while ((st = r.readCell(c)) == cyf::CyfReader::OK) got.push_back(c);
    bool ok = (st == cyf::CyfReader::END) && got.size() == cells.size();
    for (size_t i = 0; ok && i < cells.size(); ++i) ok &= cellEq(cells[i], got[i]);
    check(ok, "cell round-trip + clean EOF");
  }

  // truncation: drop the 8-byte EOF marker
  {
    std::ostringstream os(std::ios::binary);
    { cyf::CyfWriter w(os); w.writeHeader(h); for (auto& c : cells) w.writeCell(c); }
    std::string b = os.str(); b.resize(b.size() - 8);
    std::istringstream is(b, std::ios::binary);
    cyf::CyfReader r(is); CellHeader h2; r.readHeader(h2);
    cyf::CyfReader::Status st; Cell c;
    while ((st = r.readCell(c)) == cyf::CyfReader::OK) {}
    check(st == cyf::CyfReader::TRUNCATED, "truncation detected (no EOF marker)");
  }

  std::cout << "adapter round-trips (CYF-binary / CYF-text):\n";
  struct FmtCase { cyf::OutFormat fmt; const char* name; };
  for (FmtCase fc : { FmtCase{cyf::OutFormat::Binary, "binary .byf (BGZF)"},
                      FmtCase{cyf::OutFormat::Text,   "text .cyf (tab)"} }) {
    std::ostringstream os(std::ios::binary);
    { OutArchive ar(os, fc.fmt); ar(h); for (auto& c : cells) ar(c); }
    std::istringstream is(os.str(), std::ios::binary);
    InArchive ar(is);
    CellHeader h2; ar(h2);
    std::vector<Cell> got; Cell c;
    while (ar.next(c)) got.push_back(c);
    bool ok = got.size() == cells.size();
    for (size_t i = 0; ok && i < cells.size(); ++i) ok &= cellEq(cells[i], got[i]);
    check(ok, (std::string("OutArchive(") + fc.name + ") -> InArchive auto-detects + round-trips").c_str());
  }

  // Legacy cereal (.ocyf) is never written by OutArchive, but InArchive must still
  // READ it so existing files can be migrated. Write a raw cereal stream directly.
  {
    std::ostringstream os(std::ios::binary);
    { cereal::PortableBinaryOutputArchive ar(os); ar(h); for (auto& c : cells) ar(c); }
    std::istringstream is(os.str(), std::ios::binary);
    InArchive ar(is);
    CellHeader h2; ar(h2);
    std::vector<Cell> got; Cell c;
    while (ar.next(c)) got.push_back(c);
    bool ok = got.size() == cells.size();
    for (size_t i = 0; ok && i < cells.size(); ++i) ok &= cellEq(cells[i], got[i]);
    check(ok, "legacy cereal still readable by InArchive (.ocyf migration path)");
  }

  // formatForPath: extension drives the writer; .ocyf is refused.
  {
    bool ok = cyf::formatForPath("a.cyf")  == cyf::OutFormat::Text  &&
              cyf::formatForPath("a.byf") == cyf::OutFormat::Binary &&
              cyf::formatForPath("-")      == cyf::OutFormat::Binary;
    bool threw = false;
    try { cyf::formatForPath("a.ocyf"); } catch (const std::exception&) { threw = true; }
    check(ok && threw, "formatForPath: .cyf/.byf/'-' map correctly and .ocyf throws");
  }

  std::cout << "typed columns (int / string / categorical / array):\n";
  {
    const CellHeader th = buildTypedHeader();

    // header schema parses correctly
    std::vector<CyfValueType> got = th.ColumnTypes();
    std::vector<CyfValueType> want = {CyfValueType::Float, CyfValueType::Int32,
                                      CyfValueType::String, CyfValueType::Category,
                                      CyfValueType::Array};
    check(got == want, "CellHeader::ColumnTypes parses TY: into the right schema");
    check(th[4].CategoryLevels().size() == 3 && th[4].CategoryLevels()[2] == "immune",
          "Tag::CategoryLevels parses LV:");

    std::vector<cyf::CyfRow> rows = {
      typedRow(1, 12.5, 42, "cellA", 2, {2, 3, 4}),
      typedRow(2, 0.0, -7, "", 0, {}),
    };
    std::ostringstream os(std::ios::binary);
    { cyf::CyfWriter w(os); w.writeHeader(th); for (auto& r : rows) w.writeRow(r); }
    std::istringstream is(os.str(), std::ios::binary);
    cyf::CyfReader r(is);
    CellHeader h2;
    check(r.readHeader(h2), "typed readHeader");
    check(r.columnTypes() == want, "on-disk type codes match schema");
    std::vector<cyf::CyfRow> back; cyf::CyfRow row;
    cyf::CyfReader::Status st;
    while ((st = r.readRow(row)) == cyf::CyfReader::OK) back.push_back(row);
    bool ok = (st == cyf::CyfReader::END) && back.size() == rows.size();
    for (size_t i = 0; ok && i < rows.size(); ++i) ok &= rowEq(rows[i], back[i]);
    check(ok, "typed row round-trip (string/category/array preserved)");
  }

  std::cout << "self-describing flags (@FL):\n";
  {
    CellHeader h; int i = 0;
    h.appendRawTag(mk(Tag::HD_TAG, "", "VN:1.0", i++));
    h.appendRawTag(mk(Tag::MA_TAG, "DAPI", "TY:f", i++));
    cyf::addStandardFlagTags(h);
    check(h.GetFlagTags().size() == 11, "addStandardFlagTags declares 11 cflag bits");
    cyf::addStandardFlagTags(h);                       // call again
    check(h.GetFlagTags().size() == 11, "addStandardFlagTags is idempotent");

    std::ostringstream os(std::ios::binary);
    { cyf::CyfWriter w(os); w.writeHeader(h); }
    std::istringstream is(os.str(), std::ios::binary);
    cyf::CyfReader r(is); CellHeader h2;
    check(r.readHeader(h2), "header with @FL reads back");
    auto fl = h2.GetFlagTags();
    bool found_tls = false;
    for (auto& t : fl)
      if (t.id == "TLS" && t.GetField("RG") == "cflag" && t.GetField("BI") == "5") found_tls = true;
    check(fl.size() == 11 && found_tls, "@FL round-trips through CYF (TLS = cflag bit 5)");
  }

  std::cout << (failures == 0 ? "\nOK: all CYF tests passed\n" : "\nFAILED\n");
  return failures == 0 ? 0 : 1;
}
