#pragma once
//
// cyf_io.h — reference reader/writer for the CYF on-disk format (v1).
//
// This is a self-contained, explicit, little-endian codec that implements
// docs/CYF_FORMAT.md. It deliberately does NOT depend on cereal: the goal is
// for this to become the format core (a future libcyf) that any language can
// reimplement from the spec. It currently shares the in-memory Cell/CellHeader
// types with the rest of cyftools (those still carry cereal serialize() methods
// for the legacy path); when cereal is removed, only those type definitions
// change — the logic in this file does not.
//
// Layout (see docs/CYF_FORMAT.md §7):
//   magic "CYF\x01" | version u16 | reserved u16 | l_header u32 |
//   header_text[l_header] | n_dcol u32 | dcol_types u8[n_dcol] |
//   { id u64, x f32, y f32, cflag u64, pflag u64, cols f32[n_dcol] }* |
//   EOF_MARKER (8 bytes)
//
#include <cstdint>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "cyf_field.h"
#include "cell_header.h"
#include "cell_row.h"

namespace cyf {

// A fully-typed record: mandatory fields + one CyfField per data column, encoded
// per the header schema (CellHeader::ColumnTypes). This is the typed counterpart
// of the legacy all-float Cell; both are supported during the migration.
struct CyfRow {
  uint64_t id = 0;
  float    x = 0, y = 0;
  uint64_t cflag = 0, pflag = 0;
  std::vector<CyfField> values;   // one per data column, in header order
};

// On-disk constants (docs/CYF_FORMAT.md Appendix A)
static constexpr unsigned char MAGIC[4]      = {'C', 'Y', 'F', 0x01};
static constexpr unsigned char EOF_MARKER[8] = {'C', 'Y', 'F', 0x01, 'E', 'O', 'F', 0x00};
static constexpr uint16_t      FORMAT_VERSION = 1;
// The EOF marker reinterpreted as a little-endian u64: a reserved/illegal cell id.
static constexpr uint64_t      EOF_SENTINEL_ID = 0x00464F4501465943ULL;

// Canonical text rendering of a header (the @-lines), used verbatim as the
// binary header block and as the text-format header. Inverse of parseHeaderText.
std::string  renderHeaderText(const CellHeader& h);
CellHeader   parseHeaderText(const std::string& text);

// The two writable on-disk forms. Binary is always BGZF-framed (the .byf / BAM
// analog); Text is the tab-delimited .cyf (SAM analog). The legacy cereal form is
// read-only (see InArchive) and is never written.
enum class OutFormat { Text, Binary };

// Select the output format from a path's extension:
//   .cyf  -> Text      .byf -> Binary
//   .ocyf -> throws    (legacy cereal is read-only; convert to .cyf or .byf)
//   anything else (incl. "-"/stdout) -> Binary (the default)
// Reading is always content-detected (text '@', BGZF magic, CYF magic, else
// cereal), so the input extension is advisory only.
OutFormat formatForPath(const std::string& path);

// ---------------------------------------------------------------- Writer
class CyfWriter {
public:
  explicit CyfWriter(std::ostream& os) : m_os(os) {}
  ~CyfWriter() { finish(); }

  // Write preamble + verbatim header text + dcol type array. Call once, first.
  void writeHeader(const CellHeader& h);

  // Write one legacy all-float record (the column schema must be all float).
  // Throws if cols.size() disagrees with the schema, or if id hits the sentinel.
  void writeCell(const Cell& c);

  // Write one fully-typed record; each value is encoded per the header schema.
  void writeRow(const CyfRow& r);

  // Append the EOF marker. Idempotent; also called by the destructor.
  void finish();

private:
  std::ostream&             m_os;
  std::vector<CyfValueType> m_types;            // column schema captured at writeHeader
  bool                      m_all_float = true; // fast/legacy path when every column is 'f'
  bool                      m_header_written = false;
  bool                      m_finished = false;
};

// ---------------------------------------------------------------- Reader
class CyfReader {
public:
  enum Status { OK, END, TRUNCATED, BADMAGIC };

  explicit CyfReader(std::istream& is) : m_is(is) {}

  // Read + validate preamble and header. Returns false on bad magic or an
  // unsupported major version (sets status() accordingly).
  bool readHeader(CellHeader& h);

  // Read the next record into an all-float Cell (numeric columns are converted
  // to float; throws on a String/Array column, which Cell cannot hold). Returns
  // OK, END at the clean EOF marker, or TRUNCATED on a truncated stream.
  Status readCell(Cell& c);

  // Read the next record into a fully-typed CyfRow (no information loss).
  Status readRow(CyfRow& r);

  uint16_t version() const { return m_version; }
  Status   status()  const { return m_status; }
  const std::vector<CyfValueType>& columnTypes() const { return m_types; }

private:
  // shared record-prefix reader: returns false at EOF marker, throws on truncation
  bool readRecordPrefix(CyfRow& r);

  std::istream&             m_is;
  std::vector<CyfValueType> m_types;
  bool                      m_all_float = true;
  uint16_t                  m_version = 0;
  bool                      m_header_read = false;
  Status                    m_status = OK;
};

// ---------------------------------------------------------------- Detection
//
// Peek the first 4 bytes of `src` to decide whether it is a CYF stream. Because
// `src` may be a non-seekable pipe (stdin), the consumed bytes are replayed: on
// return, *out is a fresh istream that yields the original byte sequence in full
// (the peeked prefix followed by the remainder of src). Use *out for all further
// reading regardless of the detected format. Returns true if CYF magic is found.
//
// The returned istream keeps a reference to `src`; `src` must outlive *out.
bool detectCyf(std::istream& src, std::unique_ptr<std::istream>& out);

// ---------------------------------------------------------------- text form
// The CYF text format (the SAM analog): the @-line header followed by one
// tab-delimited record per line — id, x, y, cflag, pflag, then one field per data
// column (header order). Round-trips with the binary form.

class CyfTextWriter {
public:
  explicit CyfTextWriter(std::ostream& os) : m_os(os) {}
  void writeHeader(const CellHeader& h);
  void writeCell(const Cell& c);
private:
  std::ostream& m_os;
  std::size_t   m_ndcol = 0;
};

class CyfTextReader {
public:
  enum Status { OK, END, BADFORMAT };
  explicit CyfTextReader(std::istream& is) : m_is(is) {}
  bool   readHeader(CellHeader& h);
  Status readCell(Cell& c);
private:
  std::istream& m_is;
  std::size_t   m_ndcol = 0;
  std::string   m_pending;       // first record line, read while scanning the header
  bool          m_has_pending = false;
};

// True if `src` begins with '@' (the CYF text header). Replays the peeked bytes
// into *out, like detectCyf.
bool detectText(std::istream& src, std::unique_ptr<std::istream>& out);

} // namespace cyf
