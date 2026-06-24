#include "cyf_io.h"
#include "bgzf.h"

#include <algorithm>
#include <charconv>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <vector>

#include <zlib.h>

namespace cyf {

OutFormat formatForPath(const std::string& path) {
  auto ends_with = [&](const char* suf) {
    const std::string s(suf);
    return path.size() >= s.size() &&
           path.compare(path.size() - s.size(), s.size(), s) == 0;
  };
  if (ends_with(".byf")) return OutFormat::Binary;
  if (ends_with(".cyf"))  return OutFormat::Text;
  if (ends_with(".ocyf"))
    throw std::runtime_error(
        "refusing to write legacy cereal format: '" + path +
        "' (.ocyf is read-only); convert to .cyf (text) or .byf (binary)");
  return OutFormat::Binary;   // default for stdout '-' and unknown extensions
}

// ============================================================ LE primitives
// All integers/floats are written little-endian regardless of host byte order
// (docs/CYF_FORMAT.md §7), so files are portable across architectures.

static void putu16(std::ostream& os, uint16_t v) {
  unsigned char b[2] = {(unsigned char)(v), (unsigned char)(v >> 8)};
  os.write((const char*)b, 2);
}
static void putu32(std::ostream& os, uint32_t v) {
  unsigned char b[4];
  for (int i = 0; i < 4; ++i) b[i] = (unsigned char)(v >> (8 * i));
  os.write((const char*)b, 4);
}
static void putu64(std::ostream& os, uint64_t v) {
  unsigned char b[8];
  for (int i = 0; i < 8; ++i) b[i] = (unsigned char)(v >> (8 * i));
  os.write((const char*)b, 8);
}
static void putf32(std::ostream& os, float f) {
  uint32_t u;
  std::memcpy(&u, &f, 4);
  putu32(os, u);
}

static bool readN(std::istream& is, unsigned char* buf, std::size_t n) {
  is.read((char*)buf, (std::streamsize)n);
  return (std::size_t)is.gcount() == n;
}
static uint16_t getu16(const unsigned char* b) { return (uint16_t)(b[0] | (b[1] << 8)); }
static uint32_t getu32(const unsigned char* b) {
  return (uint32_t)b[0] | ((uint32_t)b[1] << 8) | ((uint32_t)b[2] << 16) | ((uint32_t)b[3] << 24);
}
static uint64_t getu64(const unsigned char* b) {
  uint64_t v = 0;
  for (int i = 0; i < 8; ++i) v |= (uint64_t)b[i] << (8 * i);
  return v;
}
static float getf32(const unsigned char* b) {
  uint32_t u = getu32(b);
  float f;
  std::memcpy(&f, &u, 4);
  return f;
}

static void putf64(std::ostream& os, double v) {
  uint64_t u; std::memcpy(&u, &v, 8); putu64(os, u);
}
static double getf64(const unsigned char* b) {
  uint64_t u = getu64(b); double d; std::memcpy(&d, &u, 8); return d;
}

// ============================================================ typed values
// Encode/decode one CyfField per a column's declared type (docs/CYF_FORMAT.md §3).

static void encodeArray(std::ostream& os, const CyfField& f) {
  os.put(f.sub);
  const bool isint = cyfIsIntegerCode(f.sub);
  const uint32_t n = isint ? (uint32_t)f.ai.size() : (uint32_t)f.ad.size();
  putu32(os, n);
  for (uint32_t k = 0; k < n; ++k) {
    if (isint) {
      const int64_t v = f.ai[k];
      switch (f.sub) {
        case 'i': putu32(os, (uint32_t)(int32_t)v); break;
        case 'I': putu32(os, (uint32_t)v); break;
        case 'l': putu64(os, (uint64_t)v); break;
        default:  putu64(os, (uint64_t)v); break;   // 'L'
      }
    } else {
      if (f.sub == 'f') putf32(os, (float)f.ad[k]); else putf64(os, f.ad[k]);
    }
  }
}

static void encodeField(std::ostream& os, CyfValueType t, const CyfField& f) {
  switch (t) {
    case CyfValueType::Int32:    putu32(os, (uint32_t)(int32_t)f.asInt64()); break;
    case CyfValueType::UInt32:   putu32(os, (uint32_t)f.asInt64());          break;
    case CyfValueType::Int64:    putu64(os, (uint64_t)f.asInt64());          break;
    case CyfValueType::UInt64:   putu64(os, (uint64_t)f.asInt64());          break;
    case CyfValueType::Float:    putf32(os, (float)f.asDouble());            break;
    case CyfValueType::Double:   putf64(os, f.asDouble());                   break;
    case CyfValueType::Category: putu32(os, f.cat);                          break;
    case CyfValueType::String:
      putu32(os, (uint32_t)f.s.size());
      os.write(f.s.data(), (std::streamsize)f.s.size());
      break;
    case CyfValueType::Array:    encodeArray(os, f);                         break;
  }
}

static bool decodeArray(std::istream& is, CyfField& out) {
  unsigned char b[8];
  unsigned char sb[1];
  if (!readN(is, sb, 1)) return false;
  const char sub = (char)sb[0];
  if (!readN(is, b, 4)) return false;
  const uint32_t n = getu32(b);
  if (cyfIsIntegerCode(sub)) {
    std::vector<int64_t> v; v.reserve(n);
    for (uint32_t k = 0; k < n; ++k) {
      if (sub == 'l' || sub == 'L') { if (!readN(is, b, 8)) return false; v.push_back((int64_t)getu64(b)); }
      else { if (!readN(is, b, 4)) return false;
             v.push_back(sub == 'i' ? (int64_t)(int32_t)getu32(b) : (int64_t)getu32(b)); }
    }
    out = CyfField::IntArray(std::move(v), sub);
  } else {
    std::vector<double> v; v.reserve(n);
    for (uint32_t k = 0; k < n; ++k) {
      if (sub == 'f') { if (!readN(is, b, 4)) return false; v.push_back(getf32(b)); }
      else            { if (!readN(is, b, 8)) return false; v.push_back(getf64(b)); }
    }
    out = CyfField::RealArray(std::move(v), sub);
  }
  return true;
}

static bool decodeField(std::istream& is, CyfValueType t, CyfField& out) {
  unsigned char b[8];
  switch (t) {
    case CyfValueType::Int32:    if (!readN(is,b,4)) return false; out = CyfField::Int((int32_t)getu32(b), t); return true;
    case CyfValueType::UInt32:   if (!readN(is,b,4)) return false; out = CyfField::Int((int64_t)getu32(b), t); return true;
    case CyfValueType::Int64:    if (!readN(is,b,8)) return false; out = CyfField::Int((int64_t)getu64(b), t); return true;
    case CyfValueType::UInt64:   if (!readN(is,b,8)) return false; out = CyfField::Int((int64_t)getu64(b), t); return true;
    case CyfValueType::Float:    if (!readN(is,b,4)) return false; out = CyfField::Real(getf32(b), t); return true;
    case CyfValueType::Double:   if (!readN(is,b,8)) return false; out = CyfField::Real(getf64(b), t); return true;
    case CyfValueType::Category: if (!readN(is,b,4)) return false; out = CyfField::Cat(getu32(b));    return true;
    case CyfValueType::String: {
      if (!readN(is,b,4)) return false;
      const uint32_t n = getu32(b);
      std::string s(n, '\0');
      if (n && !readN(is, (unsigned char*)&s[0], n)) return false;
      out = CyfField::Str(std::move(s));
      return true;
    }
    case CyfValueType::Array: return decodeArray(is, out);
  }
  return false;
}

// ============================================================ header text
// renderHeaderText / parseHeaderText are exact inverses (modulo the derived
// sort index `i`). One @-line per tag; fields are tab-separated; ID, when
// present, is emitted first. The trailing `data` blob is preserved verbatim,
// including any internal tabs, by splitting/rejoining on tab.

static uint8_t tagTypeFromCode(const std::string& code) {
  if (code == "MA") return Tag::MA_TAG;
  if (code == "CA") return Tag::CA_TAG;
  if (code == "GA") return Tag::GA_TAG;
  if (code == "PG") return Tag::PG_TAG;
  if (code == "HD") return Tag::HD_TAG;
  if (code == "SA") return Tag::SA_TAG;
  if (code == "FL") return Tag::FL_TAG;
  if (code == "IM") return Tag::IM_TAG;
  if (code == "RO") return Tag::RO_TAG;
  throw std::runtime_error("CYF: unknown header tag class '" + code + "'");
}

std::string renderHeaderText(const CellHeader& h) {
  std::ostringstream ss;
  for (const auto& t : h) {
    ss << '@' << t.PrintType();
    if (!t.id.empty()) ss << '\t' << "ID:" << t.id;
    if (!t.data.empty()) ss << '\t' << t.data;
    ss << '\n';
  }
  return ss.str();
}

CellHeader parseHeaderText(const std::string& text) {
  CellHeader h;
  std::istringstream lines(text);
  std::string line;
  int order = 0;
  while (std::getline(lines, line)) {
    if (line.empty() || line[0] != '@') continue;

    // split on tab
    std::vector<std::string> tok;
    std::string cur;
    std::istringstream ls(line);
    while (std::getline(ls, cur, '\t')) tok.push_back(cur);
    if (tok.empty()) continue;

    uint8_t type = tagTypeFromCode(tok[0].substr(1));
    std::string id;
    std::string data;
    for (std::size_t i = 1; i < tok.size(); ++i) {
      if (tok[i].rfind("ID:", 0) == 0) {
        id = tok[i].substr(3);
      } else {
        if (!data.empty()) data += '\t';
        data += tok[i];
      }
    }
    Tag tag(type, id, data);
    tag.i = order++;            // sequential file order: keeps SortTags stable
    h.appendRawTag(tag);
  }
  return h;
}

// ============================================================ CyfWriter
// Reheader-friendly layout (reserved-u16 flag = 1): the header text is padded up
// to a capacity and a BGZF block boundary is flushed right after the header, so
// the cell records start on their own BGZF block(s). `reheader` can then rewrite
// the header (within the padded capacity) and copy the compressed record blocks
// verbatim, without re-encoding any cell. Fully backward-compatible: l_header
// covers the padding, the parser skips the non-@ filler, and old readers that
// ignore the reserved field read the file unchanged.
static const uint16_t REHEADER_FRIENDLY = 1;
static uint32_t headerCapacityFor(std::size_t hlen) {
  const uint32_t unit = 64u * 1024u;                 // 64 KiB granularity; storage is cheap
  uint32_t cap = unit;
  while (cap < hlen + 1024) cap += unit;             // keep ≥ ~1 KiB slack for growth
  return cap;
}

// Serialize the binary CYF header (reheader-friendly layout) into `os` and flush a
// BGZF block boundary after it. Shared by CyfWriter::writeHeader and fastReheaderByf
// so the two always emit byte-identical header framing.
static void writeBinaryHeader(std::ostream& os, const CellHeader& h) {
  const std::vector<CyfValueType> types = h.ColumnTypes();
  std::string htext = renderHeaderText(h);
  const uint32_t cap = headerCapacityFor(htext.size());
  if (htext.size() < cap) {                          // pad with a non-@ filler line
    htext.push_back('\n');
    htext.append(cap - htext.size(), ' ');
  }
  os.write((const char*)MAGIC, 4);
  putu16(os, FORMAT_VERSION);
  putu16(os, REHEADER_FRIENDLY);                     // reserved -> layout flag
  putu32(os, (uint32_t)htext.size());
  os.write(htext.data(), (std::streamsize)htext.size());
  putu32(os, (uint32_t)types.size());
  for (CyfValueType t : types) os.put(cyfTypeCode(t));
  os.flush();                                        // close the header's BGZF block(s)
}

void CyfWriter::writeHeader(const CellHeader& h) {
  if (m_header_written) throw std::runtime_error("CYF: header already written");

  // The column schema = the value type of each data column (MA+CA), in order.
  // For a legacy header (no TY: fields) this is all Float, so the bytes below are
  // identical to the old all-float output.
  m_types = h.ColumnTypes();
  m_all_float = true;
  for (CyfValueType t : m_types)
    if (t != CyfValueType::Float) { m_all_float = false; break; }

  writeBinaryHeader(m_os, h);
  m_header_written = true;
}

void CyfWriter::writeCell(const Cell& c) {
  if (!m_header_written) throw std::runtime_error("CYF: writeCell before writeHeader");
  if (!m_all_float)
    throw std::runtime_error("CYF: writeCell (all-float) used on a typed schema; use writeRow");
  if (c.cols.size() != m_types.size())
    throw std::runtime_error("CYF: cell has " + std::to_string(c.cols.size()) +
                             " data values but header declares " + std::to_string(m_types.size()));
  if (c.id == EOF_SENTINEL_ID)
    throw std::runtime_error("CYF: cell id collides with the reserved EOF sentinel");

  putu64(m_os, c.id);
  putf32(m_os, c.x);
  putf32(m_os, c.y);
  putu64(m_os, (uint64_t)c.cflag);                  // fixed 64-bit on disk
  putu64(m_os, (uint64_t)c.pflag);
  for (float v : c.cols) putf32(m_os, v);
}

void CyfWriter::writeRow(const CyfRow& r) {
  if (!m_header_written) throw std::runtime_error("CYF: writeRow before writeHeader");
  if (r.values.size() != m_types.size())
    throw std::runtime_error("CYF: row has " + std::to_string(r.values.size()) +
                             " values but header declares " + std::to_string(m_types.size()));
  if (r.id == EOF_SENTINEL_ID)
    throw std::runtime_error("CYF: row id collides with the reserved EOF sentinel");

  putu64(m_os, r.id);
  putf32(m_os, r.x);
  putf32(m_os, r.y);
  putu64(m_os, r.cflag);
  putu64(m_os, r.pflag);
  for (std::size_t k = 0; k < m_types.size(); ++k)
    encodeField(m_os, m_types[k], r.values[k]);
}

void CyfWriter::finish() {
  if (m_finished) return;
  if (m_header_written)
    m_os.write((const char*)EOF_MARKER, 8);
  m_finished = true;
}

// ============================================================ CyfReader
bool CyfReader::readHeader(CellHeader& h) {
  unsigned char pre[12];
  if (!readN(m_is, pre, 12)) { m_status = TRUNCATED; return false; }
  if (std::memcmp(pre, MAGIC, 4) != 0) { m_status = BADMAGIC; return false; }

  m_version = getu16(pre + 4);
  if (m_version != FORMAT_VERSION) { m_status = BADMAGIC; return false; }
  // pre+6 reserved (ignored)
  uint32_t lh = getu32(pre + 8);

  std::string htext(lh, '\0');
  if (lh && !readN(m_is, (unsigned char*)&htext[0], lh)) { m_status = TRUNCATED; return false; }

  unsigned char nb[4];
  if (!readN(m_is, nb, 4)) { m_status = TRUNCATED; return false; }
  uint32_t ncol = getu32(nb);

  std::string codes(ncol, '\0');
  if (ncol && !readN(m_is, (unsigned char*)&codes[0], ncol)) { m_status = TRUNCATED; return false; }

  h = parseHeaderText(htext);

  // self-check: header's declared data columns must match the type array length
  if (h.GetDataTags().size() != ncol) { m_status = BADMAGIC; return false; }

  try {
    m_types.clear();
    m_types.reserve(ncol);
    for (char ch : codes) m_types.push_back(cyfTypeFromCode(ch));
  } catch (const std::exception&) { m_status = BADMAGIC; return false; }

  m_all_float = true;
  for (CyfValueType t : m_types)
    if (t != CyfValueType::Float) { m_all_float = false; break; }

  m_header_read = true;
  m_status = OK;
  return true;
}

// Read the mandatory record prefix (id/x/y/cflag/pflag). Returns true if a
// record was started; false with m_status == END at the EOF marker, or
// TRUNCATED if the stream ends mid-record.
bool CyfReader::readRecordPrefix(CyfRow& r) {
  unsigned char head[8];
  m_is.read((char*)head, 8);
  std::streamsize got = m_is.gcount();
  if (got == 0) { m_status = TRUNCATED; return false; }                       // ended without marker
  if (got == 8 && std::memcmp(head, EOF_MARKER, 8) == 0) { m_status = END; return false; }
  if (got != 8) { m_status = TRUNCATED; return false; }

  unsigned char rest[24];
  if (!readN(m_is, rest, 24)) { m_status = TRUNCATED; return false; }

  r.id    = getu64(head);
  r.x     = getf32(rest + 0);
  r.y     = getf32(rest + 4);
  r.cflag = getu64(rest + 8);
  r.pflag = getu64(rest + 16);
  return true;
}

CyfReader::Status CyfReader::readCell(Cell& c) {
  CyfRow pre;
  if (!readRecordPrefix(pre)) return m_status;   // END or TRUNCATED
  c.id    = pre.id;
  c.x     = pre.x;
  c.y     = pre.y;
  c.cflag = (cy_uint)pre.cflag;
  c.pflag = (cy_uint)pre.pflag;

  c.cols.resize(m_types.size());
  if (m_all_float) {
    for (std::size_t i = 0; i < m_types.size(); ++i) {
      unsigned char fb[4];
      if (!readN(m_is, fb, 4)) { m_status = TRUNCATED; return TRUNCATED; }
      c.cols[i] = getf32(fb);
    }
  } else {
    // typed columns: decode each, convert numerics to float. asFloat() throws on
    // a String/Array column (not representable in the legacy all-float Cell).
    for (std::size_t i = 0; i < m_types.size(); ++i) {
      CyfField f;
      if (!decodeField(m_is, m_types[i], f)) { m_status = TRUNCATED; return TRUNCATED; }
      c.cols[i] = f.asFloat();
    }
  }

  m_status = OK;
  return OK;
}

CyfReader::Status CyfReader::readRow(CyfRow& r) {
  if (!readRecordPrefix(r)) return m_status;   // END or TRUNCATED
  r.values.clear();
  r.values.resize(m_types.size());
  for (std::size_t i = 0; i < m_types.size(); ++i)
    if (!decodeField(m_is, m_types[i], r.values[i])) { m_status = TRUNCATED; return TRUNCATED; }

  m_status = OK;
  return OK;
}

// ============================================================ detection
namespace {
// A streambuf that first replays a small prefix, then streams the rest of an
// underlying istream in chunks. Lets us peek the magic on a non-seekable pipe
// and still hand a complete, efficient stream to whichever reader we pick.
class PrefixStreambuf : public std::streambuf {
public:
  PrefixStreambuf(std::string prefix, std::istream& src)
      : m_prefix(std::move(prefix)), m_src(src) {}

protected:
  int_type underflow() override {
    if (gptr() < egptr()) return traits_type::to_int_type(*gptr());

    if (!m_prefix_done) {
      m_prefix_done = true;
      if (!m_prefix.empty()) {
        char* base = &m_prefix[0];
        setg(base, base, base + m_prefix.size());
        return traits_type::to_int_type(*gptr());
      }
    }
    m_src.read(m_buf, (std::streamsize)sizeof(m_buf));
    std::streamsize n = m_src.gcount();
    if (n <= 0) return traits_type::eof();
    setg(m_buf, m_buf, m_buf + n);
    return traits_type::to_int_type(*gptr());
  }

private:
  std::string   m_prefix;
  bool          m_prefix_done = false;
  std::istream& m_src;
  char          m_buf[65536];
};

// An istream that owns its streambuf.
class OwningIStream : public std::istream {
public:
  explicit OwningIStream(std::unique_ptr<std::streambuf> buf)
      : std::istream(buf.get()), m_buf(std::move(buf)) {}
private:
  std::unique_ptr<std::streambuf> m_buf;
};
} // namespace

bool detectCyf(std::istream& src, std::unique_ptr<std::istream>& out) {
  char pre[4];
  src.read(pre, 4);
  std::streamsize n = src.gcount();
  std::string prefix(pre, (std::size_t)n);

  bool is_cyf = (n == 4 && std::memcmp(pre, MAGIC, 4) == 0);

  auto buf = std::unique_ptr<std::streambuf>(new PrefixStreambuf(prefix, src));
  out.reset(new OwningIStream(std::move(buf)));
  return is_cyf;
}

bool detectText(std::istream& src, std::unique_ptr<std::istream>& out) {
  char pre[1];
  src.read(pre, 1);
  std::streamsize n = src.gcount();
  bool is_text = (n == 1 && pre[0] == '@');
  auto buf = std::unique_ptr<std::streambuf>(new PrefixStreambuf(std::string(pre, (std::size_t)n), src));
  out.reset(new OwningIStream(std::move(buf)));
  return is_text;
}

// ============================================================ text form
// The shortest decimal that round-trips a float (std::to_chars, C++17).
static std::string fstr(float v) {
  char buf[32];
  auto r = std::to_chars(buf, buf + sizeof(buf), v);
  return std::string(buf, r.ptr);
}
static void rstripCR(std::string& s) { if (!s.empty() && s.back() == '\r') s.pop_back(); }

void CyfTextWriter::writeHeader(const CellHeader& h) {
  m_ndcol = h.GetDataTags().size();
  m_os << renderHeaderText(h);          // @-lines, tab-delimited, @HD first
}

void CyfTextWriter::writeCell(const Cell& c) {
  m_os << c.id << '\t' << fstr(c.x) << '\t' << fstr(c.y) << '\t'
       << (uint64_t)c.cflag << '\t' << (uint64_t)c.pflag;
  for (float v : c.cols) m_os << '\t' << fstr(v);
  m_os << '\n';
}

bool CyfTextReader::readHeader(CellHeader& h) {
  std::string headerText, line;
  while (std::getline(m_is, line)) {
    rstripCR(line);
    if (!line.empty() && line[0] == '@') { headerText += line; headerText += '\n'; }
    else { m_pending = line; m_has_pending = true; break; }   // first data record
  }
  h = parseHeaderText(headerText);
  m_ndcol = h.GetDataTags().size();
  return true;
}

CyfTextReader::Status CyfTextReader::readCell(Cell& c) {
  std::string line;
  if (m_has_pending) { line = m_pending; m_has_pending = false; }
  else { if (!std::getline(m_is, line)) return END; rstripCR(line); }
  if (line.empty()) return END;

  std::vector<std::string> f;
  std::string cur;
  std::istringstream ls(line);
  while (std::getline(ls, cur, '\t')) f.push_back(cur);
  if (f.size() < 5) return BADFORMAT;

  try {
    c.id    = std::stoull(f[0]);
    c.x     = std::stof(f[1]);
    c.y     = std::stof(f[2]);
    c.cflag = (cy_uint)std::stoull(f[3]);
    c.pflag = (cy_uint)std::stoull(f[4]);
    c.cols.clear();
    c.cols.reserve(f.size() - 5);
    for (std::size_t i = 5; i < f.size(); ++i) c.cols.push_back(std::stof(f[i]));
  } catch (const std::exception&) { return BADFORMAT; }
  return OK;
}

// ============================================================ fast reheader
// Read only the header of a (possibly BGZF) byf/cyf file into `out`.
bool readByfHeader(const std::string& file, CellHeader& out) {
  std::ifstream f(file, std::ios::binary);
  if (!f) return false;
  std::unique_ptr<std::istream> dec;
  openMaybeBgzfInput(f, dec);                  // dec = decompressed (or pass-through) bytes
  if (!dec) return false;
  CyfReader r(*dec);
  return r.readHeader(out);
}

static uint32_t getu32s(const std::string& s, std::size_t off) {
  return getu32(reinterpret_cast<const unsigned char*>(s.data()) + off);
}

// Rewrite ONLY the header of `infile` into `outfile`, copying the compressed cell
// records verbatim — no cell is ever decoded. Requires the reheader-friendly layout
// (records block-aligned, reserved flag set) and an UNCHANGED column schema. Returns
// false (caller should fall back to a full streaming rewrite) for old/un-flagged
// files, a header that straddles a block boundary, or any schema change.
bool fastReheaderByf(const std::string& infile, const std::string& outfile,
                     const CellHeader& newHeader) {
  std::ifstream in(infile, std::ios::binary);
  if (!in) return false;

  // Walk BGZF blocks, decompressing the header block(s), until the cumulative
  // uncompressed size reaches the header region (12 + l_header + 4 + ncol bytes).
  std::string hbuf;
  uint64_t cum = 0, record_offset = 0;
  uint32_t lh = 0, ncol = 0; uint16_t reserved = 0;
  bool have_lh = false, found = false;
  for (;;) {
    unsigned char bh[18];
    in.read((char*)bh, 18);
    if (in.gcount() != 18) break;
    if (bh[0] != 0x1f || bh[1] != 0x8b) break;
    const unsigned blocklen = (unsigned)(bh[16] | (bh[17] << 8)) + 1;
    if (blocklen < 26) break;
    const unsigned clen = blocklen - 18 - 8;
    std::vector<unsigned char> cdata(clen);
    if (clen && in.read((char*)cdata.data(), clen).gcount() != (std::streamsize)clen) break;
    unsigned char tail[8];
    if (in.read((char*)tail, 8).gcount() != 8) break;
    const uint32_t isize = tail[4] | (tail[5] << 8) | (tail[6] << 16) | ((uint32_t)tail[7] << 24);
    std::string un(isize, '\0');
    if (isize) {
      z_stream zs; std::memset(&zs, 0, sizeof(zs));
      if (inflateInit2(&zs, -15) != Z_OK) return false;
      zs.next_in = cdata.data(); zs.avail_in = clen;
      zs.next_out = (Bytef*)&un[0]; zs.avail_out = isize;
      const int rc = inflate(&zs, Z_FINISH); inflateEnd(&zs);
      if (rc != Z_STREAM_END) return false;
    }
    hbuf += un; cum += isize;
    if (!have_lh && hbuf.size() >= 12) {
      if (std::memcmp(hbuf.data(), MAGIC, 4) != 0) return false;
      reserved = (unsigned char)hbuf[6] | ((unsigned char)hbuf[7] << 8);
      lh = getu32s(hbuf, 8);
      have_lh = true;
    }
    if (have_lh && hbuf.size() >= (std::size_t)12 + lh + 4) {
      ncol = getu32s(hbuf, 12 + lh);
      const uint64_t hregion = 12ull + lh + 4 + ncol;
      if (cum == hregion) { record_offset = (uint64_t)in.tellg(); found = true; break; }
      if (cum > hregion) return false;        // header straddles a block boundary
    }
  }
  if (!found || reserved != REHEADER_FRIENDLY) return false;

  // The column schema must be unchanged (records depend on it byte-for-byte).
  const std::vector<CyfValueType> nt = newHeader.ColumnTypes();
  if (nt.size() != ncol) return false;
  for (uint32_t i = 0; i < ncol; ++i)
    if ((unsigned char)cyfTypeCode(nt[i]) != (unsigned char)hbuf[12 + lh + 4 + i]) return false;

  // record-block region = [record_offset, filesize - 28)  (minus the BGZF EOF marker)
  in.clear(); in.seekg(0, std::ios::end);
  const uint64_t fsize = (uint64_t)in.tellg();
  if (fsize < record_offset + 28) return false;
  uint64_t left = fsize - 28 - record_offset;

  std::ofstream out(outfile, std::ios::binary);
  if (!out) return false;
  {
    BgzfOStreambuf zbuf(out);
    std::ostream zos(&zbuf);
    writeBinaryHeader(zos, newHeader);          // new header block(s) + flush -> out
    in.clear(); in.seekg((std::streamoff)record_offset);
    std::vector<char> buf(1 << 16);
    while (left > 0) {
      const std::streamsize chunk = (std::streamsize)std::min<uint64_t>(left, buf.size());
      in.read(buf.data(), chunk);
      const std::streamsize got = in.gcount();
      if (got <= 0) return false;
      out.write(buf.data(), got);
      left -= (uint64_t)got;
    }
  }   // ~BgzfOStreambuf appends the BGZF EOF marker after the copied record blocks
  out.flush();
  return out.good();
}

} // namespace cyf
