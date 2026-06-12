#include "bgzf.h"

#include <cstdint>
#include <cstring>
#include <stdexcept>

#include <zlib.h>

namespace cyf {

// Max uncompressed bytes per block. BGZF blocks must be <= 65536 total; leave
// headroom for the gzip framing + worst-case deflate expansion.
static const unsigned BGZF_BLOCK = 0xff00; // 65280

// The standard 28-byte BGZF end-of-file marker (an empty block).
static const unsigned char BGZF_EOF[28] = {
  0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00,
  0x42, 0x43, 0x02, 0x00, 0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00,
  0x00, 0x00, 0x00, 0x00
};

// ======================================================== output
BgzfOStreambuf::BgzfOStreambuf(std::ostream& dst, int level)
    : m_dst(dst), m_level(level) {
  m_in.reserve(BGZF_BLOCK);
}

BgzfOStreambuf::~BgzfOStreambuf() { finalize(); }

// Compress `len` bytes of `data` into one BGZF block written to m_dst.
void BgzfOStreambuf::writeBlock(const char* data, unsigned len) {
  std::vector<unsigned char> comp(len + 256);
  z_stream zs;
  std::memset(&zs, 0, sizeof(zs));
  if (deflateInit2(&zs, m_level, Z_DEFLATED, -15, 8, Z_DEFAULT_STRATEGY) != Z_OK)
    throw std::runtime_error("BGZF: deflateInit2 failed");
  zs.next_in   = reinterpret_cast<Bytef*>(const_cast<char*>(data));
  zs.avail_in  = len;
  zs.next_out  = comp.data();
  zs.avail_out = static_cast<uInt>(comp.size());
  if (deflate(&zs, Z_FINISH) != Z_STREAM_END) { deflateEnd(&zs); throw std::runtime_error("BGZF: deflate failed"); }
  const unsigned clen = static_cast<unsigned>(comp.size() - zs.avail_out);
  deflateEnd(&zs);

  const uint32_t crc = crc32(crc32(0L, nullptr, 0), reinterpret_cast<const Bytef*>(data), len);
  const unsigned blocksize = 18 + clen + 8;     // header(18) + deflate + crc(4) + isize(4)
  const uint16_t bsize = static_cast<uint16_t>(blocksize - 1);

  unsigned char hdr[18] = {
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff,
    0x06, 0x00, 0x42, 0x43, 0x02, 0x00, 0x00, 0x00
  };
  hdr[16] = bsize & 0xff;
  hdr[17] = (bsize >> 8) & 0xff;
  m_dst.write(reinterpret_cast<char*>(hdr), 18);
  m_dst.write(reinterpret_cast<char*>(comp.data()), clen);

  unsigned char tail[8];
  tail[0] = crc & 0xff;        tail[1] = (crc >> 8) & 0xff;
  tail[2] = (crc >> 16) & 0xff; tail[3] = (crc >> 24) & 0xff;
  tail[4] = len & 0xff;        tail[5] = (len >> 8) & 0xff;
  tail[6] = (len >> 16) & 0xff; tail[7] = (len >> 24) & 0xff;
  m_dst.write(reinterpret_cast<char*>(tail), 8);
}

std::streambuf::int_type BgzfOStreambuf::overflow(int_type c) {
  if (c != traits_type::eof()) m_in.push_back(static_cast<char>(c));
  if (m_in.size() >= BGZF_BLOCK) {
    writeBlock(m_in.data(), static_cast<unsigned>(m_in.size()));
    m_in.clear();
  }
  return c;
}

int BgzfOStreambuf::sync() {
  if (!m_in.empty()) {
    writeBlock(m_in.data(), static_cast<unsigned>(m_in.size()));
    m_in.clear();
  }
  return 0;
}

void BgzfOStreambuf::finalize() {
  if (m_finalized) return;
  m_finalized = true;
  if (!m_in.empty()) {
    writeBlock(m_in.data(), static_cast<unsigned>(m_in.size()));
    m_in.clear();
  }
  m_dst.write(reinterpret_cast<const char*>(BGZF_EOF), sizeof(BGZF_EOF));
  m_dst.flush();
}

// ======================================================== input
BgzfIStreambuf::BgzfIStreambuf(std::istream& src) : m_src(src) {}

std::streambuf::int_type BgzfIStreambuf::underflow() {
  if (gptr() < egptr()) return traits_type::to_int_type(*gptr());

  for (;;) {
    unsigned char hdr[18];
    m_src.read(reinterpret_cast<char*>(hdr), 18);
    if (m_src.gcount() != 18) return traits_type::eof();         // clean end
    if (hdr[0] != 0x1f || hdr[1] != 0x8b)
      throw std::runtime_error("BGZF: bad block magic");

    const unsigned bsize    = static_cast<unsigned>(hdr[16] | (hdr[17] << 8)); // total - 1
    const unsigned blocklen = bsize + 1;
    if (blocklen < 26) throw std::runtime_error("BGZF: short block");
    const unsigned clen = blocklen - 18 - 8;

    std::vector<unsigned char> cdata(clen);
    if (clen && m_src.read(reinterpret_cast<char*>(cdata.data()), clen).gcount() != (std::streamsize)clen)
      throw std::runtime_error("BGZF: truncated block");
    unsigned char tail[8];
    if (m_src.read(reinterpret_cast<char*>(tail), 8).gcount() != 8)
      throw std::runtime_error("BGZF: truncated block tail");
    const uint32_t isize = tail[4] | (tail[5] << 8) | (tail[6] << 16) | (static_cast<uint32_t>(tail[7]) << 24);

    if (isize == 0) continue;   // empty block (e.g. the EOF marker) — read the next one

    m_out.resize(isize);
    z_stream zs;
    std::memset(&zs, 0, sizeof(zs));
    if (inflateInit2(&zs, -15) != Z_OK) throw std::runtime_error("BGZF: inflateInit2 failed");
    zs.next_in   = cdata.data();
    zs.avail_in  = clen;
    zs.next_out  = reinterpret_cast<Bytef*>(m_out.data());
    zs.avail_out = isize;
    const int rc = inflate(&zs, Z_FINISH);
    inflateEnd(&zs);
    if (rc != Z_STREAM_END) throw std::runtime_error("BGZF: inflate failed");

    setg(m_out.data(), m_out.data(), m_out.data() + isize);
    return traits_type::to_int_type(*gptr());
  }
}

// ======================================================== owning wrappers
namespace {
// A streambuf that replays a small prefix, then streams the rest of an istream.
class PrefixBuf : public std::streambuf {
public:
  PrefixBuf(std::string prefix, std::istream& src) : m_prefix(std::move(prefix)), m_src(src) {}
protected:
  int_type underflow() override {
    if (gptr() < egptr()) return traits_type::to_int_type(*gptr());
    if (!m_done) {
      m_done = true;
      if (!m_prefix.empty()) {
        setg(&m_prefix[0], &m_prefix[0], &m_prefix[0] + m_prefix.size());
        return traits_type::to_int_type(*gptr());
      }
    }
    m_src.read(m_buf, sizeof(m_buf));
    std::streamsize n = m_src.gcount();
    if (n <= 0) return traits_type::eof();
    setg(m_buf, m_buf, m_buf + n);
    return traits_type::to_int_type(*gptr());
  }
private:
  std::string m_prefix; bool m_done = false; std::istream& m_src; char m_buf[65536];
};

class OwningIStream : public std::istream {
public:
  explicit OwningIStream(std::unique_ptr<std::streambuf> b) : std::istream(b.get()), m_b(std::move(b)) {}
private:
  std::unique_ptr<std::streambuf> m_b;
};
class OwningOStream : public std::ostream {
public:
  explicit OwningOStream(std::unique_ptr<std::streambuf> b) : std::ostream(b.get()), m_b(std::move(b)) {}
private:
  std::unique_ptr<std::streambuf> m_b;
};
} // namespace

bool openMaybeBgzfInput(std::istream& src, std::unique_ptr<std::istream>& out) {
  char pre[2] = {0, 0};
  src.read(pre, 2);
  const std::streamsize n = src.gcount();
  const bool is_gz = (n == 2 && (unsigned char)pre[0] == 0x1f && (unsigned char)pre[1] == 0x8b);

  // a stream that yields the peeked bytes + the remainder of src
  auto prefixed = std::unique_ptr<std::streambuf>(new PrefixBuf(std::string(pre, (size_t)n), src));
  auto prefixedStream = std::make_shared<OwningIStream>(std::move(prefixed));

  if (!is_gz) {
    out.reset(new OwningIStream(std::unique_ptr<std::streambuf>(
        new PrefixBuf(std::string(pre, (size_t)n), src))));
    return false;
  }
  // BGZF: wrap the prefixed stream in a decompressing streambuf. Keep the
  // prefixed stream alive by capturing it inside the istream's buffer owner.
  struct BgzfOwner : public std::istream {
    std::shared_ptr<std::istream> under;
    std::unique_ptr<BgzfIStreambuf> buf;
    BgzfOwner(std::shared_ptr<std::istream> u)
        : std::istream(nullptr), under(std::move(u)),
          buf(new BgzfIStreambuf(*under)) { rdbuf(buf.get()); }
  };
  out.reset(new BgzfOwner(prefixedStream));
  return true;
}

std::unique_ptr<std::ostream> makeBgzfOutput(std::ostream& dst, int level) {
  return std::unique_ptr<std::ostream>(
      new OwningOStream(std::unique_ptr<std::streambuf>(new BgzfOStreambuf(dst, level))));
}

} // namespace cyf
