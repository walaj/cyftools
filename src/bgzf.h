#pragma once
//
// bgzf.h — BGZF (Blocked GZip Format) stream layer, the SAM/BAM compression
// container, implemented over zlib. A BGZF file is a series of standard gzip
// blocks (each <= 64 KiB uncompressed) carrying an extra "BC" field with the
// block size, terminated by a 28-byte empty EOF block. Because each block is
// independent, BGZF supports random access via virtual offsets (block offset,
// in-block offset) — the basis for indexing — while staying gzip-compatible.
//
// This is exposed as std::streambuf wrappers so the CYF codec is unchanged: it
// just reads from / writes to a BGZF stream. Reads are gzip-magic auto-detected
// upstream, so plain and BGZF inputs both work.
//
#include <iostream>
#include <memory>
#include <streambuf>
#include <vector>

namespace cyf {

// True if the first two bytes of `src` are the gzip/BGZF magic (1f 8b). In all
// cases `out` is set to a fresh istream replaying the peeked bytes followed by the
// rest of `src`; when BGZF, *out transparently decompresses. Use *out for reading.
// (`src` must outlive *out.)
bool openMaybeBgzfInput(std::istream& src, std::unique_ptr<std::istream>& out);

// An ostream that BGZF-compresses everything written to it into `dst`. Destroy the
// returned stream (or call its rdbuf's finalize) before closing `dst` so the last
// block and the BGZF EOF marker are flushed.
std::unique_ptr<std::ostream> makeBgzfOutput(std::ostream& dst, int level = 6);

// ---- the streambufs (exposed for advanced use / testing) ----

class BgzfOStreambuf : public std::streambuf {
public:
  BgzfOStreambuf(std::ostream& dst, int level = 6);
  ~BgzfOStreambuf() override;
  void finalize();                 // flush remaining data + write the BGZF EOF block (idempotent)
protected:
  int_type overflow(int_type c) override;
  int      sync() override;
private:
  void writeBlock(const char* data, unsigned len);
  std::ostream&     m_dst;
  int               m_level;
  std::vector<char> m_in;          // accumulates uncompressed bytes up to the block size
  bool              m_finalized = false;
};

class BgzfIStreambuf : public std::streambuf {
public:
  explicit BgzfIStreambuf(std::istream& src);
protected:
  int_type underflow() override;
private:
  std::istream&     m_src;
  std::vector<char> m_out;         // current decompressed block
};

} // namespace cyf
