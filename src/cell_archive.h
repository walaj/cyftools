#pragma once
//
// cell_archive.h — a thin layer presenting a uniform (*archive)(header) /
// (*archive)(cell) interface over the CYF formats, so the call sites scattered
// through the code do not change.
//
// Writing: OutArchive writes either text (.cyf) or BGZF binary (.byf), chosen by
// cyf::formatForPath(output). It never writes the legacy cereal form. A binary
// stream's EOF marker is emitted when the embedded CyfWriter is destroyed, so the
// existing member-destruction order (archive before its ofstream) finalizes the
// file correctly with no new explicit "close" calls.
//
// Reading: InArchive auto-detects the format from the stream's leading bytes
// (text '@', BGZF magic, CYF magic, else legacy cereal), replaying peeked bytes
// so stdin pipes work, and presents operator()(header) then next(cell) until
// false. Truncated binary streams throw rather than silently ending. The cereal
// branch is read-only and exists only to migrate old (.ocyf) files.
//
#include <memory>
#include <stdexcept>

#include <cereal/archives/portable_binary.hpp>

#include "cyf_io.h"
#include "bgzf.h"
#include "cell_header.h"
#include "cell_row.h"

// ----------------------------------------------------------------- writing
class OutArchive {
public:
  // Text = .cyf (SAM analog); Binary = .byf (BAM analog, always BGZF-framed).
  // Defaults to Binary; callers pass cyf::formatForPath(output) to honor the
  // file extension.
  explicit OutArchive(std::ostream& os, cyf::OutFormat fmt = cyf::OutFormat::Binary) {
    switch (fmt) {
      case cyf::OutFormat::Text:
        m_text = std::make_unique<cyf::CyfTextWriter>(os);
        break;
      case cyf::OutFormat::Binary:
        m_bgzf = cyf::makeBgzfOutput(os);                   // binary is always BGZF
        m_cyf  = std::make_unique<cyf::CyfWriter>(*m_bgzf);
        break;
    }
  }

  void operator()(const CellHeader& h) {
    if (m_text) m_text->writeHeader(h);
    else        m_cyf->writeHeader(h);
  }

  void operator()(const Cell& c) {
    if (m_text) m_text->writeCell(c);
    else        m_cyf->writeCell(c);
  }

private:
  // Declaration order matters for finalization: m_cyf is destroyed first (it flushes
  // the CYF EOF marker into the BGZF stream), then m_bgzf (last block + BGZF EOF).
  std::unique_ptr<std::ostream>       m_bgzf;
  std::unique_ptr<cyf::CyfWriter>     m_cyf;
  std::unique_ptr<cyf::CyfTextWriter> m_text;
};

// ----------------------------------------------------------------- reading
class InArchive {
public:
  explicit InArchive(std::istream& src) {
    // detection order: text ('@') → BGZF (gzip magic, decompress) → CYF binary → cereal
    m_is_text = cyf::detectText(src, m_l1);
    if (m_is_text) {
      m_text = std::make_unique<cyf::CyfTextReader>(*m_l1);
      return;
    }
    cyf::openMaybeBgzfInput(*m_l1, m_l2);       // m_l2 = decompressed (or passthrough) bytes
    m_is_cyf = cyf::detectCyf(*m_l2, m_stream);
    if (m_is_cyf)
      m_cyf = std::make_unique<cyf::CyfReader>(*m_stream);
    else
      m_cereal = std::make_unique<cereal::PortableBinaryInputArchive>(*m_stream);
  }

  // Read the header (must be called first).
  void operator()(CellHeader& h) {
    if      (m_text) m_text->readHeader(h);
    else if (m_cyf) {
      if (!m_cyf->readHeader(h))
        throw std::runtime_error("CYF: failed to read header (bad magic or version)");
    } else {
      (*m_cereal)(h);
    }
  }

  // Read the next cell. Returns true if one was read, false at clean end of
  // stream. Throws on a truncated/corrupt stream.
  bool next(Cell& c) {
    if (m_text) {
      cyf::CyfTextReader::Status st = m_text->readCell(c);
      if (st == cyf::CyfTextReader::OK)  return true;
      if (st == cyf::CyfTextReader::END) return false;
      throw std::runtime_error("CYF text: malformed record");
    }
    if (m_cyf) {
      cyf::CyfReader::Status st = m_cyf->readCell(c);
      if (st == cyf::CyfReader::OK)  return true;
      if (st == cyf::CyfReader::END) return false;
      throw std::runtime_error("CYF: truncated or corrupt cell stream");
    }
    try {
      (*m_cereal)(c);
      return true;
    } catch (const cereal::Exception&) {
      return false;   // cereal signals EOF by exception
    }
  }

  bool isCyf()  const { return m_is_cyf; }
  bool isText() const { return m_is_text; }

private:
  // Layering, declared before the readers that reference them (so destroyed after):
  // m_l1 (text-detect over src) -> m_l2 (BGZF/passthrough) -> m_stream (CYF-detect).
  std::unique_ptr<std::istream>                       m_l1;
  std::unique_ptr<std::istream>                       m_l2;
  std::unique_ptr<std::istream>                       m_stream;
  bool                                                m_is_cyf  = false;
  bool                                                m_is_text = false;
  std::unique_ptr<cereal::PortableBinaryInputArchive> m_cereal;
  std::unique_ptr<cyf::CyfReader>                     m_cyf;
  std::unique_ptr<cyf::CyfTextReader>                 m_text;
};
