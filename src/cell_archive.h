#pragma once
//
// cell_archive.h — a thin compatibility layer that lets cyftools read and write
// EITHER the legacy cereal stream OR the new CYF format, without changing the
// (*archive)(header) / (*archive)(cell) call sites scattered through the code.
//
// Writing: OutArchive picks cereal or CYF once at construction (see
// cyf::useCyfOutput()) and forwards operator() to the chosen backend. A CYF
// stream's EOF marker is emitted when the embedded CyfWriter is destroyed, so
// the existing member-destruction order (archive before its ofstream) finalizes
// the file correctly with no new explicit "close" calls.
//
// Reading: InArchive auto-detects the format from the stream's magic bytes
// (cyf::detectCyf, which replays the peeked bytes so stdin pipes work) and
// presents a uniform interface: operator()(header) then next(cell) until false.
// Truncated CYF streams throw rather than silently ending.
//
// This is the single seam to remove when cereal is finally dropped: delete the
// cereal branches and OutArchive/InArchive collapse onto CyfWriter/CyfReader.
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
  // Format defaults to the process-global switches; pass explicitly to override.
  // When CYF is selected and BGZF is on, the CYF byte stream is written through a
  // BGZF compressor (the .bcyf / BAM analog).
  explicit OutArchive(std::ostream& os,
                      bool use_cyf  = cyf::useCyfOutput(),
                      bool use_bgzf = cyf::useBgzfOutput()) {
    if (use_cyf) {
      if (use_bgzf) {
        m_bgzf = cyf::makeBgzfOutput(os);                  // BGZF container over os
        m_cyf  = std::make_unique<cyf::CyfWriter>(*m_bgzf);
      } else {
        m_cyf = std::make_unique<cyf::CyfWriter>(os);
      }
    } else {
      m_cereal = std::make_unique<cereal::PortableBinaryOutputArchive>(os);
    }
  }

  void operator()(const CellHeader& h) {
    if (m_cyf) m_cyf->writeHeader(h);
    else       (*m_cereal)(h);
  }

  void operator()(const Cell& c) {
    if (m_cyf) m_cyf->writeCell(c);
    else       (*m_cereal)(c);
  }

private:
  // Declaration order matters for finalization: m_cyf is destroyed first (it flushes
  // the CYF EOF marker into the BGZF stream), then m_bgzf (it writes the last block
  // + the BGZF EOF marker into os).
  std::unique_ptr<std::ostream>                        m_bgzf;
  std::unique_ptr<cereal::PortableBinaryOutputArchive> m_cereal;
  std::unique_ptr<cyf::CyfWriter>                      m_cyf;
};

// ----------------------------------------------------------------- reading
class InArchive {
public:
  explicit InArchive(std::istream& src) {
    // first peel off BGZF compression if present (gzip magic), then detect CYF vs cereal
    cyf::openMaybeBgzfInput(src, m_outer);      // m_outer = decompressed (or passthrough) bytes
    m_is_cyf = cyf::detectCyf(*m_outer, m_stream);
    if (m_is_cyf)
      m_cyf = std::make_unique<cyf::CyfReader>(*m_stream);
    else
      m_cereal = std::make_unique<cereal::PortableBinaryInputArchive>(*m_stream);
  }

  // Read the header (must be called first).
  void operator()(CellHeader& h) {
    if (m_cyf) {
      if (!m_cyf->readHeader(h))
        throw std::runtime_error("CYF: failed to read header (bad magic or version)");
    } else {
      (*m_cereal)(h);
    }
  }

  // Read the next cell. Returns true if one was read, false at clean end of
  // stream. Throws on a truncated/corrupt CYF stream.
  bool next(Cell& c) {
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

  bool isCyf() const { return m_is_cyf; }

private:
  // m_outer (BGZF/passthrough) is declared before m_stream so it is destroyed after it.
  std::unique_ptr<std::istream>                       m_outer;
  std::unique_ptr<std::istream>                       m_stream;
  bool                                                m_is_cyf = false;
  std::unique_ptr<cereal::PortableBinaryInputArchive> m_cereal;
  std::unique_ptr<cyf::CyfReader>                      m_cyf;
};
