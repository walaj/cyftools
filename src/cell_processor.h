#pragma once

// DataProcessor.h
#include <string>
#include "cell_header2.h"
#include "cell_row.h"
#include "polygon.h"
#include "cysift.h"

#include <cereal/types/vector.hpp>
#include <cereal/archives/portable_binary.hpp>

class CellProcessor {
 public:
  virtual ~CellProcessor() = default;
  
  virtual int ProcessHeader(CellHeader& header) = 0;

  virtual int ProcessLine(Cell& cell) = 0;

  void SetOutput(const std::string& file) {

    // set the output to file or stdout
    if (file == "-") {
      m_archive = std::make_unique<cereal::PortableBinaryOutputArchive>(std::cout);
    } else {
      m_os = std::make_unique<std::ofstream>(file, std::ios::binary);
      m_archive = std::make_unique<cereal::PortableBinaryOutputArchive>(*m_os);
    }

    assert(m_archive);
  }
  
  void OutputLine(const Cell& cell) {
    assert(m_archive);
    (*m_archive)(cell);
  }
  
  // common archiver objects
  std::unique_ptr<std::ofstream> m_os;
  std::unique_ptr<cereal::PortableBinaryOutputArchive> m_archive;
  
};

class LineProcessor {
 public:
  virtual ~LineProcessor() = default;
  
  virtual int ProcessHeader(CellHeader& header) = 0;
  
  virtual int ProcessLine(const std::string& line) = 0;
};


// Cut processor
class CutProcessor : public CellProcessor {
  
 public:
  
  void SetParams(const std::unordered_set<std::string>& include,
		 bool header_print, bool strict_cut) {
    m_include = include;
    m_header_print = header_print;
    m_strict_cut = strict_cut;
  }
  
  int ProcessHeader(CellHeader& header) override;
  
  int ProcessLine(Cell& cell) override;
  
 private:
  
  std::unordered_set<std::string> m_include;

  bool m_header_print;

  bool m_strict_cut = false;
  
  std::unordered_set<size_t> to_remove;
};

// Pheno processor
class PhenoProcessor : public CellProcessor {
  
 public:
  
  void SetParams(PhenoMap p, const std::string& cmd,
		 const std::string& output_file) {
    m_p = p;
    m_cmd = cmd;
    m_output_file = output_file;
  }
  
  int ProcessHeader(CellHeader& header) override;
  
  int ProcessLine(Cell& cell) override;
  
 private:
  CellHeader m_header;
  std::unordered_map<std::string, size_t> m_marker_map;
  PhenoMap m_p;
  std::string m_cmd;
  std::string m_output_file;
};


// Count processor
class CountProcessor : public CellProcessor {
  
 public:
  
  void SetParams() {}
  
  int ProcessHeader(CellHeader& header) override;
  
  int ProcessLine(Cell& cell) override;

  void PrintCount();
  
 private:
  size_t m_count = 0;
};


// Select processor
class SelectProcessor : public CellProcessor {
  
 public:

  void SetParams(uint64_t logor, uint64_t logand, bool lognot,
		 const std::string& cmd)  {
    m_or = logor;
    m_and = logand;
    m_not = lognot;
    m_cmd = cmd;
  }
  
  int ProcessHeader(CellHeader& header) override;
  
  int ProcessLine(Cell& cell) override;
  
 private:

  uint64_t m_or;
  uint64_t m_and;
  bool m_not;

  CellHeader m_header;
  
  std::string m_cmd;
  
  size_t m_flag_index;
  
};


// Log10 processor
class LogProcessor : public CellProcessor {
  
 public:

  void SetParams()  { }
  
  int ProcessHeader(CellHeader& header) override;

  int ProcessLine(Cell& cell) override;
  
 private:

  size_t m_line_number = 0;
  
  std::unordered_set<size_t> m_to_log;
  
};

class ROIProcessor : public CellProcessor {

 public:
  
  void SetParams(bool label,
		 const std::vector<Polygon>& rois)  {
    m_label = label;
    m_rois = rois;
  }
  
  int ProcessHeader(CellHeader& header) override;

  int ProcessLine(Cell& cell) override;

 private:

  std::vector<Polygon> m_rois;

  CellHeader m_header;

  bool m_label;
  
  size_t m_col_count;

  size_t x_i;
  size_t y_i;
};

class ViewProcessor : public CellProcessor { 

 public:
  
  void SetParams(bool print_header,
		 bool header_only,
		 int round) {
    
    m_print_header = print_header;
    m_header_only = header_only;
    m_round = round;
    
  }
  
  int ProcessHeader(CellHeader& header) override;

  int ProcessLine(Cell& cell) override;
  
 private:

  bool m_header_only;
  bool m_print_header;
  int m_round;
};

class BuildProcessor : public CellProcessor { 

 public:
  
  void SetParams() {}
  
  int ProcessHeader(CellHeader& header) override;

  int ProcessLine(Cell& cell) override;
  
 private:

  CellHeader m_header;

};


class CatProcessor : public CellProcessor { 

 public:
  
  void SetParams(int offset,
		 int sample) {
    
    m_offset = offset;
    m_sample = sample;
    
  }

  void SetOffset(size_t offset) { m_offset = offset; }

  void SetSample(size_t sample) { m_sample = sample; }
  
  int ProcessHeader(CellHeader& header) override;

  int ProcessLine(Cell& cell) override;

  size_t GetMaxCellID() const { return m_max_cellid; }
  
private:
  
  size_t m_cellid_index = static_cast<size_t>(-1);
  
  size_t m_offset;
  size_t m_sample;
  
  // the header to compare with (and print if needed)
  CellHeader m_master_header;
  bool m_master_set = false;

  size_t m_max_cellid = 0;

  std::vector<size_t> m_graph_indicies;
  
};

class CerealProcessor : public LineProcessor { 

 public:
  
  void SetParams(const std::string& filename,
		 const std::string& cmd) {
    m_filename = filename;
    m_cmd = cmd;
  }
  
  int ProcessHeader(CellHeader& header) override;

  int ProcessLine(const std::string& line) override;
  //int ProcessLine(const std::string& line) override;

private:

  int cellid;
  std::vector<float> vec1;

  std::string m_cmd;
  std::string m_filename;
  
  CellHeader m_header;
  
  std::unique_ptr<std::ofstream> m_os;
  std::unique_ptr<cereal::PortableBinaryOutputArchive> m_archive;
  

};

class RadialProcessor : public CellProcessor { 

 public:

  void SetParams(const std::vector<uint64_t>& inner, const std::vector<uint64_t>& outer,
		 const std::vector<uint64_t>& logor, const std::vector<uint64_t>& logand,
		 const std::vector<std::string>& label, const std::string& cmd,
		 const std::string& output_file) {

    m_cmd = cmd;
    m_output_file = output_file;
    
    // check the radial geometry parameters
    assert(inner.size());
    assert(inner.size() == outer.size());
    assert(inner.size() == logor.size());
    assert(inner.size() == logand.size());
    assert(inner.size() == label.size());    

    m_inner = inner;
    m_outer = outer;
    m_logor = logor;
    m_logand = logand;
    m_label = label;
  }
  
  
  int ProcessHeader(CellHeader& header) override;

  int ProcessLine(Cell& cell) override;
  
 private:
  
  CellHeader m_header;

  std::vector<uint64_t> m_inner, m_outer, m_logor, m_logand;
  std::vector<std::string> m_label;

  std::string m_cmd;
  std::string m_output_file;
  
};
