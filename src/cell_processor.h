#pragma once

#include "cell_header.h"
#include "cell_row.h"
#include "cysift.h"
#include "cell_selector.h"
#include "polygon.h"

#include <string>
#include <cassert>
#include <random>

#include <cereal/types/vector.hpp>
#include <cereal/archives/portable_binary.hpp>

// operation type for SelectProcessor
enum optype {
  GREATER_THAN,
  LESS_THAN,
  EQUAL_TO,
  GREATER_THAN_OR_EQUAL,
  LESS_THAN_OR_EQUAL
};

typedef std::pair<optype, float> SelectOp;
typedef std::vector<SelectOp> SelectOpVec;
typedef std::unordered_map<std::string, SelectOpVec> SelectOpMap;
typedef std::unordered_map<int, SelectOpVec> SelectOpNumMap;

const float dummy_float = 23.28323130082;

class CellProcessor {
 public:

  static const int NO_WRITE_CELL = 0;
  static const int WRITE_CELL = 1;
  static const int SAVE_CELL = 2;
  static const int SAVE_NODATA_CELL = 3;
  static const int SAVE_NODATA_NOGRAPH_CELL = 4;  
  static const int HEADER_NO_ACTION = 0;
  static const int ONLY_WRITE_HEADER = 1;
  static const int SAVE_HEADER = 2;
  
  virtual ~CellProcessor() = default;
  
  virtual int ProcessHeader(CellHeader& header) = 0;

  virtual int ProcessLine(Cell& cell) = 0;

  void SetupOutputStream() { 

    // set the output to file or stdout
    if (m_output_file == "-") {
      m_archive = std::make_unique<cereal::PortableBinaryOutputArchive>(std::cout);
    } else {
      m_os = std::make_unique<std::ofstream>(m_output_file, std::ios::binary);
      m_archive = std::make_unique<cereal::PortableBinaryOutputArchive>(*m_os);
    }
    assert(m_archive);
  }
  
  void OutputLine(const Cell& cell) const { 
    assert(m_archive);
    (*m_archive)(cell);
  }

  void SetCommonParams(const std::string& output_file,
		       const std::string& cmd,
		       bool verbose) {
    m_output_file = output_file;
    m_cmd = cmd;
    m_verbose = verbose;
  }
  
protected:

  // copy of the header to manipiulate and output
  CellHeader m_header;

  // name of the output file. empty if none
  std::string m_output_file;

  // string representation of the input command,
  // to track in PG tag in header
  std::string m_cmd;
  
  // common archiver objects
  std::unique_ptr<std::ofstream> m_os;
  std::unique_ptr<cereal::PortableBinaryOutputArchive> m_archive;

  // increase verbosity
  bool m_verbose = false;

  // count the lines as they come, for verbosity
  size_t m_count = 0; 
};

class LineProcessor {
 public:
  virtual ~LineProcessor() = default;
  
  virtual int ProcessHeader(CellHeader& header) = 0;
  
  virtual int ProcessLine(const std::string& line) = 0;
};

class ScatterProcessor : public CellProcessor {

public:
  
  void SetParams(int width, int height, int seed) {
    m_width = width;
    m_height = height;
    m_seed = seed;
  }

  int ProcessHeader(CellHeader& header) override;
  
  int ProcessLine(Cell& cell) override;

private:

  int m_width = 0;
  int m_height = 0;
  int m_seed = 42;
  
  std::mt19937 m_gen; // random number generator
  std::uniform_int_distribution<> m_wdistrib;
  std::uniform_int_distribution<> m_hdistrib;   
  
};

class HallucinateProcessor : public CellProcessor {

public:

  void SetParams(size_t n, int seed) {
    m_nphenotypes = n;
    m_seed = seed;
  }

  int ProcessHeader(CellHeader& header) override;
  
  int ProcessLine(Cell& cell) override;

private:

  size_t m_nphenotypes = 0;

  int m_seed = 42;
  
  std::unordered_set<size_t> m_to_remove;

  std::mt19937 m_gen; // random number generator
  std::uniform_int_distribution<> m_distrib;

};


class CellCountProcessor : public CellProcessor {

public:

  void SetParams() {}
  
  int ProcessHeader(CellHeader& header) override;
  
  int ProcessLine(Cell& cell) override;

  void EmitCell() const;
  
 private:

  std::vector<size_t> m_counts;
  
};

class HeadProcessor : public CellProcessor {

public:

  void SetParams(size_t n) {
    m_n = n;
  }

  int ProcessHeader(CellHeader& header) override;
  
  int ProcessLine(Cell& cell) override;

  
 private:
  
  size_t m_n = 10; // head limit
  size_t m_current_n = 0; // current count
  
};

// Cut processor
class CutProcessor : public CellProcessor {
  
 public:
  
  void SetParams(const std::unordered_set<std::string>& include) {
    m_include = include;
  }
  
  int ProcessHeader(CellHeader& header) override;
  
  int ProcessLine(Cell& cell) override;
  
 private:
  
  std::unordered_set<std::string> m_include;
  
  std::unordered_set<size_t> m_to_remove;

};

class CleanProcessor : public CellProcessor {

public:

  void SetParams(bool clean_graph, bool clean_meta, bool clean_marker) {
    
    m_clean_graph = clean_graph;
    m_clean_meta = clean_meta;
    m_clean_marker = clean_marker;
    
  }
  
  int ProcessHeader(CellHeader& header) override;
  
  int ProcessLine(Cell& cell) override;


private:

  bool m_clean_graph  = false; 
  bool m_clean_meta   = false; 
  bool m_clean_marker = false; 

  std::unordered_set<size_t> m_to_remove;
  
};

class AverageProcessor : public CellProcessor {

public:

  void SetParams() {}

  int ProcessHeader(CellHeader& header) override;
  
  int ProcessLine(Cell& cell) override;

  void EmitCell() const;
  
private:

  size_t n = 0;
  std::vector<double> sums;
  
};


// Pheno processor
class PhenoProcessor : public CellProcessor {
  
 public:
  
  void SetParams(PhenoMap p) {
    m_p = p;
  }
  
  int ProcessHeader(CellHeader& header) override;
  
  int ProcessLine(Cell& cell) override;
  
 private:

  // map connecting marker names to indicies
  std::unordered_map<std::string, size_t> m_marker_map;

  // map of all of the gates (string - pair<float,float>)
  PhenoMap m_p;

};

// Tumor processor
/*class TumorProcessor : public CellProcessor {
  
 public:
  
  void SetParams(int n, cy_uint flag, float frac) {
    m_n = n;
    m_flag = flag;
    m_frac = frac;
  }
  
  int ProcessHeader(CellHeader& header) override;
  
  int ProcessLine(Cell& cell) override;
  
 private:

  int m_n;
  cy_uint m_flag;
  float m_frac;
 
};
*/

class SummaryProcessor : public CellProcessor {

 public:
  
  void SetParams() {}
  
  int ProcessHeader(CellHeader& header) override;
  
  int ProcessLine(Cell& cell) override;

  void Print();
  
 private:
  std::vector<float> m_min;
  std::vector<float> m_max;
  std::vector<float> m_mean;  
  size_t m_count = 0;
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
  
  void SetFlagParams(CellSelector select) {
    m_select = select;
  }

  void SetFieldParams(const SelectOpMap criteria, bool or_toggle) {
    m_criteria = criteria;
    m_or_toggle = or_toggle;
  }
  
  int ProcessHeader(CellHeader& header) override;
  
  int ProcessLine(Cell& cell) override;
  
 private:

  CellSelector m_select;
  
  // field selectors
  SelectOpMap m_criteria;
  SelectOpNumMap m_criteria_int;  
  bool m_or_toggle = false;
};

// Divide processor
class DivideProcessor : public CellProcessor {

 public:

  void SetParams(const std::string& num,
		 const std::string& den,
		 float div)  {
    m_numer_string = num;
    m_denom_string = den;
    m_div_zero = div;
  }
  
  int ProcessHeader(CellHeader& header) override;

  int ProcessLine(Cell& cell) override;
  
 private:

  std::string m_numer_string; // name numerator
  std::string m_denom_string; // name denominator
  size_t m_numer = -1; // index of numerator in data column
  size_t m_denom = -1; // index of denominator in data column  
  float m_div_zero; // value to be given to divde-by-zero
  size_t m_count = 0; // for verbose
  size_t m_existing_column = -1; // if the column for the answer already exists, overwrite
};

// Log10 processor
class LogProcessor : public CellProcessor {
  
 public:

  void SetParams()  { }
  
  int ProcessHeader(CellHeader& header) override;

  int ProcessLine(Cell& cell) override;
  
 private:

  std::unordered_set<size_t> m_to_log;

  bool m_bool_warning_emitted = false;
  
};

class ROIProcessor : public CellProcessor {

 public:
  
  void SetParams(bool label,
		 const std::vector<Polygon>& rois) {
    m_label = label;
    m_rois = rois;
  }
  
  int ProcessHeader(CellHeader& header) override;

  int ProcessLine(Cell& cell) override;

 private:

  std::vector<Polygon> m_rois;

  bool m_label;

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

  // number of integers to round output to
  int m_round;
};

class BuildProcessor : public CellProcessor { 

 public:
  
  void SetParams() {}
  
  int ProcessHeader(CellHeader& header) override;

  int ProcessLine(Cell& cell) override;
  
 private:

};


class SubsampleProcessor : public CellProcessor {

public:

  void SetParams(float rate, int seed, std::unordered_set<uint32_t> samplenum) {
    m_rate = rate;
    m_seed = seed;
    m_samplenum = samplenum;
  }

  int ProcessHeader(CellHeader& header) override;

  int ProcessLine(Cell& cell) override;

private:

  std::unordered_set<uint32_t> m_samplenum;
  float m_rate = 0;

  size_t m_count = 0;
  size_t m_kept = 0;
  
  int m_seed = 42;
  std::mt19937 m_gen; // Standard mersenne_twister_engine seeded with rd()
  std::uniform_real_distribution<> m_dis; 
  
  
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

  bool m_error_emitted = false; // print error just once per file
  
  // the header to compare with (and print if needed)
  CellHeader m_master_header;
  bool m_master_set = false;

  size_t m_max_cellid = 0;

  std::vector<size_t> m_graph_indicies;
  
};

class CerealProcessor : public LineProcessor { 

 public:
  
  void SetParams(const std::string& filename,
		 const std::string& cmd,
		 uint32_t sampleid) {
    m_filename = filename;
    m_cmd = cmd;
    m_sampleid = sampleid;
  }
  
  int ProcessHeader(CellHeader& header) override;

  int ProcessLine(const std::string& line) override;
  //int ProcessLine(const std::string& line) override;

private:

  uint32_t m_cellid = 0;
  uint32_t m_sampleid = 0;  
  std::vector<float> vec1;

  std::string m_cmd;
  std::string m_filename;
  
  CellHeader m_header;
  
  std::unique_ptr<std::ofstream> m_os;
  std::unique_ptr<cereal::PortableBinaryOutputArchive> m_archive;
  

};

/*class RadialProcessor : public CellProcessor { 

 public:

  void SetParams(const std::vector<cy_uint>& inner, const std::vector<cy_uint>& outer,
		 const std::vector<cy_uint>& logor, const std::vector<cy_uint>& logand,
		 const std::vector<std::string>& label, bool normalize) {

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
    m_normalize = normalize;
  }
  
  
  int ProcessHeader(CellHeader& header) override;

  int ProcessLine(Cell& cell) override;
  
 private:
  
  std::vector<cy_uint> m_inner, m_outer, m_logor, m_logand;
  std::vector<std::string> m_label;
  bool m_normalize;
};
*/
