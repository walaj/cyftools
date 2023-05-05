#pragma once

// DataProcessor.h
#include <string>
#include "cell_header.h"
#include "polygon.h"

class CellProcessor {
 public:
  virtual ~CellProcessor() = default;
  
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
  
  int ProcessLine(const std::string& line) override;
  
 private:
  
  std::unordered_set<std::string> m_include;

  bool m_header_print;

  bool m_strict_cut = false;
  
  std::unordered_set<size_t> to_remove;
};

// Select processor
class SelectProcessor : public CellProcessor {
  
 public:

  void SetParams(uint64_t logor, uint64_t logand, bool lognot, bool print_header)  {
    m_or = logor;
    m_and = logand;
    m_not = lognot;
    m_print_header = print_header;
  }
  
  int ProcessHeader(CellHeader& header) override;
  
  int ProcessLine(const std::string& line) override;
  
 private:

  uint64_t m_or;
  uint64_t m_and;
  bool m_not;
  bool m_print_header;

  size_t m_flag_index;
  
};


// Log10 processor
class LogProcessor : public CellProcessor {
  
 public:

  void SetParams(bool print_header)  {
    m_print_header = print_header;
  }
  
  int ProcessHeader(CellHeader& header) override;
  
  int ProcessLine(const std::string& line) override;
  
 private:

  size_t m_line_number = 0;
  
  bool m_print_header;

  std::unordered_set<size_t> m_to_log;
  
};

class ROIProcessor : public CellProcessor {

 public:
  
  void SetParams(bool print_header,
		 bool label,
		 const std::vector<Polygon>& rois)  {
    m_print_header = print_header;
    m_label = label;
    m_rois = rois;
  }
  
  int ProcessHeader(CellHeader& header) override;
  
  int ProcessLine(const std::string& line) override;

 private:

  std::vector<Polygon> m_rois;

  CellHeader m_header;

  bool m_label;
  
  bool m_print_header;
  
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
  
  int ProcessLine(const std::string& line) override;
  
 private:

  bool m_header_only;
  bool m_print_header;
  int m_round;
  
};
