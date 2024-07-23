#pragma once

#include "cysift.h"
#include <string>
#include <vector>
#include <cstdint>
#include <iomanip>

#include "cell_header.h"

class Cell {

 public:

  Cell() : id(0), pflag(0), cflag(0), x(0), y(0) {}

  Cell(const std::string& row,
       int id_index, 
       int x_index,// which column is X and Y
       int y_index,
       const CellHeader& header,
       uint32_t cellid,
       uint32_t sampleid);

  template <typename T>
  void outputValue(const std::string& prefix, const T& value, bool tabPrint, int width, char delimiter) const {
    std::cout << prefix;
    if (tabPrint) {
      std::cout << std::setw(width) << value;
    } else {
      std::cout << value;
    }
    std::cout << delimiter;
  }
  
  void Print(int round, bool tabprint) const;

  void PrintWithHeader(int round,
		       bool tabprint,
		       bool header_print,
		       const CellHeader& header,
          	       bool no_print_cellid_etc) const;

  void PrintForCrevasse(const CellHeader& header) const;

  void set_cell_id(uint32_t id);

  void set_sample_id(uint32_t id);

  uint32_t get_cell_id() const;

  uint32_t get_sample_id() const;

  template <class Archive>
  void serialize(Archive & ar)
  {
    ar(id, cflag, pflag, x, y, cols); 
  }
  
  friend std::ostream& operator<<(std::ostream& os, const Cell& cell);

  uint64_t id;
  
  cy_uint pflag;
  cy_uint cflag;
  
  float x;
  float y;

  // column data
  std::vector<float> cols;

  // graph
  // std::vector<uint32_t> spatial_ids;
  // std::vector<uint32_t> spatial_dist;
  // std::vector<cy_uint>  spatial_flags;

};
