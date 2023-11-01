#pragma once

#include "cysift.h"
#include <string>
#include <vector>
#include <cstdint>

#include "cell_header.h"

class Cell {

 public:

  Cell() : id(0), pflag(0), cflag(0), x(0), y(0) {}

  Cell(const std::string& row, const CellHeader& header,
       uint32_t cellid,
       uint32_t sampleid);

  void Print(int round) const;

  void PrintWithHeader(int round, const CellHeader& header) const;  

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
