#pragma once

#include <string>
#include <vector>
#include <cstdint>

#include "cell_header2.h"

#include <cereal/archives/binary.hpp>

class Cell {

 public:

  Cell() {}

  Cell(const std::string& row, const CellHeader& header);
  
  template <class Archive>
    void serialize(Archive & ar)
    {
      ar(m_id, m_flag, m_x, m_y, m_cols, m_spatial_ids,
	 m_spatial_dist); 
    }
  
  uint32_t m_id;
  uint64_t m_flag;
  
  float m_x;
  float m_y;

  // column data
  std::vector<float> m_cols;

  // spatial graph
  std::vector<uint32_t> m_spatial_ids;
  std::vector<uint32_t> m_spatial_dist;

};
