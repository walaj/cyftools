#pragma once

#include "cysift.h"
#include <string>
#include <vector>
#include <cstdint>

#include "cell_header.h"

class Cell {

 public:

  Cell() {}

  Cell(const std::string& row, const CellHeader& header,
       uint32_t cellid,
       uint32_t sampleid);

  void Print(int round) const;

  void set_cell_id(uint32_t id);

  void set_sample_id(uint32_t id);

  uint32_t get_cell_id() const;

  uint32_t get_sample_id() const;

  template <class Archive>
  void serialize(Archive & ar)
  {
    ar(m_id, m_cell_flag, m_pheno_flag, m_x, m_y, m_cols, m_spatial_ids,
       m_spatial_dist, m_spatial_flags); 
  }
  
  friend std::ostream& operator<<(std::ostream& os, const Cell& cell);
  
  //ADD
  // should make these private?
  
  uint64_t m_id;
  
  cy_uint m_pheno_flag;
  cy_uint m_cell_flag;
  
  float m_x;
  float m_y;

  // column data
  std::vector<float> m_cols;

  // spatial graph
  std::vector<uint32_t> m_spatial_ids;
  std::vector<uint32_t> m_spatial_dist;
  std::vector<cy_uint> m_spatial_flags;

};
