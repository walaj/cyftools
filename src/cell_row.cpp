#include "cell_row.h"

#include "cell_utils.h"

Cell::Cell(const std::string& row, const CellHeader& header) {

  /*
  //
  const std::vector<std::string> tokens = tokenize_comma_delimited(row);

  size_t i = 0;
  for (const auto& t : header.GetColTags()) {
    if (t.isIDTag()) {
      m_cellid = std::stol(tokens.at(i));
    } else if (t.isFlagTag()) {
      m_flag = std::stol(tokens.at(i));
    } else if (t.isXDim()) {
      m_x = std::strtof(tokens.at(i).c_str(), nullptr);
    } else if (t.isYDim()) {
      m_y = std::strtof(tokens.at(i).c_str(), nullptr);
    } else if (t.isGraphTag()) {
      CellNode node(tokens.at(i));
      node.FillSparseFormat(m_spatial_ids, m_spatial_dist);
    } else if (t.GetName() == "sample") {
      m_sample = std::stoi(tokens.at(i));
    } else if (t.isColumnTag()) {
      m_cols.emplace_back(std::strtof(tokens.at(i).c_str(), nullptr));
    } else {
      std::cerr << " NOT FOUND " << i << " -- " << tokens.at(i) << std::endl;
    }
    i++;
  }
  */
}
