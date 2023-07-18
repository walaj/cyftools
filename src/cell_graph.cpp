#include "cell_graph.h"

#include <charconv>
#include <string_view>

CellNode::CellNode(const JNeighbors& neighbors,
		  const std::vector<cy_uint>& flags) {
  
  neighbors_ = neighbors;
  m_flags = flags;
  
}

void CellNode::sort_ascending_distance() {
  
  // Check if neighborhood and m_flags have the same size
  if (neighbors_.size() != m_flags.size()) {
    throw std::runtime_error("neighborhood and m_flags must have the same size");
  }

  // Create a vector of indices
  std::vector<size_t> indices(neighbors_.size());
  for (size_t i = 0; i != indices.size(); ++i) {
    indices[i] = i;
  }
  
  // Sort the indices based on the values in the neighbors_ vector
  std::sort(indices.begin(), indices.end(),
	    [this](size_t i1, size_t i2) {
	      return neighbors_[i1].second < neighbors_[i2].second;
	    }
	    );
  
  // Use the sorted indices to reorder the neighbors_ and m_flags vectors
  JNeighbors temp_neighbors_ = neighbors_;
  //std::vector<umappp::Neighbor<float>> temp_neighbors_ = neighbors_;
  std::vector<cy_uint> temp_m_flags = m_flags;
  for (size_t i = 0; i != indices.size(); ++i) {
    neighbors_[i] = temp_neighbors_[indices[i]];
    m_flags[i] = temp_m_flags[indices[i]];
  }
  
}

CellNode::CellNode(const std::vector<uint32_t>& ids,
		   const std::vector<uint32_t>& dist,
		   const std::vector<cy_uint>& flags) {

  assert(ids.size() == dist.size());
  assert(ids.size() == flags.size());

  m_flags = flags;
  
  for (size_t i = 0; i < ids.size(); i++) {
    neighbors_.push_back({ids.at(i), dist.at(i)});
    i++;
  }
  
}


CellNode::CellNode(const std::vector<uint32_t>& ids,
		   const std::vector<uint32_t>& dist) {

  assert(ids.size() == dist.size());

  for (size_t i = 0; i < ids.size(); i++) {
    neighbors_.push_back({ids.at(i), dist.at(i)});
    i++;
  }
}

std::ostream& operator<<(std::ostream& os, const CellNode& cn) {
  os << cn.toString();
  return os;
}

std::string CellNode::toString() const {

  std::ostringstream oss;
  auto neighbor_it = neighbors_.begin();
  if (neighbor_it != neighbors_.end()) {
    oss << neighbor_it->first << '^';
    oss << neighbor_it->second;
    ++neighbor_it;
  }
  
  if (neighbor_it != neighbors_.end()) {
    for (; neighbor_it != neighbors_.end(); ++neighbor_it) {
      oss << ';' << neighbor_it->first << '^';
      oss << neighbor_it->second;
    }
  }
  
  return oss.str();
  
}

void CellNode::OffsetNodes(size_t offset) {

  for (auto& n : neighbors_)
    n.first += offset;
  
}
