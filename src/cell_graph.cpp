#include "cell_graph.h"

#include <charconv>
#include <string_view>

CellNode::CellNode(const Neighbors& neighbors,
		  const std::vector<uint64_t>& flags) {
  
  neighbors_ = neighbors;
  m_flags = flags;
  
}


CellNode::CellNode(const std::vector<uint32_t>& ids,
		   const std::vector<uint32_t>& dist,
		   const std::vector<uint64_t>& flags) {

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
