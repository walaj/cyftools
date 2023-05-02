#ifndef CELL_GRAPH_H
#define CELL_GRAPH_H

#include "umappp/Umap.hpp"
#include "umappp/NeighborList.hpp"
#include "umappp/combine_neighbor_sets.hpp"
#include "umappp/find_ab.hpp"
#include "umappp/neighbor_similarities.hpp"
#include "umappp/optimize_layout.hpp"
#include "umappp/spectral_init.hpp"
#include "knncolle/knncolle.hpp"

//using Neighbors = std::vector<umappp::Neighbor<float>>;
using Neighbors    = std::vector<umappp::Neighbor<double>>;
using NeighborsInt = std::vector<std::pair<int, int>>;

class CellNode {
 public:
  CellNode() = default;
  
  explicit CellNode(const Neighbors& neighbors);
  
  explicit CellNode(const std::string& input_str) { parse_neighbors(input_str); }
  
  const Neighbors& get_neighbors() const { return neighbors_; }
  
  void set_neighbors(const Neighbors& neighbors) { neighbors_ = neighbors; }

  std::string toString(bool integerize) const;

  size_t size() const { return neighbors_.size(); }
  
  friend std::ostream& operator<<(std::ostream& os, const CellNode& cn);
  
 private:
  Neighbors neighbors_;

  void parse_neighbors(const std::string& input_str);

};

/*
class CellGraph {
  
 public:
  CellGraph() = default;
  explicit CellGraph(const std::vector<CellNode>& nodes) : nodes_(nodes) {}
  
  const std::vector<CellNode>& get_nodes() const { return nodes_; }
  void set_nodes(const std::vector<CellNode>& nodes) { nodes_ = nodes; }
  
  void write_to_file(const std::string& filename) const;
  static CellGraph read_from_file(const std::string& filename);

  void AddNode(const CellNode& node) { nodes_.emplace_back(node); }

  using const_iterator = std::vector<CellNode>::const_iterator;
  const_iterator begin() const { return nodes_.begin(); }
  const_iterator end() const { return nodes_.end(); }

  //std::shared_ptr<StringColumn> toStringColumn() const;
  
 private:
  std::vector<CellNode> nodes_;
  
};
*/

#endif
