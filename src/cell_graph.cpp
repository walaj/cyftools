#include "cell_graph.h"

#include <charconv>
#include <string_view>

void CellNode::parse_neighbors(const std::string& input_str) {
  std::string_view input_view(input_str);
  Neighbors neighbors;

  size_t start = 0;
  size_t end = 0;

  while ((end = input_view.find(';', start)) != std::string_view::npos) {
    std::string_view pair_view = input_view.substr(start, end - start);
    size_t caret_pos = pair_view.find('^');

    int index;
    float weight;

    std::from_chars(pair_view.data(), pair_view.data() + caret_pos, index);
    weight = std::strtof(pair_view.data() + caret_pos + 1, nullptr);

    neighbors.emplace_back(index, weight);

    start = end + 1;
  }

  // Handle the last pair if there's no trailing ';'
  if (start < input_view.size()) {
    std::string_view pair_view = input_view.substr(start);
    size_t caret_pos = pair_view.find('^');

    int index;
    float weight;

    std::from_chars(pair_view.data(), pair_view.data() + caret_pos, index);
    weight = std::strtof(pair_view.data() + caret_pos + 1, nullptr);

    neighbors.emplace_back(index, weight);
  }

  neighbors_ = neighbors;
}


CellNode::CellNode(const Neighbors& neighbors) : neighbors_(neighbors)
{
  
  

}

std::ostream& operator<<(std::ostream& os, const CellNode& cn) {
  os << cn.toString(false);
  return os;
}

// void CellNode::parse_neighbors(const std::string& input_str) {

//   std::istringstream line_stream(input_str);
//   Neighbors neighbors;
  
//   std::string pair_str;
//   while (std::getline(line_stream, pair_str, ';')) {
//     std::istringstream pair_stream(pair_str);
//     std::string index_str, weight_str;
    
//     std::getline(pair_stream, index_str, '^');
//     std::getline(pair_stream, weight_str, '^');
    
//     int index = std::stoi(index_str);
//     float weight = std::stof(weight_str);
//     neighbors.emplace_back(index, weight);
//   }
  
//   neighbors_ = neighbors;
  
// }

/*void CellGraph::write_to_file(const std::string& filename) const {
    std::ofstream outfile(filename);

    if (!outfile.is_open()) {
        std::cerr << "Failed to open the file for writing: " << filename << std::endl;
        return;
    }

    for (const auto& node : nodes_) {
        const auto& neighbors = node.get_neighbors();
        auto neighbor_it = neighbors.begin();
        if (neighbor_it != neighbors.end()) {
            outfile << neighbor_it->first << '^' << neighbor_it->second;
            ++neighbor_it;
        }
        for (; neighbor_it != neighbors.end(); ++neighbor_it) {
            outfile << ';' << neighbor_it->first << '^' << neighbor_it->second;
        }
        outfile << std::endl;
    }

    outfile.close();
}
*/

std::string CellNode::toString(bool integerize) const {

  std::ostringstream oss;
  auto neighbor_it = neighbors_.begin();
  if (neighbor_it != neighbors_.end()) {
    oss << neighbor_it->first << '^';
    if (integerize) {
      oss << static_cast<int>(neighbor_it->second);
    } else {
      oss << neighbor_it->second;
    }
    ++neighbor_it;
  }
  
  if (neighbor_it != neighbors_.end())
    for (; neighbor_it != neighbors_.end(); ++neighbor_it) {
      oss << ';' << neighbor_it->first << '^';
      if (integerize) {
	oss << static_cast<int>(neighbor_it->second);
      } else {
	oss << neighbor_it->second;
      }
    }
  return oss.str();
  
}

/*CellGraph CellGraph::read_from_file(const std::string& filename) {
    std::vector<CellNode> nodes;
    std::ifstream infile(filename);

    if (!infile.is_open()) {
        std::cerr << "Failed to open the file for reading: " << filename << std::endl;
        return CellGraph(nodes);
    }

    std::string line;
    while (std::getline(infile, line)) {
        std::istringstream line_stream(line);
        Neighbors neighbors;

        std::string pair_str;
        while (std::getline(line_stream, pair_str, ';')) {
            std::istringstream pair_stream(pair_str);
            std::string index_str, weight_str;

            std::getline(pair_stream, index_str, '^');
            std::getline(pair_stream, weight_str, '^');

            int index = std::stoi(index_str);
            float weight = std::stof(weight_str);
            neighbors.emplace_back(index, weight);
        }

        nodes.emplace_back(CellNode(neighbors));
    }

    infile.close();
    return CellGraph(nodes);
}
*/


/*std::shared_ptr<StringColumn> CellGraph::toStringColumn() const {

  auto strc = std::make_shared<StringColumn>();
  
  //std::vector<std::string> string_vec;
  //string_vec.reserve(nodes_.size());
  
  for (const auto& node : nodes_) {
    strc->PushElem(node.toString());
  }
  
  return strc;
  }*/
