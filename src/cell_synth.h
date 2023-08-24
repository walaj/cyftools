#pragma once

#include <cstddef>
#include <string> 


class CellTable;

class CellSynth {

 public:

  CellSynth(size_t w, size_t h);

  void WriteTable(const std::string& outfile);

  void Clusters(size_t num_clusters,
		size_t points_per_cluster,
		double sigma,
		size_t num_markers);

  void SetSeed(int seed) { m_seed = seed; }
  
 private:

  int m_seed = 42;
  
  size_t m_width = 1000;
  size_t m_height = 1000;
   
  CellTable * m_table = nullptr;

  
};
