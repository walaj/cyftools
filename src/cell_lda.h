#pragma once


#include <Eigen/Core>
#include <ldaplusplus/LDABuilder.hpp>
#include "cell_table.h"


class CellLDA {

 public:

  CellLDA(size_t nt, size_t nc);

  CellLDA() {}
  
  Eigen::MatrixXi X;
  Eigen::VectorXi y;

  size_t n_topics;
  size_t n_docs;
  size_t n_words;
  
  void run(CellTable& tab);
  
};
