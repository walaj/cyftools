#pragma once

#include <Eigen/Core>
#include <ldaplusplus/LDABuilder.hpp>

class CellLDA {

 public:

  CellLDA(size_t nt, size_t nc);

  CellLDA() {}
  
  Eigen::MatrixXi X;
  Eigen::VectorXi y;

  size_t n_topics;
  size_t n_classes;
  size_t n_docs;
  
  void run();
  
};
