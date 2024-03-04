#include "cell_table.h"

#ifdef HAVE_MLPACK
#include <mlpack/methods/dbscan/dbscan.hpp>
#include <mlpack/methods/gmm/gmm.hpp>
#include <mlpack/core.hpp>
#include <mlpack/methods/kmeans/kmeans.hpp>
#include <armadillo>
#endif

void CellTable::clusterDBSCAN(float epsilon,
			      size_t min_size,
			      size_t min_cluster_size
			      ) {

  validate();
  
#ifdef HAVE_MLPACK

  std::vector<double> xvec(m_x_ptr->begin(),m_x_ptr->end());
  std::vector<double> yvec(m_y_ptr->begin(),m_y_ptr->end());
  
  // Convert std::vector to Armadillo row vectors.
  arma::rowvec x_row(CellCount());
  arma::rowvec y_row(CellCount());

  size_t count = 0;
  std::vector<size_t> idx(CellCount());
  for (size_t i = 0; i < m_cflag_ptr->size(); i++) {
    
    // only cluster on marked cells
    if (IS_FLAG_SET(m_cflag_ptr->at(i), MARK_FLAG)) {
      x_row(count) = m_x_ptr->at(i);
      y_row(count) = m_y_ptr->at(i);
      idx[count] = i;
      count++;
    }
  }

  // check that there are cells to cluster
  if (count == 0) {
    std::cerr << "cysift dbscan -- no cells to cluster. Did you run cysift filter -M first?" << std::endl;
    return;
  }
  
  x_row.resize(count);
  y_row.resize(count);
  idx.resize(count);
  
  // Combine the two row vectors into a matrix.
  arma::mat dataset(2, x_row.size());
  dataset.row(0) = x_row;
  dataset.row(1) = y_row;

  if (m_verbose)
    std::cerr << "...running dbscan (mlpack) on " << AddCommas(x_row.size()) <<
      " cells. Epsilon: " << epsilon << " Min points: " <<
      min_size << " min cluster size: " << min_cluster_size << std::endl;
  
  // Perform DBSCAN clustering.
  arma::Row<size_t> assignments(idx.size(), arma::fill::zeros); // to store cluster assignments
  
  // The parameters are: epsilon (radius of neighborhood), minimum points in neighborhood
  const bool batchmode = false;
  mlpack::DBSCAN<> dbscan(epsilon, min_size, batchmode);
  dbscan.Cluster(dataset, assignments);
  
  // table of counts per cluster
  std::unordered_map<size_t, size_t> cluster_hist;
  for (const auto& a : assignments) {
    cluster_hist[a]++;
  }
  
  // store the clusters data
  FloatColPtr fc = std::make_shared<FloatCol>();
  fc->resize(CellCount()); // resize zeros as well
  for (const auto& cc : fc->getData())
    assert(cc == 0);
  
  // make the assignments
  std::unordered_set<size_t> assignment_set;
  assert(idx.size() == assignments.size());
  for (size_t i = 0; i < idx.size(); i++) {

    // don't add if not a big cluster
    if (cluster_hist[assignments[i]] < min_cluster_size) {
      fc->SetValueAt(idx[i], 0);
    // otherwise add the cluster to m_table
    } else {

      // not in a cluster
      if (assignments[i] == SIZE_MAX) {
	fc->SetValueAt(idx[i], 0);
	
      } else {
	fc->SetValueAt(idx[i], assignments[i] + 1);	
	assignment_set.insert(assignments[i]); // just to keep track of number of assignments
      }
    }
    
  }
  
  if (m_verbose)
    std::cerr << "...produced " << AddCommas(assignment_set.size()) << " clusters" << std::endl;

  Tag ttag1(Tag::CA_TAG, "dbscan_cluster","");
  assert(fc->size() == CellCount());
  AddColumn(ttag1, fc);
  
#else  
  std::cerr << "...mlpack header library not available. Need to include this during compilation" << std::endl;
#endif  
}

void CellTable::GMM_EM() {

#ifdef HAVE_MLPACK
  arma::mat dataset;
  mlpack::data::Load("data.csv", dataset, true);

  // Initialize with the default arguments.
  mlpack::GMM gmm(dataset.n_rows, 3); // 3 is the number of Gaussians in the model. Adjust as necessary.
  
  // Train the model.
  gmm.Train(dataset);
#else
  std::cerr << "Warning: Unable to run GMM without linking / including MLPACK and armadillo in build" << std::endl;
#endif 
}
