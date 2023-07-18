#pragma once

#ifdef HAVE_LDAPLUSPLUS
#include <Eigen/Core>
#include <ldaplusplus/LDABuilder.hpp>
#include <ldaplusplus/LDA.hpp>
#include <ldaplusplus/LDABuilder.hpp>
#include <ldaplusplus/NumpyFormat.hpp>
#include <ldaplusplus/events/ProgressEvents.hpp>
#endif

struct LDAModel {

#ifdef HAVE_LDAPLUSPLUS
  
  std::vector<double> alpha;
  
  std::vector<std::vector<double>> beta;

  std::vector<std::string> markers_to_run;

  std::string cmd;

  // Default constructor
  LDAModel() {}
  
  LDAModel(const Eigen::VectorXd& alpha_eigen,
	   const Eigen::MatrixXd& beta_eigen,
	   const std::vector<std::string>& markers);

  size_t getNumTopics() const { return beta.size(); } 
  
  // Method to output a ModelParameters object
  template<typename Scalar = double>
  ldaplusplus::parameters::ModelParameters<Scalar> toModelParameters() const {
    // First, convert the alpha vector to an Eigen vector
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> alpha_eigen(alpha.size());
    for (size_t i = 0; i < alpha.size(); ++i) {
      alpha_eigen(i) = static_cast<Scalar>(alpha[i]);
    }
    
    // Next, convert the beta vector of vectors to an Eigen matrix
    if (beta.empty()) {
      throw std::runtime_error("Beta cannot be empty");
    }
    Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> beta_eigen(beta.size(), beta[0].size());
    for (size_t i = 0; i < beta.size(); ++i) {
      if (beta[i].size() != beta[0].size()) {
	throw std::runtime_error("All inner vectors in beta must have the same size");
      }
      for (size_t j = 0; j < beta[i].size(); ++j) {
	beta_eigen(i, j) = static_cast<Scalar>(beta[i][j]);
      }
    }
    
    // Finally, construct a ModelParameters object from the Eigen vector and matrix
    ldaplusplus::parameters::ModelParameters<Scalar> params(std::move(alpha_eigen), std::move(beta_eigen));
    
    return params;
  }
  
  template<class Archive>
  void serialize(Archive & ar) {
     ar(cereal::make_nvp("alpha", alpha), 
	cereal::make_nvp("beta", beta), 
	cereal::make_nvp("markers", markers_to_run),
	cereal::make_nvp("cmd", cmd));
  }

  // make sure its initialized
  bool initialized() const {
    return alpha.size() > 0;
  }

  // json archive
  void JSONArchive(const std::string& output) const {
    std::ofstream os(output);
    cereal::JSONOutputArchive oarchive(os);
    oarchive(cereal::make_nvp("LDAModel", *this));    
    return;
  }
  
#endif  
};


