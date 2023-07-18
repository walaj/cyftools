#pragma once

#include "cell_column.h"
#include "cell_header.h"
#include "cell_processor.h"
#include "cysift.h"

#include <cereal/types/vector.hpp>
#include <cereal/archives/portable_binary.hpp>
#include <cereal/archives/json.hpp>

// forward declarations
class Polygon;

#ifdef HAVE_TIFFLIB
#include "tiff_writer.h"
#endif

#ifdef HAVE_LDAPLUSPLUS
#include <Eigen/Core>
#include <ldaplusplus/LDABuilder.hpp>
#include <ldaplusplus/LDA.hpp>
#include <ldaplusplus/LDABuilder.hpp>
#include <ldaplusplus/NumpyFormat.hpp>
#include <ldaplusplus/events/ProgressEvents.hpp>
#endif

#ifdef HAVE_KDTREE
#include "KDTree.hpp"

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



class CellTable {
  
public:

  //////
  // Constructors
  //////
  CellTable(const std::string& cmd, size_t threads = 1)
    : m_cmd(cmd), m_threads(threads) {}
  
  CellTable() {}

  //////
  // Getters/setters
  //////
  void setCmd(const std::string& cmd); 
  void setThreads(size_t threads) { m_threads = threads; }
  void setVerbose(bool verbose) { m_verbose = verbose; }
  const CellHeader& getHeader() const { return m_header; }

  //////
  // Table IO
  //////
  void BuildTable(const std::string& file);

  void StreamTableCSV(LineProcessor& proc, const std::string& file);

  int StreamTable(CellProcessor& proc, const std::string& file);
  
  void SetupOutputWriter(const std::string& file);

  void OutputTable();
  
  void HDF5Write(const std::string& file) const;


  
  //////
  // Basic table operations
  //////
  size_t size() const { return m_table.size(); }

  void sortxy(bool reverse);

  void sort(const std::string& field, bool reverse);
  
  void AddColumn(const Tag& tag, ColPtr value);

  bool HasColumn(const std::string& col) const;

  friend std::ostream& operator<<(std::ostream& os, const CellTable& table);

  bool ContainsColumn(const std::string& name) const;

  size_t CellCount() const;

  void PlotASCII(int width, int height) const;
  
  //////
  // Spatial ops
  //////
  void BuildKDTree();

  void UMAP(int num_neighbors);
  
  void UMAPPlot(const std::string& file, int width, int height) const;
  
  void KNN_spatial(int num_neighbors, int dist);  

  void Delaunay(const std::string& pdf_delaunay,
		const std::string& pdf_voronoi,
	        int limit);

  int RadialDensityKD(std::vector<cy_uint> inner, std::vector<cy_uint> outer,
		      std::vector<cy_uint> logor, std::vector<cy_uint> logand,
		      std::vector<std::string> label);
  
  void TumorCall(int num_neighbors, float frac,
		 cy_uint orflag, cy_uint andflag, cy_uint dist);
  

  //////
  // Numeric ops
  //////
  void Log10();

  void PrintPearson(bool csv, bool sort) const;

  //////
  // Image ops
  //////
#ifdef HAVE_TIFFLIB
  void Convolve(TiffWriter& otoif, int boxwidth, float microns_per_pixel);
#endif
  
  // ML ops
  void GMM_EM();
  
  //////
  // Subsetting / filtering ops
  //////
  void Subsample(int n, int s);

  void Crop(float xlo, float xhi, float ylo, float yhi);

  void SubsetROI(const std::vector<Polygon> &polygons);

  void Remove(const std::string& token);
  
  void Cut(const std::set<std::string>& tokens);

  void select(cy_uint on, cy_uint off);

  //////
  // Phenotyping ops
  //////
  PhenoMap phenoread(const std::string& filename) const;
  
  void phenotype(const std::unordered_map<std::string, std::pair<float,float>>& thresh);

  //////
  // Latent Dirichlet Allocation ops
  //////
  void LDA_load_model(const std::string& model_file);

  void LDA_write_model(const std::string& model_out) const;
  
  void LDA_create_model(const std::vector<std::string>& marker_cols,
			size_t n_topics,
			size_t n_iterations);

  void LDA_score_cells(const std::string& pdffile); 
  
 private:

  unordered_map<string, ColPtr> m_table;

  // the cysift command
  std::string m_cmd;

  // IO
  std::unique_ptr<std::ofstream> m_os;
  std::unique_ptr<cereal::PortableBinaryOutputArchive> m_archive;

  CellHeader m_header;

  LDAModel m_ldamodel;

#ifdef HAVE_KDTREE
  KDTree m_kdtree;
#endif
  
  // for verbose
  size_t m_count = 0;
  
  // params
  bool m_verbose = false;
  size_t m_threads = 1;
  
  // internal member functions
#ifndef DOXYGEN_SHOULD_SKIP_THIS

  void initialize_cols();
  
  void add_cell_to_table(const Cell& cell, bool nodata, bool nograph);

  void print_correlation_matrix(const std::vector<std::pair<std::string, const ColPtr>>& data,
				const std::vector<std::vector<float>>& correlation_matrix, bool sort) const;
  
#endif    
};

