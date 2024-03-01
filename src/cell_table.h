#pragma once

#include <set>

#include "cell_column.h"
#include "cell_header.h"
#include "cell_processor.h"
#include "cysift.h"
#include "cell_synth.h"

#ifdef HAVE_KNNCOLLE
#include "knncolle/knncolle.hpp"
#endif

// to-do - back this out, shouldn't have this here
#include "cell_processor.h"

#include <cereal/types/vector.hpp>
#include <cereal/archives/portable_binary.hpp>
#include <cereal/archives/json.hpp>

#ifdef HAVE_MLPACK
#include <mlpack/core.hpp>
#include <mlpack/methods/range_search/range_search.hpp>
#endif

// forward declarations
class Polygon;
class LDAModel;
class TiffWriter;
class CellSelector;

class CellTable {
  
public:

  //////
  // Constructors
  //////
  CellTable(const std::string& cmd, size_t threads = 1)
    : m_cmd(cmd), m_threads(threads) {}
  
  CellTable() {}

  // new table, with zerod data
  CellTable(size_t num_cells);
  
  //////
  // Getters/setters
  //////
  void setCmd(const std::string& cmd); 
  void setThreads(size_t threads) { m_threads = threads;}
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
  // Plotting
  //////
  int PlotPNG(const std::string& file,
	      float scale_factor,
	      const std::string& module,
	      const std::string& roifile,
	      const std::string& title	      
	      ) const;
  
  //////
  // Basic table operations
  //////
  void validate() const;
  
  size_t size() const { return m_table.size(); }

  void sortxy(bool reverse);

  void sort(const std::string& field, bool reverse);
  
  void AddColumn(const Tag& tag, FloatColPtr value);

  void DeleteColumn(const std::string& col);
  
  void AddICPXYColumn(ColPtr value,
		      const std::string& type);
  
  
  bool HasColumn(const std::string& col) const;

  friend std::ostream& operator<<(std::ostream& os, const CellTable& table);

  size_t CellCount() const;

  void PlotASCII(int width, int height) const;

  void Collapse();
  
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
		      std::vector<std::string> label, std::vector<int> normalize_local,
		      std::vector<int> normalize_global);
  
  void AnnotateCall(int num_neighbors, float frac, cy_uint dist, cy_uint flag_to_set);
  
  void TumorMargin(float dist, cy_uint tumor_flag, cy_uint margin_flag);

  void IslandFill(size_t n, int flag_from, bool invert_from,
		  int flag_to, bool invert_to);
  
  void MoranI(const std::vector<cy_uint>& flags);

  //////
  // Clustering ops
  //////
  void clusterDBSCAN(float epsilon,
		     size_t min_size,
		     size_t min_cluster_size
		     ); 
  
  //////
  // Numeric ops
  //////
  void PrintPearson(bool csv, bool sort) const;

  void PrintJaccardSimilarity(bool csv, bool sort, bool subset_score) const;
  
  //////
  // Image ops
  //////
#ifdef HAVE_TIFFLIB
  void Convolve(TiffWriter* otoif, int boxwidth, float microns_per_pixel);
#endif
  
  // ML ops
  void GMM_EM();

  //////
  // Null model ops
  //////
  void ScramblePflag(int seed, bool lock_flags, bool phenotype_only);
  
  //////
  // Subsetting / filtering ops
  //////
  void Subsample(int n, int s);

  void Crop(float xlo, float xhi, float ylo, float yhi);

  void SubsetROI(const std::vector<Polygon> &polygons);

  void Remove(const std::string& token);
  
  void Cut(const std::set<std::string>& tokens);

  void Select(CellSelector select, 
	      SelectOpMap criteria,
	      bool or_toggle,
	      float radius);

  void ClearFlag(const int flag);
  
  void FlagToFlag(const bool clear_flag_to,
		  const int flag_from, bool flag_from_negative, 
		  const int flag_to, bool flag_to_negative);
  
  //////
  // Phenotyping ops
  //////
  PhenoMap phenoread(const std::string& filename) const;
  
  void phenotype(const std::unordered_map<std::string, std::pair<float,float>>& thresh);

  //////
  // Latent Dirichlet Allocation ops
  //////
#ifdef HAVE_LDAPLUSPLUS  
  void LDA_load_model(const std::string& model_file);

  void LDA_write_model(const std::string& model_out) const;
  
  void LDA_create_model(const std::vector<std::string>& marker_cols,
			size_t n_topics,
			size_t n_iterations,
			int seed);

  void LDA_score_cells(const std::string& pdffile,
		       int topic_highlight,
		       float cont_cutoff); 

  //////
  // Kristin neighborhood analysis
  /////
  const std::vector<std::vector<double>> create_inverse_distance_weights();
  
#endif

 private:

  unordered_map<string, FloatColPtr> m_table;

  FloatColPtr m_x_ptr;
  FloatColPtr m_y_ptr;
  IntColPtr m_pflag_ptr;
  IntColPtr m_cflag_ptr;
  IDColPtr m_id_ptr;
  
  // the cysift command
  std::string m_cmd;

  // IO
  std::unique_ptr<std::ofstream> m_os;
  std::unique_ptr<cereal::PortableBinaryOutputArchive> m_archive;

  CellHeader m_header;

  LDAModel* m_ldamodel = nullptr;

#ifdef HAVE_MLPACK
  mlpack::RangeSearch<mlpack::EuclideanDistance> * ml_kdtree = nullptr;
#endif


  
  // for verbose
  size_t m_count = 0;
  
  // params
  bool m_verbose = false;
  size_t m_threads = 1;

  // neighbors structure

  
  // selected write
  std::unordered_set<int> m_cells_to_write;
  
  // internal member functions
#ifndef DOXYGEN_SHOULD_SKIP_THIS

  void initialize_cols();
  
  void add_cell_to_table(const Cell& cell, bool nodata, bool nograph);

  void print_correlation_matrix(const std::vector<std::string>& labels,
				const std::vector<std::vector<float>>& correlation_matrix, bool sort) const;
    
  int dfs(int cellIndex,
	  std::vector<int>& component_label,
	  int currentLabel,
	  const std::vector<std::vector<size_t>>& neighbors) const;

  knncolle::VpTree<knncolle::distances::Euclidean, int, float> build_vp_tree() const {
    static const std::vector<bool> empty_vector;
    return build_vp_tree(empty_vector);
  }
  
  knncolle::VpTree<knncolle::distances::Euclidean, int, float> build_vp_tree(const std::vector<bool>& ix) const;
  
#endif    
};

