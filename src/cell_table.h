#pragma once

#include <set>

#include "cell_column.h"
#include "cell_header.h"
#include "cell_processor.h"
#include "cysift.h"

#include <cereal/types/vector.hpp>
#include <cereal/archives/portable_binary.hpp>
#include <cereal/archives/json.hpp>

// forward declarations
class Polygon;
class LDAModel;
class KDTree;
class TiffWriter;


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
  // Plotting
  //////
  int PlotPNG(const std::string& file) const;
  
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
		      std::vector<std::string> label, std::vector<int> normalize);
  
  void TumorCall(int num_neighbors, float frac,
		 cy_uint orflag, cy_uint andflag, cy_uint dist);
  

  //////
  // Numeric ops
  //////
  void PrintPearson(bool csv, bool sort) const;

  //////
  // Image ops
  //////
#ifdef HAVE_TIFFLIB
  void Convolve(TiffWriter* otoif, int boxwidth, float microns_per_pixel);
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

  void LDA_score_cells(const std::string& pdffile,
		       int topic_highlight,
		       float cont_cutoff); 
  
 private:

  unordered_map<string, ColPtr> m_table;

  // the cysift command
  std::string m_cmd;

  // IO
  std::unique_ptr<std::ofstream> m_os;
  std::unique_ptr<cereal::PortableBinaryOutputArchive> m_archive;

  CellHeader m_header;

  LDAModel* m_ldamodel = nullptr;

#ifdef HAVE_KDTREE
  KDTree* m_kdtree = nullptr;
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

