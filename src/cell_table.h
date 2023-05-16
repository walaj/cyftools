#pragma once
#include <functional>

#include "csv.h"
#include "cell_column.h"
#include "polygon.h"
#include "cell_header.h"
#include "cell_processor.h"

// Define the function wrapper type
typedef std::function<bool(const std::string& in, std::string& out)> LineStreamerWrapper;

class CellTable {
  
public:
  
  CellTable() {}

  CellTable(CellRowFunc func);

  void BuildTable(const std::string& file);

  void StreamTable(CellProcessor& proc, const std::string& file);
  
  // add columns
  void AddColumn(const std::string& key, ColPtr column); 

  void AddGraphColumn(const Tag& tag,
		     GraphColPtr value);

  void AddFlagColumn(const Tag& tag,
		     FlagColPtr value,
		     bool overwrite);

  void AddMetaColumn(const std::string& key, ColPtr value);

  // set params
  void SetPrecision(size_t n);
  
  void SetVerbose() { m_verbose = true; }

  void SetThreads(size_t threads) { m_threads = threads; }

  void SetPrintHeader() { m_print_header = true; }

  void SetHeaderOnly() { m_header_only = true; }
  
  IntColPtr GetIDColumn() const;
  
  friend std::ostream& operator<<(std::ostream& os, const CellTable& table);
  
  bool ContainsColumn(const std::string& name) const;

  size_t CellCount() const;

  int RadialDensity(uint64_t inner, uint64_t outer, uint64_t logor, uint64_t logand,
		    const std::string& label);
  
  const CellHeader& GetHeader() const;
  
  // display
  void PlotASCII(int width, int height) const;
  
  // IO
  void PrintHeader() const;

  void PrintTable(bool header) const;

  void ConvertColumns();
  
  // numeric operations
  void Log10();

  void PrintPearson(bool csv, bool sort) const;

  // graph ops
  void KNN_marker(int num_neighbors);

  void KNN_spatial(int num_neighbors, int dist);  

  // filtering
  void Subsample(int n, int s);

  void Crop(float xlo, float xhi, float ylo, float yhi);

  void SubsetROI(const std::vector<Polygon> &polygons);

  void Remove(const std::string& token);
  
  void Cut(const std::set<std::string>& tokens);

  void select(uint64_t on, uint64_t off);

  std::unordered_map<std::string, std::pair<float,float>> phenoread(const std::string& filename) const;
  
  void phenotype(const std::unordered_map<std::string, std::pair<float,float>>& thresh);
  
 private:
  
  unordered_map<string, ColPtr> m_table;
  
  string x;
  string y;

  CellHeader m_header;

  // params
  bool m_verbose = false;
  bool m_header_only = false;
  bool m_print_header = false;
  size_t m_threads = 1;
  
  // internal member functions
#ifndef DOXYGEN_SHOULD_SKIP_THIS
  
  bool read_csv_line__(io::LineReader& reader, CellRow& values) const;
  
  void header_read__(const std::string& header_line);
  
  CellRow add_row_to_table__(const CellRow& values);

  CellRow select_row_from_table__(const CellRow& values, uint64_t on,
					     uint64_t off);

  int read_one_line__(const std::string& line, CellRow& values) const;

  void check_header__() const;

  void print_correlation_matrix(const std::vector<std::pair<std::string, const ColPtr>>& data,
				const std::vector<std::vector<float>>& correlation_matrix, bool sort) const;

  
  void process_csv_file__(const std::string& file, const std::function<CellRow(const CellRow&)>& func);


  void column_to_row_major(std::vector<double>& data, int nobs, int ndim) const;
  
#endif    
};

struct RadialSelector {

  explicit RadialSelector(const std::string& line) {

    std::vector<std::string> tokens = split(line, ',');
    
    if (tokens.size() != 5) {
      throw std::runtime_error("There must be exactly 5 tokens.");
    }

    int_data.resize(4);
    for (int i = 0; i < 4; i++) {
      try {
	int_data[i] = std::stoi(tokens[i]);
      } catch (const std::invalid_argument &e) {
	throw std::runtime_error("The first 4 tokens must be integers.");
      }
    }

    label = tokens[4];
    
  }

  // split tokens
  std::vector<std::string> split(const std::string &input, char delimiter) {
    std::vector<std::string> tokens;
    std::istringstream stream(input);
    std::string token;
    
    while (std::getline(stream, token, delimiter)) {
      tokens.push_back(token);
    }
    
    return tokens;
  }

  friend std::ostream& operator<<(std::ostream& os, const RadialSelector& rs) {
    for (size_t i = 0; i < rs.int_data.size(); ++i) {
      os << rs.int_data[i];
      if (i < rs.int_data.size() - 1) {
	os << ",";
      }
    }
    os << "," << rs.label;
    return os;
  }
  
  std::vector<uint64_t> int_data;
  std::string label;
  
};
