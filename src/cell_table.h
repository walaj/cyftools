#ifndef CELL_TABLE_H
#define CELL_TABLE_H

#include <functional>

#include "csv.h"
#include "cell_column.h"
#include "polygon.h"
#include "cell_header.h"

// Define the function wrapper type
typedef std::function<bool(const std::string& in, std::string& out)> LineStreamerWrapper;

class CellTable {
  
public:
  
  CellTable() {}

  CellTable(const char* file, bool verbose, bool header_only, bool convert);

  CellTable(const char* file, CellRowFunc func, bool verbose, bool header_only);

  CellTable(const char* file, bool verbose, bool header_only, uint64_t on, uint64_t off);

  // Stream select
  void StreamSelect(uint64_t on, uint64_t off);

  bool StreamSelect(const std::string& line_in, std::string& line_out,
		    uint64_t on, uint64_t off);
  
  void Streamer(bool print_header, const LineStreamerWrapper& func);  
  
  void AddColumn(const std::string& key, std::shared_ptr<Column> column); 

  void AddGraphColumn(const Tag& tag,
		      const std::shared_ptr<GraphColumn> value);

  void AddFlagColumn(const Tag& tag,
		     const std::shared_ptr<FlagColumn> value,
		     bool overwrite);
  
  void SetPrecision(size_t n);
  
  void SetVerbose() { m_verbose = true; }

  std::shared_ptr<NumericColumn<uint64_t>> GetIDColumn() const;
  
  friend std::ostream& operator<<(std::ostream& os, const CellTable& table);
  
  bool ContainsColumn(const std::string& name) const;

  size_t CellCount() const;

  int RadialDensity(uint64_t inner, uint64_t outer, uint64_t on, uint64_t off,
		    const std::string& label, bool verbose);
  
  const CellHeader& GetHeader() const;
  
  // insertion
  void AddMetaColumn(const std::string& key, const std::shared_ptr<Column> value);
  
  // display
  void PlotASCII(int width, int height) const;
  
  // IO
  void PrintHeader() const;

  void PrintTable(bool header) const;

  // numeric operations
  void Log10();

  void PrintPearson(bool csv, bool sort) const;

  // graph ops
  void KNN_marker(int num_neighbors, bool verbose, int threads);

  void KNN_spatial(int num_neighbors, int dist, bool verbose, int threads);  

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

  bool m_verbose = false;
  
  unordered_map<string, shared_ptr<Column>> m_table;
  
  string x;
  string y;

  CellHeader m_header;
  
  // internal member functions
#ifndef DOXYGEN_SHOULD_SKIP_THIS
  //void read_markers_json__(const char* markers_file);
  
  void verbose_line_read__(int count) const;
  
  bool read_csv_line__(io::LineReader& reader, CellRow& values) const;
  
  void header_read__(const std::string& header_line);
  
  CellRow add_row_to_table__(const CellRow& values);

  //CellRow add_row_to_table_by_index__(const CellRow& values, size_t ind, size_t dim);

  CellRow select_row_from_table__(const CellRow& values, uint64_t on,
					     uint64_t off);

  int read_one_line__(const std::string& line, CellRow& values) const;

  void check_header__() const;

  //void print_correlation_matrix(const std::vector<std::pair<int, std::string>>& data, const std::vector<std::vector<double>>& correlation_matrix, bool sort) const;

  void print_correlation_matrix(const std::vector<std::pair<std::string, const std::shared_ptr<Column>>>& data,
				const std::vector<std::vector<float>>& correlation_matrix, bool sort) const;


  void process_csv_file__(const char* file,
			  const std::function<CellRow(const CellRow&)>& func,
			  bool verbose, bool header_only);

  void convert_columns__();

  void column_to_row_major(std::vector<double>& data, int nobs, int ndim) const;
  
#endif    
};

#endif
