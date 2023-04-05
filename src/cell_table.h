#ifndef CELL_TABLE_H
#define CELL_TABLE_H

#include "json_reader.h"
#include "cell_column.h"
#include "polygon.h"
#include "cell_header.h"

class CellTable {
  
public:
  
  CellTable() {}

  CellTable(int test);
  
  CellTable(const char* file, bool verbose, bool header_only);

  CellTable(const char* file, CellRowFunc func, bool verbose, bool header_only);
  
  void AddColumn(const std::string& key, std::shared_ptr<Column> column); 

  void AddGraphColumn(const Tag& tag,
		      const std::shared_ptr<StringColumn> value);
  
  void SetPrecision(size_t n);
  
  //CellTable(const char* file, const char* markers_file, bool verbose);
  
  friend std::ostream& operator<<(std::ostream& os, const CellTable& table);
  
  bool ContainsColumn(const std::string& name) const;

  size_t CellCount() const;

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

  // filtering
  void Subsample(int n, int s);

  void Crop(float xlo, float xhi, float ylo, float yhi);

  void SubsetROI(const std::vector<Polygon> &polygons);

  void Remove(const std::string& token);
  
  void Cut(const std::set<std::string>& tokens);


  
 private:
  
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

  int read_one_line__(const std::string& line, CellRow& values) const;

  void check_header__() const;

  //void print_correlation_matrix(const std::vector<std::pair<int, std::string>>& data, const std::vector<std::vector<double>>& correlation_matrix, bool sort) const;

  void print_correlation_matrix(const std::vector<std::pair<std::string, const std::shared_ptr<Column>>>& data,
				const std::vector<std::vector<float>>& correlation_matrix, bool sort) const;


  void process_csv_file__(const char* file,
			  const std::function<CellRow(const CellRow&)>& func,
			  bool verbose, bool header_only);
  
#endif    
};

#endif
