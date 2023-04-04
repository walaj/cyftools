#ifndef CELL_TABLE_H
#define CELL_TABLE_H

#include "json_reader.h"
#include "cell_column.h"
#include "polygon.h"
#include "cell_header.h"

class CellTable {
  
public:
  
  CellTable() {}

  CellTable(const char* file, bool verbose);
  
  void AddColumn(const std::string& key, std::shared_ptr<Column> column); 

  void SetPrecision(size_t n);
  
  CellTable(const char* file, const char* markers_file, bool verbose);
  
  friend std::ostream& operator<<(std::ostream& os, const CellTable& table);
  
  bool ContainsColumn(const std::string& name) const;

  size_t CellCount() const;

  // insertion
  void AddMetaColumn(const std::string& key, const std::shared_ptr<Column> value);
  
  // display
  void PlotASCII(int width, int height) const;
  
  // IO
  void PrintHeader() const;

  void PrintTable() const;

  // numeric operations
  void Log10();

  void PrintPearson(bool csv) const;
  
  // filtering
  void Subsample(int n, int s);

  void Crop(float xlo, float xhi, float ylo, float yhi);

  void SubsetROI(const std::vector<Polygon> &polygons);
  
 private:
  
  unordered_map<string, shared_ptr<Column>> m_table;
  
  vector<string> col_order;
  
  set<string> markers;
  set<string> meta;
  
  string x;
  string y;

  CellHeader m_header;
  
  // internal member functions
#ifndef DOXYGEN_SHOULD_SKIP_THIS
  void read_markers_json__(const char* markers_file);
  
  void verbose_line_read__(int count) const;
  
  bool read_csv_line__(io::LineReader& reader, CellRow& values) const;
  
  void header_read__(const std::string& header_line);
  
  void add_row_to_table__(const CellRow& values);

  void read_one_line__(const std::string& line, CellRow& values) const;

#endif    
};

#endif
