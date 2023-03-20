#include "cell_table.h"

#include <random>

void CellTable::AddColumn(const std::string& key, std::shared_ptr<Column> column) {
  m_table[key] = column;
}

void CellTable::SetPrecision(size_t n) {
  for (auto& k : m_table) {
    k.second->SetPrecision(n);
  }
}

void CellTable::Log10() {
  for (auto& k : m_table) {
    if (markers.count(k.first)) {
      if (k.second->GetType() == ColumnType::INT) {
	k.second = k.second->CopyToFloat();
      }
      k.second->Log10();
    }
  }
}

CellTable::CellTable(const char* file, const char* markers_file, bool verbose) {

  // read markers first
  read_markers_json__(markers_file);

  // make the csv reader
  io::LineReader reader(file, stdin);
      
  // reader the header
  header_read__(reader.next_line());
      
  // container to store just one line
  CellRow values(col_order.size()); // m_table.size());
      
  // for verbose
  size_t count = 0;
      
  // read the data
  while (read_csv_line__(reader, values)) {

    // will throw an error if detects type mismatch
    add_row_to_table__(values);

    // verbose output
    if (verbose) 
      verbose_line_read__(count);
    count++;
	
  }
      
}


void CellTable::read_markers_json__(const char* markers_file) {
      
  // read the markers first
  std::ifstream       mfile(markers_file);
  std::string         line;
      
  // read in the markers
  std::string mm = std::string(markers_file);
  JsonReader json_reader(mm);
  json_reader.ReadData();
      
  x = json_reader.GetX();
  y = json_reader.GetY();
  assert(!x.empty());
  assert(!y.empty());  
      
  markers = json_reader.GetMarkers();
  meta = json_reader.GetMetaCells();  
      
}

void CellTable::verbose_line_read__(int count) const {
  
  // verbose output
  if (count % 100000 == 0) {
    std::cerr << "...reading line " << AddCommas(count) << std::endl;
  }
  
  return;
}

bool CellTable::read_csv_line__(io::LineReader& reader,
				CellRow& values) const {
      
  // Read the data lines
  while (char* line = reader.next_line()) {
    char *value_str = strtok(line, ",");
	
    size_t n = 0;
    while (value_str) {
	  
      std::string val(value_str);
      std::istringstream iss(val);
      float float_val;
      int int_val;
      char decimal_point;
      char extra_char;

      if (iss >> int_val && !(iss >> decimal_point)) {

	// Integer
	values[n] = int_val;
      } else {

	iss.clear();
	iss.str(val);
	if (iss >> float_val && !(iss >> extra_char)) {
	  // Float
	  values[n] = float_val;
	} else {
	  // String
	  values[n] = val;
	}
      }
	  
      value_str = strtok(nullptr, ",");
      n++;
    }
	
    return true;
  }
      
  return false;
}

bool CellTable::ContainsColumn(const std::string& name) const {
  return m_table.count(name) > 0;
}

void CellTable::header_read__(const std::string& header_line) {
      
  std::istringstream header_stream(header_line);
  std::string header_item;
      
  // store the markers, meta etc in an ordered way
  // this is needed to link the vec's to the right
  // table element, since unordered_map doesn't store
  // by order of insertion
  while (std::getline(header_stream, header_item, ',')) {
	
    // remove any special characters 
    header_item.erase(std::remove(header_item.begin(), header_item.end(), '\r'), header_item.end());
	
    // store the original order
    col_order.push_back(header_item);
  }
}


std::ostream& operator<<(std::ostream& os, const CellTable& table) {
  for (const auto& key : table.col_order) {
    auto col_ptr = table.m_table.at(key);
    os << key << " -- " << col_ptr->toString() << std::endl;
  }
  return os;
}


void CellTable::add_row_to_table__(const CellRow& values) {
      
  // make sure that the row is expected length
  // if col_order.size() is zero, you may not have read the header first
  if (values.size() != col_order.size()) {
    throw std::runtime_error("Row size and expected size (from col_order) must be the same.");
  }
      
  for (size_t i = 0; i < col_order.size(); i++) {
    const std::string& col_name = col_order[i];
    const std::variant<int, float, std::string>& value = values.at(i);
	
    if (ContainsColumn(col_name)) {
	  
      // NEW WAY 
      auto col_ptr = m_table[col_name];
      ColumnType col_type = col_ptr->GetType();
      switch (col_type) {
      case ColumnType::INT:
	if (std::holds_alternative<int>(value)) {
	  static_cast<NumericColumn<int>*>(col_ptr.get())->AddElem(std::get<int>(value));
	}
	break;
      case ColumnType::FLOAT:
	if (std::holds_alternative<float>(value)) {
	  static_cast<NumericColumn<float>*>(col_ptr.get())->AddElem(std::get<float>(value));
	}
	break;
      case ColumnType::STRING:
	if (std::holds_alternative<std::string>(value)) {
	  static_cast<StringColumn*>(col_ptr.get())->AddElem(std::get<std::string>(value));
	}
	break;
      default:
	// handle unknown column type
	std::cerr << "Error: column " << col_name << " has unknown type" << std::endl;
      }

    } else {
      // create new column with appropriate type
      if (std::holds_alternative<int>(value)) {
	m_table[col_name] = std::make_shared<NumericColumn<int>>(std::get<int>(value));
      } else if (std::holds_alternative<float>(value)) {
	m_table[col_name] = std::make_shared<NumericColumn<float>>(std::get<float>(value));
      } else if (std::holds_alternative<std::string>(value)) {
	m_table[col_name] = std::make_shared<StringColumn>(std::get<std::string>(value));
      } else {
	// handle unknown column type
	std::cerr << "Error: column " << col_name << " has unknown type" << std::endl;
      }
	  
    }
  }
      
}

void CellTable::Subsample(int n, int s) {
  
  // Create a random number generator with the provided seed
  std::default_random_engine generator(s);
  
  // Find the number of rows in the table
  if (m_table.empty()) {
    return;
  }
  
  size_t num_rows = CellCount();
  
  // Make sure n is within the range [0, num_rows]
  n = std::min(n, static_cast<int>(num_rows));
  
  // Generate a vector of row indices
  std::vector<size_t> row_indices(num_rows);
  for (size_t i = 0; i < num_rows; ++i) {
    row_indices[i] = i;
  }
  
  // Shuffle the row indices
  std::shuffle(row_indices.begin(), row_indices.end(), generator);
  
  // Select the first n indices after shuffling
  row_indices.resize(n);
  
  // sort the indicies so order is not scrambled
  std::sort(row_indices.begin(), row_indices.end());
  
  // Subsample the table using the selected row indices
  for (auto &kv : m_table) {
    kv.second->SubsetColumn(row_indices);
    // new_column[i] = kv.second[row_indices[i]];
      
      // std::visit(
      // 		 [&row_indices, n](auto &vec) {
      // 		   using ValueType = typename std::decay_t<decltype(vec)>::value_type;
      // 		   std::vector<ValueType> new_column(n);
      // 		   for (size_t i = 0; i < n; ++i) {
      // 		     new_column[i] = vec[row_indices[i]];
      // 		   }
      // 		   vec = std::move(new_column);
      // 		 },
      // 		 kv.second);
    }

    return;
}

size_t CellTable::CellCount() const {

  std::optional<size_t> prev_size;
  size_t n = 0;
  
  // loop the table and find the maximum length
  for (const auto &c : m_table) {
    size_t current_size = c.second->size();
    
    if (prev_size.has_value() && current_size != prev_size.value()) {
      std::cerr << "Warning: Column sizes do not match. Column: " << c.first << std::endl;
    }
    
    prev_size = current_size;
    n = std::max(n, current_size);
  }
  
  return n;
}

void CellTable::PrintHeader() const {

  // Write the header (keys)
  size_t count = 0;
  for (auto& c: col_order) {
    count++;
    auto it = m_table.find(c);
    
    if (it == m_table.end()) {
      throw std::runtime_error("Trying to print header name without corresponding table data");
      continue;
    }
    
    std::cout << it->first;
    //if (std::next(it) != table.end()) {
    if (count != col_order.size())
      std::cout << ",";
  }

  return;
  
}

void CellTable::PrintTable() const {
  
  PrintHeader();
  std::cout << std::endl;
  
  // Create a lookup table with pointers to the Column values
  std::vector<const std::shared_ptr<Column>> lookup_table;
  for (const auto& c : col_order) {
    auto it = m_table.find(c);
    if (it != m_table.end()) {
      lookup_table.push_back(it->second);
    }
  }

  // Write the data (values)
  size_t numRows = CellCount();

  size_t count = 0;
  for (size_t row = 0; row < numRows; ++row) {
    
    for (const auto& c : lookup_table) {
      c->PrintElem(row);
      if (count < lookup_table.size())
	std::cout << ",";
      count++;
    }
    std::cout << std::endl;
  }
}



void CellTable::Crop(float xlo, float xhi, float ylo, float yhi) {
  
  // Find the x and y columns in the table
  auto x_it = m_table.find(x);
  auto y_it = m_table.find(y);

  // total number of rows of table
  size_t nc = CellCount();

  // error checking
  if (x_it == m_table.end() || y_it == m_table.end()) 
    throw std::runtime_error("x or y column not found in the table");
  
  std::vector<size_t> valid_indices;
  for (size_t i = 0; i < nc; ++i) {
    
    float x_value = x_it->second->GetNumericElem(i);
    float y_value = y_it->second->GetNumericElem(i);    
    
    if (x_value >= xlo && x_value <= xhi &&
	y_value >= ylo && y_value <= yhi) {
      valid_indices.push_back(i);
    }
  }
  
  // Update the table to include only the valid_indices
  for (auto &kv : m_table) {
    kv.second->SubsetColumn(valid_indices);
  }
}



void CellTable::PrintPearson(bool csv) const {

  // collect the marker data in one structure
  std::vector<std::pair<std::string, const std::shared_ptr<Column>>> data;  
  for (const auto &t : m_table) {
    if (markers.count(t.first)) {
      data.push_back({t.first, t.second});
    }
  }
  
  // get the correlation matrix
  size_t n = data.size();
  std::vector<std::vector<float>> correlation_matrix(n, std::vector<float>(n, 0));
  for (size_t i = 0; i < n; i++) {
    for (size_t j = i + 1; j < n; j++) {
      correlation_matrix[i][j] = data[i].second->Pearson(*data[j].second);
    }
  }
  
  // if csv, print to csv instead of triangle table
  if (csv) {
    for (int i = 0; i < n; i++) {
      for (int j = i + 1; j < n; j++) {
	std::cout << data[i].first << "," << data[j].first << 
	  "," << correlation_matrix[i][j] << std::endl;
      }
    }
    return;
  }
  
  // find the best spacing
  size_t spacing = 8;
  for (const auto& c : markers)
    if (c.length() > spacing)
      spacing = c.length() + 2;
  
  // Print x-axis markers
  std::cout << std::setw(spacing) << " ";
  for (size_t i = 0; i < n; i++) {
    std::cout << std::setw(spacing) << data[i].first;
  }
  std::cout << std::endl;

  // print the matrix
  for (size_t i = 0; i < n; i++) {
    
    // Print y-axis marker
    std::cout << std::setw(spacing) << data[i].first;
    
    for (size_t j = 0; j < n; j++) {
      if (j < i) {
	std::cout << std::setw(spacing) << std::fixed << std::setprecision(4) << correlation_matrix[j][i];	
      } else {
	std::cout << std::setw(spacing) << " ";
      }
    }
    std::cout << std::endl;
  }
  
}
