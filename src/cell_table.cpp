#include "cell_table.h"
#include "cell_graph.h"

#include <random>

const CellHeader& CellTable::GetHeader() const {
  return m_header;
}

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
    if (m_header.hasMarker(k.first)) {
      //    if (markers.count(k.first)) {
      if (k.second->GetType() == ColumnType::INT) {
	k.second = k.second->CopyToFloat();
      }
      k.second->Log10();
    }
  }
}

CellTable::CellTable(const char* file, CellRowFunc func, bool verbose, bool header_only) {
  process_csv_file__(file, func, verbose, header_only);
}

CellTable::CellTable(const char* file, bool verbose, bool header_only) {
  process_csv_file__(file, std::bind(&CellTable::add_row_to_table__, this, std::placeholders::_1), verbose, header_only);
}

void CellTable::verbose_line_read__(int count) const {
  
  // verbose output
  if (count % 100000 == 0) {
    std::cerr << "...reading line " << AddCommas(count) << std::endl;
  }
  
  return;
}

int CellTable::read_one_line__(const std::string& line,
				CellRow& values) const {
  
  char *cstr = new char[line.length() + 1];
  std::strcpy(cstr, line.c_str());
  char *value_str = strtok(cstr, ",");

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

  delete[] cstr;

  return n;
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

std::ostream& operator<<(std::ostream& os, const CellTable& table) {
  for (const auto& key : table.m_header.GetColOrder()) {
    auto col_ptr = table.m_table.at(key);
    std::string ctype;
    if (table.m_header.hasMarker(key))
      ctype = "Marker";
    else if (table.m_header.hasMeta(key))
      ctype = "Meta";
    else if (key == table.x || key == table.y)
      ctype = "Dim";
    else if (key == table.m_header.GetID())
      ctype = "ID";
    else
      ctype = "\nUNKNOWN COLUMN TYPE";
    os << key << " -- " << ctype << " -- " << col_ptr->toString() << std::endl;
  }
  return os;
}

CellRow CellTable::add_row_to_table__(const CellRow& values) {

  int col_count = m_header.ColumnCount();
  
  // make sure that the row is expected length
  // if m_header.col_order.size() is zero, you may not have read the header first
  if (values.size() != col_count) {
    throw std::runtime_error("Row size and expected size (from m_header.ColumnCount) must be the same.");
  }

  for (size_t i = 0; i < col_count; i++) {
    const std::string& col_name = m_header.GetColumnTag(i).GetName(); 
    const std::variant<int, float, std::string>& value = values.at(i);
	
    if (ContainsColumn(col_name)) {

      // NEW WAY 
      auto col_ptr = m_table[col_name];
      ColumnType col_type = col_ptr->GetType();

      switch (col_type) {
      case ColumnType::INT:
	if (std::holds_alternative<int>(value)) {
	  static_cast<NumericColumn<int>*>(col_ptr.get())->PushElem(std::get<int>(value));
	}
	break;
      case ColumnType::FLOAT:
	if (std::holds_alternative<float>(value)) {
	  static_cast<NumericColumn<float>*>(col_ptr.get())->PushElem(std::get<float>(value));
	} else if (std::holds_alternative<int>(value)) {
	  static_cast<NumericColumn<float>*>(col_ptr.get())->PushElem(static_cast<float>(std::get<int>(value)));
	}
	/*	if (std::holds_alternative<float>(value) || std::holds_alternative<int>(value)) {
		static_cast<NumericColumn<float>*>(col_ptr.get())->PushElem(std::get<float>(value));
	  }*/
	break;
      case ColumnType::STRING:
	if (std::holds_alternative<std::string>(value)) {
	  static_cast<StringColumn*>(col_ptr.get())->PushElem(std::get<std::string>(value));
	}
	break;
      default:
	// handle unknown column type
	std::cerr << "Error: column " << col_name << " has unknown type" << std::endl;
      }

    } else {
      // create new column with appropriate type
      if (std::holds_alternative<int>(value)) {
	//m_table[col_name] = std::make_shared<NumericColumn<int>>(std::get<int>(value));
	std::shared_ptr<NumericColumn<int>> new_ptr = std::make_shared<NumericColumn<int>>();
	new_ptr->PushElem(std::get<int>(value));
	m_table[col_name] = new_ptr; 
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

  return CellRow{}; // Return an empty CellRow
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
      std::cerr << "Warning: Column sizes do not match. Column: " <<
	c.first << " prev size " << prev_size.value() << " current_size " << current_size << std::endl;
    }
    
    prev_size = current_size;
    n = std::max(n, current_size);
  }
  
  return n;
}


void CellTable::Cut(const std::set<std::string>& tokens) {
  
  m_header.Cut(tokens);

  // cut the table itself
  std::unordered_map<std::string, std::shared_ptr<Column>> new_m_table;
  for (const auto& token : tokens) {
    auto it = m_table.find(token);
    if (it != m_table.end()) {
      new_m_table.emplace(token, it->second);
    } else {
      std::cerr << "Warning: Token '" << token << "' not found in m_table" << std::endl;
    }
  }
  
  m_table.swap(new_m_table);

}


void CellTable::PrintTable(bool header) const {

  if (header)
    m_header.Print();
  
  // Create a lookup table with pointers to the Column values
  std::vector<const std::shared_ptr<Column>> lookup_table;
  for (const auto& c : m_header.tags) {
    auto it = m_table.find(c.GetName());
    if (it != m_table.end()) {
      lookup_table.push_back(it->second);
    }
  }

  // Write the data (values)
  size_t numRows = CellCount();
  
  for (size_t row = 0; row < numRows; ++row) {
    size_t count = 1;
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

void CellTable::PrintPearson(bool csv, bool sort) const {

  // collect the marker data in one structure
  std::vector<std::pair<std::string, const std::shared_ptr<Column>>> data;
  for (const auto &t : m_table) {
    if (m_header.hasMarker(t.first)) {
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
    print_correlation_matrix(data, correlation_matrix, sort);
    /*

    
    for (int i = 0; i < n; i++) {
      for (int j = i + 1; j < n; j++) {
	std::cout << data[i].first << "," << data[j].first << 
	  "," << correlation_matrix[i][j] << std::endl;
      }
      }*/
    return;
  }
  
  // find the best spacing
  size_t spacing = 8;
  for (const auto& c : m_header.markers_)
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


void CellTable::PlotASCII(int width, int height) const { 

  if (width <= 0 || height <= 0) {
    std::cerr << "Warning: No plot generated, height and/or width <= 0" << std::endl;
    return;
  }

  // find the max size of the original cell table
  float x_min = m_table.at(x)->Min();
  float x_max = m_table.at(x)->Max();
  float y_min = m_table.at(y)->Min();
  float y_max = m_table.at(y)->Max();
  
  // scale it to fit on the plot
  size_t nc = CellCount();
  std::vector<std::pair<int, int>> scaled_coords;
  for (size_t i = 0; i < nc; i++) {
    float x_ = m_table.at(x)->GetNumericElem(i);
    float y_ = m_table.at(y)->GetNumericElem(i);    
    
    int x = static_cast<int>((x_ - x_min) / (x_max - x_min) * (width  - 2)) + 1;
    int y = static_cast<int>((y_ - y_min) / (y_max - y_min) * (height - 2)) + 1;
    scaled_coords.push_back({x, y});  //non-rotated    
					//scaled_coords.push_back({y, x}); // rotated
  }

  // make the plot pixel grid
  std::vector<std::string> grid(height, std::string(width, ' ')); // non-rotated
  //std::vector<std::string> grid(width, std::string(height, ' '));  // rotated

  // Draw the border
  for (int i = 0; i < width; i++) {
    grid[0][i] = '-';
    grid[height - 1][i] = '-';
  }
  for (int i = 1; i < height - 1; i++) {
    grid[i][0] = '|';
    grid[i][width - 1] = '|';
  }

  grid[0][0] = '+';
  grid[0][width - 1] = '+';
  grid[height - 1][0] = '+';
  grid[height - 1][width - 1] = '+';

  // Plot the cells
  for (const auto &coord : scaled_coords) {
    grid[coord.second][coord.first] = 'o';
  }
  
  // Print the grid
  for (const auto &row : grid) {
    std::cout << row << std::endl;
  }
  
}

void CellTable::AddGraphColumn(const Tag& tag,
			       const std::shared_ptr<StringColumn> value) {

  if (value->size() != CellCount()) {
    throw std::runtime_error("Adding graph column of incorrect size");
  }
  
  // check if already exists in the table
  if (m_table.find(tag.GetName()) != m_table.end()) {
    throw std::runtime_error("Adding graph column already in table");
  }
  
  // insert as a meta key
  m_header.addTag(tag); 

  // add to the table
  m_table[tag.GetName()] = value;
  
}
			       

void CellTable::AddMetaColumn(const std::string& key, const std::shared_ptr<Column> value) {

  // warn if creating a ragged table, but proceed anyway
  if (value->size() != CellCount()) {
    throw std::runtime_error("Adding meta column of incorrect size");
  }
  
  // check if already exists in the table
  if (m_table.find(key) != m_table.end()) {
    throw std::runtime_error("Adding meta column already in table");
  }

  // insert as a meta key
  m_header.addTag(Tag("CA",key));

  // add the table
  m_table[key] = value;
  
  return;
  
}

void CellTable::SubsetROI(const std::vector<Polygon> &polygons) {

  size_t nc = CellCount();
  
  // Initialize the "sparoi" column
  std::shared_ptr<NumericColumn<int>> new_data = std::make_shared<NumericColumn<int>>();
  new_data->reserve(nc);
  
  //Column new_data = std::vector<int>(nc, 0);
  //AddMetaColumn("sparoi", new_data); 
  
  // Loop the table and check if the cell is in the ROI
  for (size_t i = 0; i < nc; i++) {
    float x_ = m_table.at(x)->GetNumericElem(i);
    float y_ = m_table.at(y)->GetNumericElem(i);    

    // Loop through all polygons and check if the point is inside any of them
    for (const auto &polygon : polygons) {
      // if point is in this polygon, add the polygon id number to the roi
      if (polygon.PointIn(x_,y_)) {
	std::cerr << " ADDING POINT " << x_ << "," << y_ << std::endl;
	new_data->SetNumericElem(polygon.Id, i);
	
	//std::visit([i, id = polygon.Id](auto& vec) { vec[i] = id; }, table["sparoi"]);
	
        // uncomment below if want to prevent over-writing existing
        // break;
      }
    }
  }
}

void CellTable::check_header__() const {

  assert(m_header.validate());
  
}

void CellTable::KNN_marker(int num_neighbors, bool verbose, int threads) {

  int ndim = m_header.markers_.size();
  int nobs = CellCount();

  if (verbose)
    std::cerr << "...finding K nearest-neighbors on " << AddCommas(nobs) << " cells" << std::endl;  

  // column major the marker data
  std::vector<float> concatenated_data;
  
  for (const auto& key : m_header.markers_) {
    auto it = m_table.find(key);
    if (it != m_table.end()) {
      auto column_ptr = it->second;
      
      // Check if the column is a NumericColumn
      if (auto numeric_column_ptr = std::dynamic_pointer_cast<NumericColumn<float>>(column_ptr)) {
	const auto& data = numeric_column_ptr->getData();
	concatenated_data.insert(concatenated_data.end(), data.begin(), data.end());
      }
    }
  }

  if (verbose)
    std::cerr << "...setting up KNN graph on " << AddCommas(concatenated_data.size()) << " points" << std::endl;
  knncolle::AnnoyEuclidean<int, float> searcher(ndim, nobs, concatenated_data.data());

  umappp::NeighborList<float> output(nobs);

  if (verbose)
    std::cerr << "...building KNN graph" << std::endl;

  auto cell_graph = std::make_shared<CellGraph>();
  
  #pragma omp parallel for num_threads(threads)
  for (size_t i = 0; i < nobs; ++i) {
    if (i % 5000 == 0 && verbose)
      std::cerr << "...working on cell " <<
	AddCommas(i) << " with thread " <<
	omp_get_thread_num() << " K " << num_neighbors << std::endl;
    cell_graph->AddNode(CellNode(searcher.find_nearest_neighbors(i, num_neighbors)));
  }

  //
  Tag gtag("GA","knn");
  gtag.addValue("NN", std::to_string(num_neighbors));
  AddGraphColumn(gtag, cell_graph->toStringColumn());

}

void CellTable::print_correlation_matrix(const std::vector<std::pair<std::string, const std::shared_ptr<Column>>>& data,
					 //std::vector<std::pair<int, std::string>>& data,
					 const std::vector<std::vector<float>>& correlation_matrix, bool sort) const {
    int n = data.size();
    std::vector<std::tuple<int, int, double>> sorted_correlations;

    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            sorted_correlations.push_back(std::make_tuple(i, j, correlation_matrix[i][j]));
        }
    }

    if (sort) {
        std::sort(sorted_correlations.begin(), sorted_correlations.end(),
            [](const std::tuple<int, int, double>& a, const std::tuple<int, int, double>& b) {
                return std::get<2>(a) < std::get<2>(b);
            });
    }

    for (const auto& item : sorted_correlations) {
        int i = std::get<0>(item);
        int j = std::get<1>(item);
        double corr = std::get<2>(item);

        std::cout << data[i].first << "," << data[j].first << "," << corr << std::endl;
    }
}

void CellTable::process_csv_file__(const char* file,
				   const std::function<CellRow(const CellRow&)>& func,
				   bool verbose, bool header_only) {

  // make the csv reader
  io::LineReader reader(file, stdin);
  
  // for verbose
  size_t count = 0;
  
  // Read and process the lines
  std::string line;
  bool header_read = false;
  char* next_line_ptr;
  while ((next_line_ptr = reader.next_line()) != nullptr) {
    line = std::string(next_line_ptr);
    if (line.size() == 0) {
      break;
    }
    // If the line starts with '@', parse it as a Tag and add it to m_header
    if (line[0] == '@') {
      Tag tag(line);
      m_header.addTag(tag);
    } else {

      if (!header_read) {
	x = m_header.GetX();
	y = m_header.GetY();
	header_read = true;
	
	check_header__();
      }

      if (header_only) {
	break;
      }

      // container to store just one line
      int col_count = m_header.ColumnCount();
      CellRow values(col_count); 

      // read one line
      int num_elems = read_one_line__(line, values);
      if (num_elems != col_count) {
	throw std::runtime_error("Error on line: " + line + "\n\tread " +
				 std::to_string(num_elems) + " expected " +
				 std::to_string(col_count));
      }
      
      func(values);
      
      // will throw an error if detects type mismatch
      //add_row_to_table__(values);
      
      // verbose output
      if (verbose)
	verbose_line_read__(count);
      count++;
    }
  }
}
