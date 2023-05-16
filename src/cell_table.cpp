#include <random>
#include <cstdlib>

#include "cell_utils.h"
#include "cell_table.h"
#include "cell_graph.h"

static size_t debugr = 0;

const CellHeader& CellTable::GetHeader() const {
  return m_header;
}

void CellTable::AddColumn(const std::string& key, ColPtr column) {
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
      if (k.second->GetType() == ColumnType::INT) {
	k.second = k.second->CopyToFloat();
      }
      k.second->Log10();
    }
  }
}

void CellTable::BuildTable(const std::string& file) {
  process_csv_file__(file, [this](const CellRow& values) {
    return this->add_row_to_table__(values);
  });
}

bool CellTable::ContainsColumn(const std::string& name) const {
  return m_table.count(name) > 0;
}

std::ostream& operator<<(std::ostream& os, const CellTable& table) {
  
  for (const auto& key : table.m_header.GetColTags()) {
    const std::string tagn = key.GetName();
    auto col_ptr = table.m_table.at(tagn);
    std::string ctype;
    if (table.m_header.hasMarker(tagn))
      ctype = "Marker";
    else if (table.m_header.hasMeta(tagn))
      ctype = "Meta";
    else if (table.m_header.hasFlag(tagn))
      ctype = "Flag";
    else if (table.m_header.hasGraph(tagn))
      ctype = "Graph";
    else if (tagn == table.x || tagn == table.y)
      ctype = "Dim";
    else if (tagn == table.m_header.GetID())
      ctype = "ID";
    else
      ctype = "\nUNKNOWN COLUMN TYPE";
    os << tagn << " -- " << ctype << " -- " << col_ptr->toString() << std::endl;
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
    const std::variant<uint64_t, float, string>& value = values.at(i);
    
    // column already added
    if (ContainsColumn(col_name)) {

      // NEW WAY 
      auto col_ptr = m_table[col_name];
      ColumnType col_type = col_ptr->GetType();

      switch (col_type) {
      case ColumnType::INT:
	if (std::holds_alternative<uint64_t>(value)) 
	  static_cast<IntCol*>(col_ptr.get())->PushElem(std::get<uint64_t>(value));
	break;

      case ColumnType::FLOAT:

	if (std::holds_alternative<float>(value)) {
	  static_cast<FloatCol*>(col_ptr.get())->PushElem(std::get<float>(value));
	} else if (std::holds_alternative<uint64_t>(value)) {
	  static_cast<FloatCol*>(col_ptr.get())->PushElem(static_cast<float>(std::get<uint64_t>(value)));
	}
	break;
	
      case ColumnType::STRING:
	if (std::holds_alternative<std::string>(value))
	  static_cast<StringColumn*>(col_ptr.get())->PushElem(std::get<std::string>(value));
	break;
      default:
	std::cerr << "Error: column " << col_name << " has unknown type" << std::endl;
      }

    } else {

      // lookup what type of tag this is
      const Tag& tag = m_header.GetColTags().at(i);
      
      // create new column with appropriate type
      if (std::holds_alternative<uint64_t>(value)) {
	//      if (std::holds_alternative<uint64_t>(value) && !tag.isMarkerTag()) {	
	IntColPtr new_ptr = std::make_shared<IntCol>();
	new_ptr->PushElem(std::get<uint64_t>(value));
	m_table[col_name] = new_ptr; 
	//      } else if (std::holds_alternative<float>(value) || tag.isMarkerTag()) {
      } else if (std::holds_alternative<float>(value)) {	

	// enforce that all marker columns are float
	//if (std::holds_alternative<uint64_t>(value))
	//  m_table[col_name] = std::make_shared<FloatCol>(static_cast<float>(std::get<uint64_t>(value)));
	//else
	m_table[col_name] = std::make_shared<FloatCol>(std::get<float>(value));
	
      } else if (std::holds_alternative<std::string>(value)) {
	m_table[col_name] = std::make_shared<StringColumn>(std::get<std::string>(value));
      } else {
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

void CellTable::PrintTable(bool header) const {

  if (m_verbose)
    std::cerr << "...outputting table" << std::endl;
  
  if (header)
    m_header.Print();
  
  // Create a lookup table with pointers to the Column values
  std::vector<ColPtr> lookup_table;
  for (const auto& c : m_header.GetColTags()) {
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
  std::vector<std::pair<std::string, const ColPtr>> data;
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
			       GraphColPtr value) {

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
			       

void CellTable::AddMetaColumn(const std::string& key,
			      ColPtr value) {

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
  IntColPtr new_data = std::make_shared<IntCol>();
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

void CellTable::KNN_spatial(int num_neighbors, int dist) {

  // number of cells
  int nobs = CellCount();

  if (m_verbose)
    std::cerr << "...finding K nearest-neighbors (spatial) on " << AddCommas(nobs) << " cells" << std::endl;  

  // column major the coordinate data
  std::vector<double> concatenated_data;

  int ndim = 0;
  if (m_table.find(m_header.x_) != m_table.end()) {
    const auto& data = std::dynamic_pointer_cast<FloatCol>(m_table.at(m_header.x_))->getData();
    if (m_verbose) 
      std::cerr << "...adding " << AddCommas(data.size()) << " points on x" << std::endl;
    concatenated_data.insert(concatenated_data.end(), data.begin(), data.end());
    ndim++;
  }
  if (m_table.find(m_header.y_) != m_table.end()) {
    const auto& data = std::dynamic_pointer_cast<FloatCol>(m_table.at(m_header.y_))->getData();
    if (m_verbose) 
      std::cerr << "...adding " << AddCommas(data.size()) << " points on y" << std::endl;
    concatenated_data.insert(concatenated_data.end(), data.begin(), data.end());
    ndim++;    
  }
  if (m_table.find(m_header.z_) != m_table.end()) {
    const auto& data = std::dynamic_pointer_cast<FloatCol>(m_table.at(m_header.z_))->getData();
    if (m_verbose) 
      std::cerr << "...adding " << AddCommas(data.size()) << " points on z" << std::endl;
    concatenated_data.insert(concatenated_data.end(), data.begin(), data.end());
    ndim++;    
  }

  // convert to row major?
  column_to_row_major(concatenated_data, nobs, ndim);

  // get the cell id colums
  auto id_ptr = GetIDColumn();
  
  if (m_verbose)
    std::cerr << "...setting up KNN graph (spatial) on " << AddCommas(nobs) << " points" << std::endl;

  umappp::NeighborList<float> output(nobs);

  if (m_verbose)
    std::cerr << "...building KNN (spatial) graph" << std::endl;

  // store the final graph
  auto graph = std::make_shared<GraphColumn>();
  graph->resize(nobs);

  // track number of cases where all N nearest neighbors are within the
  // distance cutoff, implying that there are likely additional cells within
  // the desired radius that are cutoff
  size_t lost_cell = 0;

  // initialize the tree. Can choose from different algorithms, per knncolle library
  knncolle::VpTreeEuclidean<int, double> searcher(ndim, nobs, concatenated_data.data());
  //knncolle::AnnoyEuclidean<int, double> searcher(ndim, nobs, concatenated_data.data());  
  
#pragma omp parallel for num_threads(m_threads)
  for (size_t i = 0; i < nobs; ++i) {
    if (i % 50000 == 0 && m_verbose)
      std::cerr << "...working on cell " <<
	AddCommas(i) << " with thread " <<
	omp_get_thread_num() << " K " << num_neighbors << " Dist: " << dist <<
	std::endl;
    
    Neighbors neigh = searcher.find_nearest_neighbors(i, num_neighbors);

    // set the node pointer to CellID rather than 0-based
    for (auto& cc : neigh) {
      cc.first = id_ptr->GetNumericElem(cc.first);
    }
    
    // remove less than distance
    if (dist > 0) {
      size_t osize = neigh.size();
      Neighbors neigh_trim;
      for (const auto& nnn : neigh) {
	if (nnn.second < dist)
	  neigh_trim.push_back(nnn);
      }
      
      if (osize == neigh_trim.size()) {
	      #pragma omp critical
	{
	lost_cell++;
	if (lost_cell % 500 == 0)
	  std::cerr << "osize " << osize << " Lost cell " << AddCommas(lost_cell) << " of " << AddCommas(nobs) << std::endl;
	}
      }

      #pragma omp critical
      {
	graph->SetValueAt(i, CellNode(neigh_trim));
      }
      
    // no distance limitation
    } else {
#pragma omp critical
      {
	graph->SetValueAt(i, CellNode(neigh));
      }
    }
  }// end for

  if (m_verbose)
    std::cerr << "...done with graph construction" << std::endl;

  // force it to print distances as integers to save space
  graph->SetIntegerize(true);
  
  //
  Tag gtag("GA","spat");
  gtag.addValue("NN", std::to_string(num_neighbors));
  AddGraphColumn(gtag, graph);
}

void CellTable::KNN_marker(int num_neighbors) {

  int ndim = m_header.markers_.size();
  int nobs = CellCount();

  if (m_verbose)
    std::cerr << "...finding K nearest-neighbors on " << AddCommas(nobs) << " cells" << std::endl;  

  // column major the marker data
  std::vector<double> concatenated_data;
  
  for (const auto& key : m_header.markers_) {
    auto it = m_table.find(key);
    if (it != m_table.end()) {
      auto column_ptr = it->second;
      
      // Check if the column is a NumericColumn
      if (auto numeric_column_ptr = std::dynamic_pointer_cast<FloatCol>(column_ptr)) {
	const auto& data = numeric_column_ptr->getData();
	concatenated_data.insert(concatenated_data.end(), data.begin(), data.end());
      }
    }
  }

  // convert to row major?
  column_to_row_major(concatenated_data, nobs, ndim);

  // get the cell id colums
  auto id_ptr = GetIDColumn();
  
  if (m_verbose)
    std::cerr << "...setting up KNN graph on " << AddCommas(concatenated_data.size()) << " points" << std::endl;

  umappp::NeighborList<double> output(nobs);

  if (m_verbose)
    std::cerr << "...building KNN (marker) graph" << std::endl;

  auto graph = std::make_shared<GraphColumn>();
  graph->resize(nobs);

  // initialize the tree. Can choose from different algorithms, per knncolle library
  knncolle::VpTreeEuclidean<int, double> searcher(ndim, nobs, concatenated_data.data());
  //knncolle::AnnoyEuclidean<int, double> searcher(ndim, nobs, concatenated_data.data());  
  
  #pragma omp parallel for num_threads(m_threads)
  for (size_t i = 0; i < nobs; ++i) {
    if (i % 50000 == 0 && m_verbose)
      std::cerr << "...working on cell " <<
	AddCommas(i) << " with thread " <<
	omp_get_thread_num() << " K " << num_neighbors << std::endl;

    Neighbors neigh = searcher.find_nearest_neighbors(i, num_neighbors);

    // set the node pointer to CellID rather than 0-based
    for (auto& cc : neigh) {
      cc.first = id_ptr->GetNumericElem(cc.first);
    }
    
#pragma omp critical
    {
      graph->SetValueAt(i, CellNode(neigh));
    }
  }

  //
  Tag gtag("GA","knn");
  gtag.addValue("NN", std::to_string(num_neighbors));
  AddGraphColumn(gtag, graph);

}

void CellTable::print_correlation_matrix(const std::vector<std::pair<std::string, const ColPtr>>& data,
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

void CellTable::process_csv_file__(const std::string& file, const std::function<CellRow(const CellRow&)>& func) {

  // csv reader
  std::unique_ptr<io::LineReader> reader;

  // Initialize the reader depending on the input string
  if (file == "-") {
    // Read from stdin
    reader = std::make_unique<io::LineReader>("", stdin);
  } else {
    // Read from file
    reader = std::make_unique<io::LineReader>(file);
  }

  // for m_verbose
  size_t count = 0;
  
  // Read and process the lines
  std::string line;
  bool header_read = false;
  char* next_line_ptr;
  while ((next_line_ptr = reader->next_line()) != nullptr) {
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

      // if header only, we're done
      if (m_header_only) {
	break;
      }
      
      // container to store just one line
      int col_count = m_header.ColumnCount();
      CellRow values(col_count); 

      // read one line
      int num_elems = read_one_line_to_cellrow(line, values, m_header);
      if (num_elems != col_count) {
	throw std::runtime_error("Error on line: " + line + "\n\tread " +
				 std::to_string(num_elems) + " expected " +
				 std::to_string(col_count));
      }

      func(values);
      
      // m_verbose output
      if (m_verbose && count % 500000 == 0) {
	std::cerr << "...read line " << AddCommas(count) << std::endl;
      }
      count++;
      
    }
  }
}

void CellTable::ConvertColumns() {

  // convert flag columns
  for (const auto& f : m_header.flag_) {

    auto it = m_table.find(f);
    
    // make sure we have data to convert
    if (it == m_table.end())
      throw std::runtime_error("CellTable convert columns error, expected flag column not found");

    ColPtr column = it->second;
    ColumnType ct = column->GetType();
    StringColPtr sc;
    IntColPtr nc;

    switch (ct) {
    case ColumnType::FLAG: {
      std::cerr << "Flag column already converted" << std::endl; break;
    } case ColumnType::STRING: {
	throw std::runtime_error("Flag column should be read as int, not string");
	break;
      } case ColumnType::INT: {
	  nc =  std::dynamic_pointer_cast<IntCol>(column);
	  if (!nc)
	    throw std::runtime_error("Column is not a NumericColumn as expected");
	  m_table[f] = std::make_shared<FlagColumn>(nc);
	  break;
	}  default:
      throw std::runtime_error("Trying to convert non-string column to flag");
    }
    
  }
  
  // convert graph columns
  if (m_verbose)
    std::cerr << "...converting graph columns" << std::endl;
  
  for (const auto& f : m_header.graph_) {
    
    auto it = m_table.find(f);
    
    // make sure we have data to convert
    if (it == m_table.end())
      throw std::runtime_error("CellTable convert columns error, expected graph column not found");

    ColPtr column = it->second;
    ColumnType ct = column->GetType();
    StringColPtr sc;
    
    switch (ct) {
    case ColumnType::GRAPH: {
      std::cerr << "Graph column already converted" << std::endl; break;
    } case ColumnType::STRING: {
      sc = std::dynamic_pointer_cast<StringColumn>(column);
      if (!sc)
	throw std::runtime_error("Column is not a StringColumn as expected");
      m_table[f] = std::make_shared<GraphColumn>(sc, m_threads);
      break;
    } default:
      throw std::runtime_error("Trying to convert non-string column to graph");
    }

  }
}

void CellTable::select(uint64_t on, uint64_t off) {

  FlagColPtr fc = std::dynamic_pointer_cast<FlagColumn>(m_table["cellflag"]);

  assert(fc->size());
  std::vector<size_t> s;

  for (size_t i = 0; i < fc->size(); i++) {

    if (fc->TestFlag(on, off, i))
      s.push_back(i);
  }

  for (auto& m : m_table) {
    m.second->SubsetColumn(s);
  }
  
  return;
}

void CellTable::phenotype(const std::unordered_map<std::string, std::pair<float,float>>& thresh) {

  FlagColPtr fc = std::make_shared<FlagColumn>();
  fc->resize(CellCount());
  
  for (const auto& b : thresh) {
    
    if (m_table.count(b.first) == 0) {
      continue;
    }
    
    if (m_table.at(b.first)->GetType() != ColumnType::FLOAT) {
      std::cerr << "Error: Expected marker from threshold table not a float marker: " << b.first << std::endl;
      continue;
    }

    FloatColPtr nc = std::dynamic_pointer_cast<FloatCol>(m_table.at(b.first));

    int index = m_header.whichMarkerColumn(b.first);
        
    for (size_t i = 0; i < nc->size(); i++) {
      if (nc->GetNumericElem(i) >= b.second.first &&
	  nc->GetNumericElem(i) <= b.second.second) {
	fc->SetFlagOn(index, i);
      }
    }
  }

  Tag ftag("FA","cellflag");
  AddFlagColumn(ftag, fc, true);
}

std::unordered_map<std::string, std::pair<float,float>> CellTable::phenoread(const std::string& filename) const {
  
  std::unordered_map<std::string, std::pair<float, float>> data;
  std::ifstream file(filename);
  
  if (!file.is_open()) {
    std::cerr << "Error opening file: " << filename << std::endl;
    return data;
  }
  
  std::string line;
  while (std::getline(file, line)) {
    std::istringstream lineStream(line);
        std::string key;
        float value1, value2;
	
        if (std::getline(lineStream, key, ',')) {
	  lineStream >> value1;
	  if (lineStream.get() == ',') {
                lineStream >> value2;
                data[key] = std::make_pair(value1, value2);
	  }
        }
  }
  
  file.close();
  return data;
  
}


void CellTable::AddFlagColumn(const Tag& tag,
			      FlagColPtr value,
			      bool overwrite) {

  if (value->size() != CellCount()) {
    throw std::runtime_error("Adding flag column of incorrect size");
  }

 
  
  // check if already exists in the table
  if (m_table.find(tag.GetName()) != m_table.end()) {
    if (!overwrite)
      throw std::runtime_error("Adding flag column already in table");
  }
  
  // insert as a meta key
  if (m_table.find(tag.GetName()) == m_table.end())
    m_header.addTag(tag); 

  // add to the table
  m_table[tag.GetName()] = value;
  
}

void CellTable::column_to_row_major(std::vector<double>& data, int nobs, int ndim) const {

  double* temp = new double[data.size()];
  
  for (int row = 0; row < nobs; ++row) {
    for (int col = 0; col < ndim; ++col) {
      int old_index = col * nobs + row;
            int new_index = row * ndim + col;
            temp[new_index] = data[old_index];
    }
  }
  
  for (size_t i = 0; i < data.size(); ++i) {
    data[i] = temp[i];
  }
  
  delete[] temp;
  
}

IntColPtr CellTable::GetIDColumn() const {

  for (const auto& t: m_header.GetColTags()) {
    if (t.isIDTag()) {
      return std::dynamic_pointer_cast<IntCol>(m_table.at(t.GetName()));
    }
  }

  assert(false);
  return IntColPtr();
}

void CellTable::StreamTable(CellProcessor& proc, const std::string& file) {

  // csv reader
  std::unique_ptr<io::LineReader> reader;

  // Initialize the reader depending on the input string
  if (file == "-") {
    // Read from stdin
    reader = std::make_unique<io::LineReader>("", stdin);
  } else {
    // Read from file
    reader = std::make_unique<io::LineReader>(file);
  }
  
  // Read and process the lines
  std::string line;
  bool header_read = false;
  size_t flag_index = -1;
  char* next_line_ptr;

  // for verbose
  size_t count = 0;
  
  while ((next_line_ptr = reader->next_line()) != nullptr) {

    // get the line. Skip if empty
    line = std::string(next_line_ptr);
    if (line.size() == 0) {
      break;
    }
    
    // If the line starts with '@', parse it as a Tag and add it to m_header
    if (line[0] == '@') {
      
      // make sure there are no header lines in the middle of file
      if (header_read)
	throw std::runtime_error("Misformed file: header lines should all be at top of file");
      
      // add the tag to the header
      Tag tag(line);
      m_header.addTag(tag);
      
    } else {
      
      // label that the header has already been read
      if (!header_read) {
	header_read = true;

	// process the header
	// if it gives non-zero output, it means stop there
	if (proc.ProcessHeader(m_header))
	  return;
      }

      // process lines
      proc.ProcessLine(line);

      // m_verbose output
      if (m_verbose && count % 50000 == 0) {
	std::cerr << "...read line " << AddCommas(count) << std::endl;
      }
      count++;

      
    }
  }
}

int CellTable::RadialDensity(uint64_t inner, uint64_t outer, uint64_t logor, uint64_t logand,
			       const std::string& label) {
  if (m_table.find("cellflag") == m_table.end() && (logor || logand)) {
    std::cerr << "Warning: need to phenotype cells first" << std::endl;
    return -1;
  }

  if (m_table.find("spat") == m_table.end()) {
    std::cerr << "Warning: need to make spatial knn graph first" << std::endl;
    return -1;
  }

  // check if already exists in the table
  if (m_table.find(label) != m_table.end()) {
    throw std::runtime_error("Adding density meta column " + label + " already in table");
  }

  // get the flag column
  shared_ptr<FlagColumn> fc = std::dynamic_pointer_cast<FlagColumn>(m_table["cellflag"]);

  // get the graph column
  shared_ptr<GraphColumn> gc = std::dynamic_pointer_cast<GraphColumn>(m_table["spat"]);

  // check sizes are good
  assert(gc->size() == fc->size());
  
  // store the densities
  //shared_ptr<FloatCol> dc = std::make_shared<FloatCol>();
  shared_ptr<FloatCol> dc = std::make_shared<FloatCol>();  
  dc->resize(gc->size());

  dc->SetPrecision(2);

  // get the cell id colums
  auto id_ptr = GetIDColumn();
  assert(id_ptr->size() == gc->size());

  // flip it so that the id points to the zero-based index
  std::unordered_map<int, size_t> inverse_lookup;
  for (size_t i = 0; i < id_ptr->size(); ++i) {
    int value = id_ptr->GetNumericElem(i);
    inverse_lookup[value] = i;
  }
  
  if (m_verbose)
    std::cerr << "...radial density: starting loop" << std::endl;
  
  // loop the cells
  //  #pragma omp parallel for num_threads(m_threads)
  for (size_t i = 0; i < gc->size(); i++) {
    
    float cell_count = 0;

    const CellNode& node = gc->GetNode(i);

    const Neighbors& neigh = node.get_neighbors();
    
    // loop the nodes connected to each cell
    for (const auto& n : neigh) {

      assert(inverse_lookup.count(n.first));
      int cellindex = inverse_lookup[n.first];

      //debug
      //assert(!fc->TestFlagAndOr(logor, logand, cellindex));
      
      // test if the connected cell meets the flag criteria
      // n.first is cell_id of connected cell to this cell
      if ((!logor && !logand) || fc->TestFlagAndOr(logor, logand, cellindex)) {

	// if it meets flag criteria, test if if it's in the radius bounds
	if (n.second >= inner && n.second <= outer) {
	  cell_count++;	  
	}
	
      }
    }

    // calculate the density
    float area = static_cast<float>(outer) * static_cast<float>(outer) * 3.1415926535;
    area = area - static_cast<float>(inner) * static_cast<float>(inner) * 3.1415926535;
    
    if (neigh.size()) {
      dc->SetNumericElem(cell_count * 1000000 / area, i); // density per 1000 square pixels
    } else {
      dc->SetNumericElem(0, i); 
    }

    // veborse
    if (m_verbose && i % 50000 == 0)
      std::cerr << "...radial density: looping on " << AddCommas(i) << " and found density " << dc->GetNumericElem(i) << std::endl;
    
  }

  if (m_verbose)
    std::cerr << "...adding the density column" << std::endl;
  
  // add the column
  Tag dtag("CA",label);
  
  // insert as a meta key
  m_header.addTag(dtag); 

  // add densities to the table
  m_table[dtag.GetName()] = dc;

  return 0;
}


