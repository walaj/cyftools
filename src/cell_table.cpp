#include <random>
#include <cstdlib>

#include "cell_utils.h"
#include "cell_table.h"
#include "cell_graph.h"

#define CELL_BUFFER_LIMIT 1
static std::vector<Cell> cell_buffer;
static size_t debugr = 0;

#pragma omp threadprivate(cell_buffer)

const CellHeader& CellTable::GetHeader() const {
  return m_header;
}

bool CellTable::ContainsColumn(const std::string& name) const {
  return m_table.count(name) > 0;
}

std::ostream& operator<<(std::ostream& os, const CellTable& table) {

  os << "CellID -- " << table.m_table.at("id")->toString() << std::endl;
  os << "Flag -- "   << table.m_table.at("flag")->toString() << std::endl;
  os << "X -- "      << table.m_table.at("x")->toString() << std::endl;
  os << "Y -- "      << table.m_table.at("y")->toString() << std::endl;      
  
  for (const auto& t : table.m_header.GetDataTags()) {
    
    auto col_ptr = table.m_table.at(t.id);
    std::string ctype;

    switch (t.type) {
      case Tag::MA_TAG: ctype = "Marker"; break;
      case Tag::GA_TAG: ctype = "Graph"; break;
      case Tag::CA_TAG: ctype = "Meta"; break;
      default: ctype = "UNKNOWN"; break;
    }
    
    os << t.id << " -- " << ctype << " -- " << col_ptr->toString() << std::endl;
  }
  return os;
}

void CellTable::add_cell_to_table(const Cell& cell) {

  // add fixed data
  static_cast<IntCol*>(m_table["id"].get())->PushElem(cell.m_id);
  static_cast<FlagColumn*>(m_table["flag"].get())->PushElem(CellFlag(cell.m_flag));
  static_cast<FloatCol*>(m_table["x"].get())->PushElem(cell.m_x);
  static_cast<FloatCol*>(m_table["y"].get())->PushElem(cell.m_y);

  // add the info data
  size_t i = 0;
  for (const auto& t : m_header.GetDataTags()) {
    static_cast<FloatCol*>(m_table[t.id].get())->PushElem(cell.m_cols.at(i));
    i++;
  }

  // add the graph data
  CellNode node(cell.m_spatial_ids, cell.m_spatial_dist);
  static_cast<GraphColumn*>(m_table["spat"].get())->PushElem(node);
  
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
	c.first << " prev size " << prev_size.value() <<
	" current_size " << current_size << std::endl;
    }
    
    prev_size = current_size;
    n = std::max(n, current_size);
  }
  
  return n;
}

void CellTable::PrintTable(bool header) const {

  //  for (const auto& m : m_table) {
  // std::cerr << m.first << " - " << m.second->size() << std::endl;
  //}
  
  if (m_verbose)
    std::cerr << "...outputting table" << std::endl;
  
  if (header)
    m_header.Print();

  // graph data
  auto g_ptr = m_table.find("spat");
  
  // Create a lookup table with pointers to the Column values
  std::vector<ColPtr> lookup_table;
  for (const auto& c : m_header.GetDataTags()) {
    auto it = m_table.find(c.id);
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

    // print the graph
    if (g_ptr != m_table.end() && g_ptr->second->size()) {
      std::cerr << ",";
      g_ptr->second->PrintElem(row);
    }
    
    std::cout << std::endl;
  }

}

void CellTable::SetupOutputWriter(const std::string& file) {

  // set the output to file or stdout
  if (file == "-") {
    m_archive = std::make_unique<cereal::PortableBinaryOutputArchive>(std::cout);
  } else {
    m_os = std::make_unique<std::ofstream>(file, std::ios::binary);
    m_archive = std::make_unique<cereal::PortableBinaryOutputArchive>(*m_os);
  }
  
}

void CellTable::OutputTable() const {

  assert(m_archive);
  
  // archive the header
  (*m_archive)(m_header);

  // create the cells and print
  size_t numRows = CellCount();

  auto id_ptr = m_table.at("id");
  auto flag_ptr = m_table.at("flag");
  auto x_ptr = m_table.at("x");
  auto y_ptr = m_table.at("y");
  auto g_ptr = m_table.find("spat");

  std::vector<ColPtr> col_ptr;
  for (const auto& t : m_header.GetDataTags()) {
    col_ptr.push_back(m_table.at(t.id));
  }

  for (size_t i = 0; i < numRows; i++) {

    Cell cell;
    
    cell.m_id   = static_cast<IntCol*>(id_ptr.get())->GetNumericElem(i);
    cell.m_flag = static_cast<IntCol*>(flag_ptr.get())->GetNumericElem(i);
    cell.m_x    = static_cast<FloatCol*>(x_ptr.get())->GetNumericElem(i);
    cell.m_y    = static_cast<FloatCol*>(y_ptr.get())->GetNumericElem(i);

    for (const auto& c : col_ptr) {
      cell.m_cols.push_back(static_cast<FloatCol*>(c.get())->GetNumericElem(i));
    }
    
    // fill the Cell graph
    if (g_ptr != m_table.end()) {
      const CellNode& n = static_cast<GraphColumn*>(g_ptr->second.get())->GetNode(i);
      n.FillSparseFormat(cell.m_spatial_ids, cell.m_spatial_dist);
    }
    
    // write it
    (*m_archive)(cell);
    
  }

}

void CellTable::Crop(float xlo, float xhi, float ylo, float yhi) {
  
  // Find the x and y columns in the table
  auto x_it = m_table.find("x");
  auto y_it = m_table.find("y");

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
  for (const auto& t : m_header.GetDataTags()) {
    if (t.type == Tag::MA_TAG) {
      auto ptr = m_table.find(t.id);
      assert(ptr != m_table.end());
      data.push_back({ptr->first, ptr->second});
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
  for (const auto& t : m_header.GetDataTags()) {
    if (t.type == Tag::MA_TAG) {
      if (t.id.length() > spacing)
	spacing = t.id.length() + 2;
    }
  }
  
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
  auto x_ptr = m_table.at("x");
  auto y_ptr = m_table.at("y");
  
  float x_min = x_ptr->Min();
  float x_max = x_ptr->Max();
  float y_min = y_ptr->Min();
  float y_max = y_ptr->Max();
  
  // scale it to fit on the plot
  size_t nc = CellCount();
  std::vector<std::pair<int, int>> scaled_coords;
  for (size_t i = 0; i < nc; i++) {
    float x_ = x_ptr->GetNumericElem(i);
    float y_ = y_ptr->GetNumericElem(i);    
    
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

void CellTable::AddColumn(const Tag& tag,
			  ColPtr value) {
  
  if (value->size() != CellCount()) {
    throw std::runtime_error("Adding column of incorrect size");
  }
  
  // check if already exists in the table
  if (m_table.find(tag.id) != m_table.end()) {
    //throw std::runtime_error("Adding column already in table");
    std::cerr << "Warning: Overwriting existing column " << tag.id << std::endl;
  }
  
  // insert as a meta key
  m_header.addTag(tag);

  // add to the table
  m_table[tag.id] = value;
  
}

/*void CellTable::AddMetaColumn(const std::string& key,
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
  
  }*/

void CellTable::SubsetROI(const std::vector<Polygon> &polygons) {

  size_t nc = CellCount();
  
  // Initialize the "sparoi" column
  FloatColPtr new_data = std::make_shared<FloatCol>();
  new_data->reserve(nc);
  
  // Loop the table and check if the cell is in the ROI
  for (size_t i = 0; i < nc; i++) {
    float x_ = m_table.at("x")->GetNumericElem(i);
    float y_ = m_table.at("y")->GetNumericElem(i);    

    // Loop through all polygons and check if the point is inside any of them
    for (const auto &polygon : polygons) {
      // if point is in this polygon, add the polygon id number to the roi
      if (polygon.PointIn(x_,y_)) {
	new_data->SetNumericElem(static_cast<float>(polygon.Id), i);
        // uncomment below if want to prevent over-writing existing
        // break;
      }
    }
  }
}
 
void CellTable::KNN_spatial(int num_neighbors, int dist) {
  
  // number of cells
  int nobs = CellCount();
  
  if (m_verbose)
    std::cerr << "...finding K nearest-neighbors (spatial) on " << AddCommas(nobs) << " cells" << std::endl;  
   
   // column major the coordinate data
   std::vector<float> concatenated_data;
   
   const int ndim = 2;

   auto x_ptr = m_table.at("x");
   auto y_ptr = m_table.at("y");
   
   const auto& x_data = std::dynamic_pointer_cast<FloatCol>(x_ptr)->getData(); // just a ptr, not a copy
   if (m_verbose) 
     std::cerr << "...adding " << AddCommas(x_data.size()) << " points on x" << std::endl;
   concatenated_data.insert(concatenated_data.end(), x_data.begin(), x_data.end());
   
   const auto& y_data = std::dynamic_pointer_cast<FloatCol>(y_ptr)->getData();
   if (m_verbose) 
     std::cerr << "...adding " << AddCommas(y_data.size()) << " points on y" << std::endl;
   concatenated_data.insert(concatenated_data.end(), y_data.begin(), y_data.end());
   
   // convert to row major?
   column_to_row_major(concatenated_data, nobs, ndim);
   
   // get the cell id colums
   //auto id_ptr = GetIDColumn();
   auto id_ptr = m_table.at("id");
   
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
  //knncolle::VpTree<knncolle::distances::Euclidean, int, float> searcher(ndim, nobs, concatenated_data.data());
  //knncolle::AnnoyEuclidean<int, float> searcher(ndim, nobs, concatenated_data.data());
  knncolle::Kmknn<knncolle::distances::Euclidean, int, float> searcher(ndim, nobs, concatenated_data.data());  
    
  // archive the header
  assert(m_archive);
  (*m_archive)(m_header);

  // create the cells and print
  size_t numRows = CellCount();

  // setup for converting to Cell
  auto flag_ptr = m_table.at("flag");
  std::vector<ColPtr> col_ptr;
  for (const auto& t : m_header.GetDataTags()) {
    col_ptr.push_back(m_table.at(t.id));
  }
  
  // Create a thread-local storage for each thread to hold its cells.
  
#pragma omp parallel for num_threads(m_threads)
  for (size_t i = 0; i < nobs; ++i) {

    // verbose printing
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

    CellNode node;
    
    // remove less than distance
    if (dist > 0) {
      size_t osize = neigh.size();
      Neighbors neigh_trim;
      for (const auto& nnn : neigh) {
	if (nnn.second < dist)
	  neigh_trim.push_back(nnn);
      }

      // print a warning if we trimmed off too many neighbors
      if (osize == neigh_trim.size()) {  
    #pragma omp critical
	{
	  lost_cell++;
	  if (lost_cell % 500 == 0)
	    std::cerr << "osize " << osize << " Lost cell " << AddCommas(lost_cell) << " of " << AddCommas(nobs) << std::endl;
	}
      }
      
      neigh = neigh_trim;

    }

    // add the cell flags
    std::vector<uint64_t> flag_vec(neigh.size());
    for (size_t j = 0; j < neigh.size(); j++) {
      flag_vec[j] = static_cast<IntCol*>(flag_ptr.get())->GetNumericElem(neigh.at(j).first);
    }

    assert(flag_vec.size() == neigh.size());
    
    node = CellNode(neigh, flag_vec);
    
    // build the Cell object
    ////////////////////////
    Cell cell;
    cell.m_id   = static_cast<IntCol*>(id_ptr.get())->GetNumericElem(i);
    cell.m_flag = static_cast<IntCol*>(flag_ptr.get())->GetNumericElem(i);
    cell.m_x    = static_cast<FloatCol*>(x_ptr.get())->GetNumericElem(i);
    cell.m_y    = static_cast<FloatCol*>(y_ptr.get())->GetNumericElem(i);

    // fill the Cell data columns
    for (const auto& c : col_ptr) {
      cell.m_cols.push_back(static_cast<FloatCol*>(c.get())->GetNumericElem(i));
    }

    // fill the Cell graph data
    node.FillSparseFormat(cell.m_spatial_ids, cell.m_spatial_dist);
    assert(flag_vec.size() == cell.m_spatial_ids.size());
    cell.m_spatial_flags = flag_vec;

    // Add the cell to the thread-local buffer.
    //cell_buffer.push_back(cell);

    //if (cell_buffer.size() >= CELL_BUFFER_LIMIT) {
#pragma omp critical
    {
      	(*m_archive)(cell);
      //for (const auto& buffered_cell : cell_buffer) {
      //	(*m_archive)(buffered_cell);
      //}
    }
    //cell_buffer.clear();
    //}
    
  }// end for
  
  // After the loop, dump any remaining cells in the buffer.
#pragma omp critical
  {
    for (const auto& buffered_cell : cell_buffer) {
      (*m_archive)(buffered_cell);
    }
  }
  cell_buffer.clear();
  
  if (m_verbose)
    std::cerr << "...done with graph construction" << std::endl;
  
  return;
  
  Tag gtag(Tag::GA_TAG, "spat", "NN:" + std::to_string(num_neighbors));
  //gtag.addValue("NN", std::to_string(num_neighbors));
  AddColumn(gtag, graph);
}

void CellTable::KNN_marker(int num_neighbors) {

  return;

  /*
  // get the number of markers
  int ndim = 0;
  for (const auto& t : m_header.GetDataTags())
    if (t.type == Tag::MA_TAG)
      ndim++;

  // number of cells
  int nobs = CellCount();
  
  if (m_verbose)
    std::cerr << "...finding K nearest-neighbors on " << AddCommas(nobs) << " cells" << std::endl;  
  
  // column major the marker data
  std::vector<double> concatenated_data;

  // concatenate the data
  for (const auto& t : m_header.GetDataTags()) {

    // skip if not a marker tag
    if (t.type != Tag::MA_TAG)
      continue;
    
    auto it = m_table.find(t.id);
    if (it != m_table.end()) {
      auto column_ptr = it->second;
      auto numeric_column_ptr = std::dynamic_pointer_cast<FloatCol>(column_ptr);
      assert(numeric_column_ptr);
      
      // Check if the column is a NumericColumn
      //if (auto numeric_column_ptr = std::dynamic_pointer_cast<FloatCol>(column_ptr)) {
      const auto& data = numeric_column_ptr->getData();
      concatenated_data.insert(concatenated_data.end(), data.begin(), data.end());
      //     }
    }
  }
  
  // convert to row major?
  column_to_row_major(concatenated_data, nobs, ndim);

  // get the cell id colums
  auto id_ptr = m_table.at("id");  
  //auto id_ptr = GetIDColumn();
  
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
  Tag gtag(Tag::GA_TAG, "knn", "NN:" + std::to_string(num_neighbors));  
  AddColumn(gtag, graph);
  */
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

void CellTable::initialize_cols() {

  // initialize cell id
  m_table["id"] = std::make_shared<IntCol>();
  m_table["flag"] = std::make_shared<FlagColumn>();  
  m_table["x"] = std::make_shared<FloatCol>();
  m_table["y"] = std::make_shared<FloatCol>();

  m_table["spat"] = std::make_shared<GraphColumn>();
  
  // initialize other data
  for (const auto& t : m_header.GetDataTags()) {

    // check that already doesn't exist
    auto m_ptr = m_table.find(t.id);
    if (m_ptr != m_table.end()) {
      throw std::runtime_error("Can't initialize " + t.id + " since already in table");
    }

    // make the empty float column
    m_table[t.id] = std::make_shared<FloatCol>();
  }

  
}


/*void CellTable::select(uint64_t on, uint64_t off) {

  FlagColPtr fc = std::dynamic_pointer_cast<FlagColumn>(m_table["flag"]);
  assert(fc);
  
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

    int index = m_header.WhichColumn(b.first, Tag::MA_TAG);
        
    for (size_t i = 0; i < nc->size(); i++) {
      if (nc->GetNumericElem(i) >= b.second.first &&
	  nc->GetNumericElem(i) <= b.second.second) {
	fc->SetFlagOn(index, i);
      }
    }
  }
  
  //Tag ftag("FA","cellflag");
  //AddFlagColumn(ftag, fc, true);
  m_table["flag"] = fc;
}
*/


/*void CellTable::AddFlagColumn(const Tag& tag,
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
  

  }*/


/*IntColPtr CellTable::GetIDColumn() const {

  for (const auto& t: m_header.GetColTags()) {
    if (t.isIDTag()) {
      return std::dynamic_pointer_cast<IntCol>(m_table.at(t.GetName()));
    }
  }

  assert(false);
  return IntColPtr();
  }*/

void CellTable::StreamTable(CellProcessor& proc, const std::string& file) {

  bool build_table_memory =false;
  
  std::istream *inputStream = nullptr;
  std::unique_ptr<std::ifstream> fileStream;
  
  // set input from file or stdin
  if (file == "-") {
    inputStream = &std::cin;
  } else {
    fileStream = std::make_unique<std::ifstream>(file, std::ios::binary);
    inputStream = fileStream.get();
  }

  cereal::PortableBinaryInputArchive inputArchive(*inputStream);

  // First read the CellHeader
  try {
    inputArchive(m_header);
  } catch (const std::bad_alloc& e) {
    // Handle bad_alloc exception
    std::cerr << "Memory allocation failed during deserialization: " << e.what() << std::endl;
    return;  // or handle the error appropriately for your program
  } catch (const cereal::Exception& e) {
    // Handle exception if any error occurs while deserializing header
    std::cerr << "Error while deserializing header: " << e.what() << std::endl;
    return;  // or handle the error appropriately for your program
  }

  // process the header.
  // if , don't print rest
  int val = proc.ProcessHeader(m_header);
  if (val == 1) { // just exit, all we need is header
    return;
  } else if (val == 2) { // we want to build the table
    initialize_cols();
  }
  
  // now read the Cell objects
  while (true) {
    Cell cell;
    try {
      m_count++;            
      if (m_verbose && (m_count % 500000 == 0 || m_count == 1))
	std::cerr << "...reading cell " << AddCommas(m_count) << std::endl;
      inputArchive(cell);
    } catch (const cereal::Exception& e) {
      // Catch exception thrown on EOF or any other errors
      break;
    }

    // process the cell and output if needed (returns 1)
    int val = proc.ProcessLine(cell);
    if (val == 1) {
      proc.OutputLine(cell);
    } else if (val == 2) {
      build_table_memory = true;
      add_cell_to_table(cell);
    }

  }

  if (!build_table_memory)
    return;
  
  // if graph is entirely empty, just remove it
  // but it should at least exist at this point, with empty nodes
  shared_ptr<GraphColumn> gc = std::dynamic_pointer_cast<GraphColumn>(m_table["spat"]);
  bool empty = true;
  for (size_t i = 0; i < gc->size(); i++) {
    if (gc->GetNode(i).size()) {
      empty = false;
      break;
    }
  }
  if (empty)
    m_table.erase("spat");
  
  // clean up table
  for (auto it = m_table.begin(); it != m_table.end();) {
    if (!it->second->size()) {
      it = m_table.erase(it);
    } else {
      ++it;
    }
  }
  
  
}
  
void CellTable::StreamTableCSV(LineProcessor& proc, const std::string& file) {

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

  // get read to read file
  std::string line;
  bool header_read = false;
  size_t flag_index = -1;
  char* next_line_ptr;

  // for verbose
  size_t count = 0;
  
  // Read and process the lines
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

int CellTable::RadialDensity(std::vector<uint64_t> inner, std::vector<uint64_t> outer,
			     std::vector<uint64_t> logor, std::vector<uint64_t> logand,
			     std::vector<std::string> label) {

  // check the radial geometry parameters
  assert(inner.size());
  assert(inner.size() == outer.size());
  assert(inner.size() == logor.size());
  assert(inner.size() == logand.size());
  assert(inner.size() == label.size());    

  if (m_verbose)  {
    for (size_t i = 0; i < inner.size(); i++) {
      std::cerr << "In " << inner.at(i) << " Out " << outer.at(i) <<
	" log " << logor.at(i) << " logand " << logand.at(i) << " label " << label.at(i) << std::endl;
    }
  }
  
  // get the flag column
  shared_ptr<FlagColumn> fc = std::dynamic_pointer_cast<FlagColumn>(m_table["flag"]);
  assert(fc);
  assert(fc->size());

  // get the graph column
  shared_ptr<GraphColumn> gc = std::dynamic_pointer_cast<GraphColumn>(m_table["spat"]);
  assert(gc);
  assert(gc->size());

  // check sizes are good
  assert(gc->size() == fc->size());
  assert(gc->size() == CellCount());

  // check if already exists in the table
  for (const auto& l : label) {
    if (m_table.find(l) != m_table.end()) {
      throw std::runtime_error("Adding density meta column " + l + " already in table");
    }
  }
  
  
  // store the densities
  std::vector<std::shared_ptr<FloatCol>> dc(inner.size());
  for(auto& ptr : dc) {
    ptr = std::make_shared<FloatCol>();
    ptr->resize(gc->size());
  }

  //shared_ptr<FloatCol> dc = std::make_shared<FloatCol>();  
  //dc->resize(gc->size());
  //dc->SetPrecision(2);

  // get the cell id colums
  auto id_ptr = m_table.at("id"); 
  assert(id_ptr->size() == gc->size());

  // flip it so that the id points to the zero-based index
  std::unordered_map<uint32_t, size_t> inverse_lookup;
  if (m_verbose)
    std::cerr << "...done loading, doing the flip" << std::endl;
  uint32_t max_cell_id = 0;
  for (size_t i = 0; i < id_ptr->size(); ++i) {
    uint32_t value = id_ptr->GetNumericElem(i);
    inverse_lookup[value] = i;
    if (value > max_cell_id)
      max_cell_id = value;
  }

  // convert to vector for speed
  //std::cerr << " Max cell id " << AddCommas(max_cell_id) << std::endl;
  //std::vector<uint32_t> inverse_lookup_v(max_cell_id);
  //for (const auto& i : inverse_lookup)
  //  inverse_lookup_v[i.first] = i.second;
  
  if (m_verbose)
    std::cerr << "...radial density: starting loop" << std::endl;
  
  // loop the cells
  #pragma omp parallel for num_threads(m_threads)
  for (size_t i = 0; i < gc->size(); i++) {

    // initialize the counts for each radial condition
    std::vector<float> cell_count(inner.size());

    const CellNode& node = gc->GetNode(i);

    const Neighbors& neigh = node.get_neighbors();

    // loop the nodes connected to each cell
    for (const auto& n : neigh) {
      
      uint32_t cellindex = inverse_lookup[n.first];
      
      // test if the connected cell meets the flag criteria
      // n.first is cell_id of connected cell to this cell
      for (size_t j = 0; j < inner.size(); j++) {
	
	// both are 0, so take all cells OR it meets flag criteria
	if ( (!logor[j] && !logand[j]) || fc->TestFlagAndOr(logor[j], logand[j], cellindex)) {
	  
	  // then increment cell count if cell in bounds
	  cell_count[j] += n.second >= inner[j] && n.second <= outer[j];
	  
	}
      }
    }
    
    // calculate the density
    std::vector<float> area(inner.size());
    for (size_t j = 0; j < inner.size(); ++j) {
      float outerArea = static_cast<float>(outer[j]) * static_cast<float>(outer[j]) * 3.1415926535f;
      float innerArea = static_cast<float>(inner[j]) * static_cast<float>(inner[j]) * 3.1415926535f;
      area[j] = outerArea - innerArea;
    }

    // do the density calculation for each condition
    // remember, i is iterator over cells, j is over conditions
    for (size_t j = 0; j < area.size(); ++j) {
      if (!neigh.empty()) {
	float value = cell_count[j] * 1000000 / area[j]; // density per 1000 square pixels
	dc[j]->SetNumericElem(value, i);
      } else {
	dc[j]->SetNumericElem(0, i); 
      }
    }
    
    if (m_verbose && i % 50000 == 0)
      std::cerr << "...radial density: looping on " << AddCommas(i) << " and found density[0] " << dc.at(0)->GetNumericElem(i) << std::endl;
    
  } // end the main cell loop
  
  if (m_verbose)
    std::cerr << "...adding the density column" << std::endl;

  for (size_t i = 0; i < label.size(); ++i) {
    
    // form the data tag
    Tag dtag(Tag::CA_TAG, label[i], "");

    AddColumn(dtag, dc[i]);
    
    // insert column info 
    //m_header.addTag(dtag); 
    
    // add densities to the table
    //m_table[dtag.id] = dc[i];
  }

  /*  
  // add the column
  Tag dtag(Tag::CA_TAG, label, "");
  
  // insert as a meta key
  m_header.addTag(dtag); 

  // add densities to the table
  m_table[dtag.id] = dc;
  */
  return 0;
}
