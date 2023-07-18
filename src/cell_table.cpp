#include "cell_table.h"
#include "csv.h"
#include "color_map.h"

#include <random>
#include <array>
#include <cstdlib>
#include <delaunator.hpp>

#ifdef HAVE_BOOST
#include <boost/functional/hash.hpp>
#endif

#ifdef HAVE_TIFFLIB
#include "tiff_writer.h"
#endif

#ifdef HAVE_CAIRO
#include "cairo/cairo.h"
#include "cairo/cairo-pdf.h"
#endif

#include "cell_utils.h"
#include "cell_graph.h"

#ifdef HAVE_HDF5
#include <H5Cpp.h>
#endif

// cgal is needed only for Voronoi diagram output
#ifdef HAVE_CGAL
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/draw_triangulation_2.h>
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Delaunay_triangulation_2<K> DelaunayData;
typedef K::Point_2 Point;
#endif



#ifdef HAVE_MLPACK
#include <mlpack/methods/gmm/gmm.hpp>
#include <mlpack/core.hpp>
#include <mlpack/methods/kmeans/kmeans.hpp>
#include <armadillo>
#endif

struct JPoint {
  
  float x;
  float y;

  JPoint(float mx, float my) : x(mx), y(my) {}
  
  std::string print() const { return std::to_string(x) + "," + std::to_string(y); }

  bool operator==(const JPoint& other) const {
    return x == other.x && y == other.y;
  }
  
};

namespace std {
    template<> struct hash<JPoint> {
        std::size_t operator()(const JPoint& p) const {
            std::size_t hx = std::hash<float>{}(p.x);
            std::size_t hy = std::hash<float>{}(p.y);

#ifdef HAVE_BOOST
            std::size_t seed = 0;
            boost::hash_combine(seed, hx);
            boost::hash_combine(seed, hy);
            return seed;
#else
            return hx ^ (hy << 1);  // Shift hy 1 bit to the left and XOR it with hx.
#endif
        }
    };
}

// hash table structures for Delaunay and Voronoi (to keep from duplicating lines)
struct pair_hash {
    template <class T1, class T2>
    std::size_t operator () (const std::pair<T1,T2> &pair) const {
        return std::hash<T1>()(pair.first) ^ std::hash<T2>()(pair.second);
    }
};

struct jline_eq {
    bool operator() (const std::pair<JPoint, JPoint>& lhs, const std::pair<JPoint, JPoint>& rhs) const {
        return (lhs.first == rhs.first && lhs.second == rhs.second) ||
               (lhs.first == rhs.second && lhs.second == rhs.first);
    }
};

bool CellTable::ContainsColumn(const std::string& name) const {
  return m_table.count(name) > 0;
}

void CellTable::HDF5Write(const std::string& file) const {

#ifdef HAVE_HDF5
  
  ///////
  // VAR (markers)
  ///////
  // setup for the "var" data (marker string names)
  std::vector<std::string> var_strings; // Use this to manage the lifetime of the strings
  std::vector<const char*> var;
  for (const auto& t : m_header.GetMarkerTags()) {
    var_strings.push_back(t.id);
    if (m_table.find(t.id) == m_table.end())
      std::cerr << " can't find " << t.id << std::endl;
  }
  
  // Then populate var using the strings in var_strings
  for (const auto& s : var_strings) {
    var.push_back(s.c_str());
  }

  /// DIMS
  hsize_t var_dims[1] = {var.size()};  
  hsize_t obs_dims[1] = {CellCount()};
  hsize_t x_dims[2] = {CellCount(), var_strings.size()}; // Number of rows and columns
  hsize_t empty_dims[1] = {0};
  
  ///////
  // OBS
  ///////
  // setup for the "var" data (marker string names)
  std::vector<std::string> obs_strings; // Use this to manage the lifetime of the strings
  std::vector<const char*> obs;
  obs_strings.push_back("x");
  obs_strings.push_back("y");
  obs_strings.push_back("cflag");
  obs_strings.push_back("pflag");  
  for (const auto& t : m_header.GetMetaTags()) {
    obs_strings.push_back(t.id);
    if (m_table.find(t.id) == m_table.end())
      std::cerr << " can't find " << t.id << std::endl;
  }
  
  // Then populate obs using the strings in obs_strings
  for (const auto& s : obs_strings) {
    obs.push_back(s.c_str());
  }

  ///////
  // OBS DATA
  ///////
  // Convert meta data to contiguous array
  std::vector<float> flat_data_meta;
  for (const auto& t : obs_strings) {
    const auto& ptr = m_table.find(std::string(t));
    assert(ptr != m_table.end());
    for (size_t i = 0; i < ptr->second->size(); i++) {
      flat_data_meta.push_back(ptr->second->GetNumericElem(i));
    }
  }

  
  ///////
  // X
  ///////
  // Convert X data to contiguous array
  // NOTE that in python, data is read as column major, so we need to do the same here
  std::vector<float> flat_data;
  for (const auto& t : var_strings) {
    const auto& ptr = m_table.find(std::string(t));
    assert(ptr != m_table.end());
    for (size_t i = 0; i < ptr->second->size(); i++) {
      flat_data.push_back(ptr->second->GetNumericElem(i));
    }
  }

  ///////
  // CELL ID
  ///////
  // setup for cellid
  std::vector<const char*> id_data;
  std::vector<std::unique_ptr<char[]>> unique_id_data;
  const auto& ptr = m_table.find("id");
  for (size_t i = 0; i < ptr->second->size(); i++) {
    std::string str = "cellid_" + std::to_string(static_cast<uint32_t>(ptr->second->GetNumericElem(i)));
    std::unique_ptr<char[]> cstr(new char[str.length() + 1]);
    //std::strcpy(cstr.get(), str.c_str());
    std::strncpy(cstr.get(), str.c_str(), str.size() + 1);
    id_data.push_back(cstr.get());
    unique_id_data.push_back(std::move(cstr));
  }

  //// FILE
  ///////
  // just in time setup file output
  H5::StrType strdatatype(H5::PredType::C_S1, H5T_VARIABLE);  
  H5::H5File h5file(file, H5F_ACC_TRUNC);

  // setup groups
  const H5std_string OBS_GROUP_NAME("obs");
  H5::Group group(h5file.createGroup(OBS_GROUP_NAME));

  const H5std_string VAR_GROUP_NAME("var");
  H5::Group vargroup(h5file.createGroup(VAR_GROUP_NAME));

  const H5std_string UNS_GROUP_NAME("uns");
  H5::Group unsgroup(h5file.createGroup(UNS_GROUP_NAME));

  //////
  // WRITE

  // write obs (float data)
  assert(obs_strings.size() == obs.size());
  H5::DataSpace obs_dataspace(1, obs_dims);
  for (const auto& o : obs_strings) {
    H5::DataSet obs_dataset(group.createDataSet(o, H5::PredType::NATIVE_FLOAT, obs_dataspace));
    std::vector<float> odata;
    const auto ptr = m_table.find(o);
    odata.reserve(CellCount());
    for (size_t i = 0; i < CellCount(); i++)
      odata.push_back(ptr->second->GetNumericElem(i));
    obs_dataset.write(odata.data(), H5::PredType::NATIVE_FLOAT);
  }

  // Assume that `obsGroup` is an H5::Group corresponding to the 'obs' group,
  // and `colName` is a std::string containing the name of your column.
  //H5::StrType strType(0, H5T_VARIABLE); // 0 means size of the string is variable
  hsize_t odim[1] = {obs.size()};
  H5::DataSpace obsnames_dataspace(1, odim);
  H5::Attribute columnOrderAttribute = group.createAttribute("column-order", strdatatype, obsnames_dataspace);
  columnOrderAttribute.write(strdatatype, obs.data());

  write_hdf5_dataframe_attributes(group);
  write_hdf5_dataframe_attributes(vargroup);  
  
  // write the _index data for obs
  H5::DataSet obsnames_dataset(group.createDataSet("_index", strdatatype, obs_dataspace));  
  obsnames_dataset.write(unique_id_data.data(), strdatatype);
  
  // write the _index data va
  H5::DataSpace var_dataspace(1, var_dims);
  H5::DataSet varnames_dataset(vargroup.createDataSet("_index", strdatatype, var_dataspace));
  varnames_dataset.write(var.data(), strdatatype);

  // and the var column order
  H5::DataSpace empty_dataspace = H5::DataSpace(1, empty_dims);
  H5::Attribute varcolumnOrderAttribute = vargroup.createAttribute("column-order", strdatatype, empty_dataspace);
  //varcolumnOrderAttribute.write(strdatatype, var.data());
  
  // write obs (index)
  //H5::DataSet obs_dataset(group.createDataSet("_index", strdatatype, obs_dataspace));
  //std::vector<const char*> cell_id;
  //obs_dataset.write(unique_id_data.data(), strdatatype);

  // write uns
  H5::DataSpace uns_dataspace(1, var_dims);
  H5::DataSet uns_dataset(unsgroup.createDataSet("all_markers", strdatatype, uns_dataspace));
  uns_dataset.write(var.data(), strdatatype);
  
  /*hsize_t obs_dims[1] = {obs.size()};
  H5::DataSpace obs_dataspace(1, obs_dims);
  H5::DataSet obs_dataset = h5file.createDataSet("/obs", strdatatype, obs_dataspace);  
  obs_dataset.write(obs_data.data(), strdatatype);

  // write obs (numeric)
  hsize_t obs_data_dims[2] = {CellCount(), obs.size()};
  H5::DataSpace obs_data_dataspace(1, obs_data_dims);
  H5::DataSet obs_data_dataset = h5file.createDataSet("/obs_data", H5::PredType::NATIVE_FLOAT, obs_dataspace);
  obs_data_dataset.write(flat_data_meta.data(), H5::PredType::NATIVE_FLOAT);
  */
  
  // write X (numeric)
  H5::DataSpace dataspace(2, x_dims);
  H5::DataSet dataset = h5file.createDataSet("/X", H5::PredType::NATIVE_FLOAT, dataspace);
  dataset.write(flat_data.data(), H5::PredType::NATIVE_FLOAT);

  h5file.close();
  
#else
  std::cerr << "Warning: Not able to read/write HDF5 file without including / linking HD5 library with build (and need -DHAVE_HDF5)" << std::endl;
#endif
  
}

std::ostream& operator<<(std::ostream& os, const CellTable& table) {

  os << "CellID -- " << table.m_table.at("id")->toString() << std::endl;
  os << "Pheno Flag -- "   << table.m_table.at("pflag")->toString() << std::endl;
  os << "Cell Flag -- "   << table.m_table.at("cflag")->toString() << std::endl;  
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

void CellTable::add_cell_to_table(const Cell& cell, bool nodata, bool nograph) {

  // add fixed data
  static_cast<IntCol*>(m_table["id"].get())->PushElem(cell.m_id);
  static_cast<IntCol*>(m_table["pflag"].get())->PushElem(cell.m_pheno_flag);
  static_cast<IntCol*>(m_table["cflag"].get())->PushElem(cell.m_cell_flag);  
  static_cast<FloatCol*>(m_table["x"].get())->PushElem(cell.m_x);
  static_cast<FloatCol*>(m_table["y"].get())->PushElem(cell.m_y);

  // add the info data
  if (!nodata) {
    size_t i = 0;
    for (const auto& t : m_header.GetDataTags()) {
      static_cast<FloatCol*>(m_table[t.id].get())->PushElem(cell.m_cols.at(i));
      i++;
    }
  }
  
  // add the graph data
  if (!nograph) {
    CellNode node(cell.m_spatial_ids, cell.m_spatial_dist);
    static_cast<GraphColumn*>(m_table["spat"].get())->PushElem(node);
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

void CellTable::SetupOutputWriter(const std::string& file) {

  // set the output to file or stdout
  if (file == "-") {
    m_archive = std::make_unique<cereal::PortableBinaryOutputArchive>(std::cout);
  } else {
    m_os = std::make_unique<std::ofstream>(file, std::ios::binary);
    m_archive = std::make_unique<cereal::PortableBinaryOutputArchive>(*m_os);
  }
  
}

#ifdef HAVE_TIFFLIB
void CellTable::Convolve(TiffWriter* otif, int boxwidth, float microns_per_pixel) {

  int xwidth, ywidth;
  assert(otif);
  TIFFGetField(otif->get(), TIFFTAG_IMAGEWIDTH, &xwidth);
  TIFFGetField(otif->get(), TIFFTAG_IMAGELENGTH, &ywidth);  
  
  TiffImage tif(xwidth, ywidth, 16);
  
  // Store points in an unordered set for fast lookup.
  FloatColPtr x_ptr = dynamic_pointer_cast<FloatCol>(m_table.at("x"));
  FloatColPtr y_ptr = dynamic_pointer_cast<FloatCol>(m_table.at("y"));  
  
  if (m_verbose)
    std::cerr << "...generating the integral sum image. Cell table has " << AddCommas(x_ptr->size()) << " cells" << std::endl;
  vector<vector<int>> dp(xwidth+1, vector<int>(ywidth+1, 0));  
  for (size_t i = 0; i < x_ptr->size(); i++) {
    int xp = static_cast<int>(x_ptr->GetNumericElem(i) / microns_per_pixel);
    int yp = static_cast<int>(y_ptr->GetNumericElem(i) / microns_per_pixel);
    if (xp > xwidth) {
      std::cerr << "calculated pixel xpos " << xp << " out of bounds of image width " << xwidth << " with microns_per_pixel " << microns_per_pixel << std::endl;
      assert(false);
    }
    if (yp > ywidth) {
      std::cerr << "calculated pixel ypos " << yp << " out of bounds of image width " << ywidth << " with microns_per_piyel " << microns_per_pixel << std::endl;
      assert(false);
    }
    
    dp[xp+1][yp+1] = 1;    
  }

  // Create a prefix sum matrix
  for (int i = 1; i <= xwidth; i++) {
    for (int j = 1; j <= ywidth; j++) {
      dp[i][j] += dp[i-1][j] + dp[i][j-1] - dp[i-1][j-1];
    }
  }

  // perform the convolution
  int bw2 = boxwidth / 2;
#pragma omp parallel for num_threads(m_threads)
  for (int i = 0; i < xwidth; i++) {
    if (m_verbose && i % 10000 == 0)
      std::cerr << "...convolving on image row " << AddCommas(i) << " of " << AddCommas(xwidth) << std::endl;
    for (int j = 0; j < ywidth; j++) {
      int x1 = max(i - bw2, 0);
      int y1 = max(j - bw2, 0);
      int x2 = min(i + bw2, xwidth-1);
      int y2 = min(j + bw2, ywidth-1);
      int count = dp[x2+1][y2+1] - dp[x1][y2+1] - dp[x2+1][y1] + dp[x1][y1];
      if (count > 65536)
	std::cerr << " warning count overflow of 16 bit range: " << count << std::endl;
      tif.SetPixel(i, j, PIXEL_GRAY, static_cast<uint16_t>(count));
    }
  }
  
  otif->Write(tif);

  return;
}
#endif

void CellTable::sortxy(bool reverse) {

  // fill the coordinate vector
  FloatColPtr x_ptr = dynamic_pointer_cast<FloatCol>(m_table.at("x"));
  FloatColPtr y_ptr = dynamic_pointer_cast<FloatCol>(m_table.at("y"));  
  
  std::vector<size_t> indices(x_ptr->size()); // index vector
  std::iota(indices.begin(), indices.end(), 0); // fill with 0, 1, ..., n-1

  if (m_verbose)
    std::cerr << "...sorting " << AddCommas(CellCount()) << " cells" << std::endl;
  
  // sort indices based on comparing values in (x,y)
  if (reverse) {
    std::sort(indices.begin(), indices.end(),
	      [&x_ptr, &y_ptr](size_t i1, size_t i2) {
		return std::hypot(x_ptr->GetNumericElem(i1), y_ptr->GetNumericElem(i1)) > 
		  std::hypot(x_ptr->GetNumericElem(i2), y_ptr->GetNumericElem(i2));
	      });
  } else {
    std::sort(indices.begin(), indices.end(),
	      [&x_ptr, &y_ptr](size_t i1, size_t i2) {
		return std::hypot(x_ptr->GetNumericElem(i1), y_ptr->GetNumericElem(i1)) <
		  std::hypot(x_ptr->GetNumericElem(i2), y_ptr->GetNumericElem(i2));
	      });
  }

  if (m_verbose)
    std::cerr << "...done sorting" << std::endl;
  
  // Now you can use this 'indices' vector to rearrange your 'another_column' vector
  for (auto t : m_table) {
     if (t.second->size())
      t.second->Order(indices);
   }
  
}

void CellTable::sort(const std::string& field, bool reverse) {

  auto it = m_table.find(field);
  if (it == m_table.end()) {
    std::cerr << "Warning: Sort field " << field << " not available" << std::endl;
    return;
  }

  if (it->second->GetType() != ColumnType::INT && 
      it->second->GetType() != ColumnType::FLOAT) {
    std::cerr << "Warning: Currently only able to sort on numeric fields: " << field << std::endl;
    return;
  }
    
  // fill the coordinate vector
  FloatColPtr c_ptr = dynamic_pointer_cast<FloatCol>(it->second);
  
  std::vector<size_t> indices(c_ptr->size()); // index vector
  std::iota(indices.begin(), indices.end(), 0); // fill with 0, 1, ..., n-1

  if (m_verbose)
    std::cerr << "...sorting " << AddCommas(CellCount()) << " cells" << std::endl;
  
  // sort indices based on comparing values in (x,y)
  if (reverse) {
    std::sort(indices.begin(), indices.end(),
	      [&c_ptr](size_t i1, size_t i2) {
		return c_ptr->GetNumericElem(i1) > c_ptr->GetNumericElem(i2);
	      });
  } else {
    std::sort(indices.begin(), indices.end(),
	      [&c_ptr](size_t i1, size_t i2) {
		return c_ptr->GetNumericElem(i1) < c_ptr->GetNumericElem(i2);
	      });
  }
  
  if (m_verbose)
    std::cerr << "...done sorting" << std::endl;
  
  // Now you can use this 'indices' vector to rearrange your 'another_column' vector
  for (auto t : m_table) {
     if (t.second->size())
      t.second->Order(indices);
   }
}

void CellTable::Delaunay(const std::string& pdf_delaunay,
		const std::string& pdf_voronoi,
	        int limit) {

#ifndef HAVE_BOOST
  std::cerr << "Warning: Will run but may be signifiacntly slower without including Boost library (uses Boost hash function for speed)" << std::endl;
#endif
  
  // fill the coordinate vector
  FloatColPtr x_ptr = dynamic_pointer_cast<FloatCol>(m_table.at("x"));
  FloatColPtr y_ptr = dynamic_pointer_cast<FloatCol>(m_table.at("y"));  

  float xmax = 0;
  float ymax = 0;  

  size_t ncells = CellCount();
  std::vector<double> coords(ncells * 2);
  
  for (size_t i = 0; i < ncells; i++) {
    coords[i*2  ] = x_ptr->GetNumericElem(i);
    coords[i*2+1] = y_ptr->GetNumericElem(i);
    if (coords[i*2] > xmax)
      xmax = coords[i*2];
    if (coords[i*2+1] > ymax)
      ymax = coords[i*2+1];
  }

#ifdef HAVE_CAIRO
  int width = xmax;
  int height = ymax;
#endif
  
  // construct the graph
  if (m_verbose) 
    std::cerr << "...constructing Delaunay triangulation on " << AddCommas(ncells) << " cells" << std::endl;
  delaunator::Delaunator d(coords);

  float micron_per_pixel = 0.325f;
  int limit_sq = limit <= 0 ? INT_MAX : limit * limit;
  
  // hash set to check if point already made
  std::unordered_set<std::pair<JPoint, JPoint>, pair_hash, jline_eq> lines;
  
  // reserve memory to avoid dynamic reallocation
  lines.reserve(d.triangles.size() / 3);
  
  size_t skip_count = 0;
  size_t draw_count = 0;
  
  // Find which lines to keep, by de-duplicating and removing lines > size limit
  for(size_t i = 0; i < d.triangles.size(); i+=3) {
    
    if (i/3 % 100000 == 0 && m_verbose)
      std::cerr << "...drawing triangle " << AddCommas(i/3) << " of " << AddCommas(d.triangles.size() / 3) << " for " <<
	AddCommas(ncells) << " cells" << std::endl;
    
    float x0 = static_cast<float>(d.coords[2 * d.triangles[i  ]]);
    float y0 = static_cast<float>(d.coords[2 * d.triangles[i  ] + 1]);
    float x1 = static_cast<float>(d.coords[2 * d.triangles[i+1]]);
    float y1 = static_cast<float>(d.coords[2 * d.triangles[i+1] + 1]);
    float x2 = static_cast<float>(d.coords[2 * d.triangles[i+2]]);
    float y2 = static_cast<float>(d.coords[2 * d.triangles[i+2] + 1]);

    std::array<std::pair<JPoint, JPoint>, 3> line_array = {      
      {std::make_pair(JPoint(x0, y0), JPoint(x1, y1)), 
       std::make_pair(JPoint(x1, y1), JPoint(x2, y2)), 
       std::make_pair(JPoint(x2, y2), JPoint(x0, y0))}
    };

    // deduplicate and remove Delaunay connections that are too short
    for (auto& line : line_array) {
      float dx = line.first.x - line.second.x;
      float dy = line.first.y - line.second.y;
      if (lines.find(line) == lines.end() && dx*dx + dy*dy <= limit_sq) {
        lines.insert(line);
        ++draw_count;
      } else {
        ++skip_count;
      }
    }
  }
  
  // store the adjaceny list
  std::unordered_map<JPoint, std::vector<JPoint>> adjList;  
  for (const auto& line : lines) {
    adjList[line.first].push_back(line.second);
    adjList[line.second].push_back(line.first); // assuming undirected graph
  }
  
  // build the adjaceny map
  std::unordered_set<JPoint> visited;
  std::unordered_map<JPoint, int> pointToComponentId;
  int currentComponentId = 1;
  
  // Iterative DFS using a stack
  std::stack<JPoint> stack;
  for (const auto& point : adjList) {
    if (visited.find(point.first) == visited.end()) {
      stack.push(point.first);
      while (!stack.empty()) {
	JPoint currentPoint = stack.top();
	stack.pop();
	if (visited.find(currentPoint) == visited.end()) {
	  visited.insert(currentPoint);
	  pointToComponentId[currentPoint] = currentComponentId;
	  for (const auto& neighbor : adjList[currentPoint]) {
	    stack.push(neighbor);
	  }
	}
      }
      currentComponentId++;
    }
  }
  
  // setup columns to store the delaunay components
  std::shared_ptr<IntCol> d_label = std::make_shared<IntCol>();
  std::shared_ptr<IntCol> d_size = std::make_shared<IntCol>();
  
  // count the number of nodes for each component
  std::unordered_map<int, size_t> dcount;
  for (const auto& c : pointToComponentId) {
    dcount[c.second]++;
  }
  
  // fill the data into the columns
  for (size_t i = 0; i < x_ptr->size(); i++) {
    //Point p = {x_ptr->GetNumericElem(i),y_ptr->GetNumericElem(i)};
    JPoint p = {x_ptr->GetNumericElem(i),y_ptr->GetNumericElem(i)};    
    
    // find the point in the component labels
    auto idd = pointToComponentId.find(p);

    // point not found, so delaunay edges already deleted
    if (idd == pointToComponentId.end()) {
      d_label->PushElem(0);
      d_size->PushElem(1);

    // point found
    } else {
      // set the label
      d_label->PushElem(idd->second);
      
      // set the label count
      assert(dcount.find(idd->second) != dcount.end());
      d_size->PushElem(dcount[idd->second]);
    }
  }
  
  // form the data tag
  Tag dtag_label(Tag::CA_TAG, "delaunay_component", "");
  AddColumn(dtag_label, d_label);
  Tag dtag_size(Tag::CA_TAG, "delaunay_count", "");
  AddColumn(dtag_size, d_size);
  
  // draw the Delaunay PDF
  if (!pdf_delaunay.empty()) {

    if (m_verbose)
      std::cerr << "...setting up for outputing Delaunay triangulation to PDF: " << pdf_delaunay << std::endl;

#ifdef HAVE_CAIRO
    
    cairo_surface_t *surface = cairo_pdf_surface_create(pdf_delaunay.c_str(), width, height);
    cairo_t *cr = cairo_create(surface);
    
    // Set the color of the lines to black.
    cairo_set_source_rgba(cr, 0, 0, 0, 1);
    cairo_set_line_width(cr, 0.3);

    // draw the actual lines
    for (const auto& line : lines) {
      //cairo_move_to(cr, line.first.x(), line.first.y());
      //      cairo_line_to(cr, line.second.x(), line.second.y());
      cairo_move_to(cr, line.first.x, line.first.y);
      cairo_line_to(cr, line.second.x, line.second.y);
      
    }

    cairo_stroke(cr); // Stroke all the lines
    
    cairo_new_path(cr); // Start a new path

    // setup a color map
    std::unordered_map<int, Color> color_map;
    for (auto c : pointToComponentId) {
      if (color_map.find(c.second) == color_map.end())
	color_map[c.second] = {rand() % 256,
			       rand() % 256,
			       rand() % 256};
    }

    // draw the points (cells) colored by component
    for (const auto& pair : pointToComponentId) {
      const JPoint& point = pair.first;
      int componentId = pair.second;
      Color c = color_map[componentId];

      cairo_set_source_rgb(cr, c.red/255.0, c.green/255.0, c.blue/255.0);
      cairo_arc(cr, point.x, point.y, 1.0, 0.0, 2.0 * M_PI);
      cairo_fill(cr);
    }
    
    if (m_verbose) 
      std::cerr << "...finalizing PDF of Delaunay triangulation on " << AddCommas(draw_count) <<
	" unique lines, skipping " << AddCommas(skip_count) << " lines for having length < " << limit << std::endl;
    
    // Clean up
    cairo_destroy(cr);
    cairo_surface_destroy(surface);
#else
    std::cerr << "Warning: Cairo PDF library needs to be linked during cysift build, no PDF will be output." << std::endl;
#endif
  }

  ///////////////
  // VORONOI
  ///////////////
  if (!pdf_voronoi.empty()) {

    if (m_verbose)
      std::cerr << "...setting up for outputing Voronoi diagram to PDF: " << pdf_voronoi << std::endl;

#if defined(HAVE_CAIRO) && defined(HAVE_CGAL)
    
    std::vector<Point> points(ncells);
    for (size_t i = 0; i < ncells; ++i) {
      points[i] = Point(x_ptr->GetNumericElem(i), y_ptr->GetNumericElem(i));
    }
    
    DelaunayData dt;
    dt.insert(points.begin(), points.end());

    // Create the Cairo surface and context
    cairo_surface_t *surface = cairo_pdf_surface_create(pdf_voronoi.c_str(), width, height); 
    cairo_t *cr = cairo_create(surface);

    // draw the original points
    for (const auto& p : points) {
      cairo_set_source_rgba(cr, 0.1, 0.0, 0.0, 1.0); // Set color to red
      cairo_arc(cr, p.x(), p.y(), 2.0, 0.0, 2.0 * M_PI);
      cairo_fill(cr);
    }

    // Set the color of the lines to black.
    cairo_set_source_rgba(cr, 0, 0, 0, 1);
    cairo_set_line_width(cr, 0.1);
    
    // draw the Voronoi edges    
    for (DelaunayData::Finite_edges_iterator eit = dt.finite_edges_begin(); eit != dt.finite_edges_end(); ++eit) {
      
      // The line segment between the circumcenters of the two triangles adjacent to each edge is an edge in the Voronoi diagram.
      DelaunayData::Face_handle f1 = eit->first;
      int adjacent_index = dt.mirror_index(f1, eit->second);
      DelaunayData::Face_handle f2 = f1->neighbor(adjacent_index);
      Point voronoi_edge_start = dt.circumcenter(f1);
      Point voronoi_edge_end = dt.circumcenter(f2);
      
      // Draw the Voronoi edge
      //cairo_move_to(cr, voronoi_edge_start.x(), voronoi_edge_start.y());
      //cairo_line_to(cr, voronoi_edge_end.x(), voronoi_edge_end.y());
    }

    //cairo_stroke(cr);
    
    // Finish the PDF
    if (m_verbose) 
      std::cerr << "...finalizing PDF of Voronoi diagram" << std::endl;
    
    cairo_show_page(cr);
    cairo_destroy(cr);
    cairo_surface_destroy(surface);
#else
    std::cerr << "Warning: Cairo PDF library and CGAL geometry libraries needs to be included/linked during cysift build for Voronoi output. Otherwise no PDF will be output." << std::endl;
#endif    
  }

  
  return;
  
}

void CellTable::OutputTable() {

  assert(m_archive);

  m_header.SortTags();
  
  // archive the header
  (*m_archive)(m_header);

  // create the cells and print
  size_t numRows = CellCount();

  auto id_ptr = m_table.at("id");
  auto cflag_ptr = m_table.at("cflag");
  auto pflag_ptr = m_table.at("pflag");  
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
    cell.m_cell_flag = static_cast<IntCol*>(cflag_ptr.get())->GetNumericElem(i);
    cell.m_pheno_flag = static_cast<IntCol*>(pflag_ptr.get())->GetNumericElem(i);
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

bool CellTable::HasColumn(const std::string& col) const {
  return m_table.find(col) != m_table.end();
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
  m_table["pflag"] = std::make_shared<IntCol>();
  m_table["cflag"] = std::make_shared<IntCol>();    
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


int CellTable::StreamTable(CellProcessor& proc, const std::string& file) {

  std::istream *inputStream = nullptr;
  std::unique_ptr<std::ifstream> fileStream;
  
  // set input from file or stdin
  if (file == "-") {
    inputStream = &std::cin;
  } else {
    fileStream = std::make_unique<std::ifstream>(file, std::ios::binary);
    if (!fileStream->good()) {
      std::cerr << "Error opening: " << file<< " - file may not exist" << std::endl;
      return 1;
    }
    inputStream = fileStream.get();
  }

  cereal::PortableBinaryInputArchive inputArchive(*inputStream);

  // First read the CellHeader
  try {
    inputArchive(m_header);
  } catch (const std::bad_alloc& e) {
    // Handle bad_alloc exception
    std::cerr << "Memory allocation failed during deserialization: " << e.what() << std::endl;
    return 1;  // or handle the error appropriately for your program
  } catch (const cereal::Exception& e) {
    // Handle exception if any error occurs while deserializing header
    std::cerr << "Error while deserializing header: " << e.what() << std::endl;
    return 1;  // or handle the error appropriately for your program
  }

  // process the header.
  // if , don't print rest
  int val = proc.ProcessHeader(m_header);
  if (val == CellProcessor::ONLY_WRITE_HEADER) { // just exit, all we need is header
    return 0;
  } else if (val == CellProcessor::SAVE_HEADER) { // we want to build the table
    initialize_cols();
  } else {
    ;  // the processor already handled the header, nothing to do here
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
    if (val == CellProcessor::WRITE_CELL) {
      proc.OutputLine(cell);
    } else if (val == CellProcessor::SAVE_CELL) {
      add_cell_to_table(cell, false, false);
    } else if (val == CellProcessor::SAVE_NODATA_CELL) {
      add_cell_to_table(cell, true, false);
    } else if (val == CellProcessor::SAVE_NODATA_NOGRAPH_CELL) {
      add_cell_to_table(cell, true, true);
    } else if (val == CellProcessor::NO_WRITE_CELL) {
      ; // do nothing
    } else {
      assert(false);
    }

  }

  return 0;
  
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


void CellTable::setCmd(const std::string& cmd) {

  m_cmd = cmd;
  m_header.addTag(Tag(Tag::PG_TAG, "", cmd));
    
}


void CellTable::GMM_EM() {

#ifdef HAVE_MLPACK
  arma::mat dataset;
  mlpack::data::Load("data.csv", dataset, true);

  // Initialize with the default arguments.
  mlpack::GMM gmm(dataset.n_rows, 3); // 3 is the number of Gaussians in the model. Adjust as necessary.
  
  // Train the model.
  gmm.Train(dataset);
#else
  std::cerr << "Warning: Unable to run GMM without linking / including MLPACK and armadillo in build" << std::endl;
#endif 
}

