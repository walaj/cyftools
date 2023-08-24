#include "cell_table.h"
#include "csv.h"
#include "color_map.h"
#include "cell_selector.h"
#include "cell_utils.h"
#include "cell_graph.h"

#include <random>
#include <cstdlib>
#include <delaunator.hpp>

#ifdef HAVE_MLPACK
#include <mlpack/core.hpp>
#include <mlpack/methods/dbscan/dbscan.hpp>
#endif

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

CellTable::CellTable(size_t num_cells) {

  // make the ID
  m_id_ptr = std::make_shared<IntCol>();
  m_id_ptr->resize(num_cells);

  // instantiate with ids
  for (size_t i = 0; i < num_cells; i++) {
    uint64_t mid = static_cast<uint64_t>(0) << 32 | (i+1); // 0 is sample id
    m_id_ptr->SetNumericElem(mid, i);
  }
  //m_table["id"] = id;
  
  // make the flag columns
  m_pflag_ptr = std::make_shared<IntCol>();
  m_cflag_ptr = std::make_shared<IntCol>();
  m_pflag_ptr->resize(num_cells);
  m_cflag_ptr->resize(num_cells);  

  // make the x and y columns
  m_x_ptr = std::make_shared<FloatCol>();
  m_y_ptr = std::make_shared<FloatCol>();
  m_x_ptr->resize(num_cells);
  m_y_ptr->resize(num_cells);  
}

/*void CellTable::HDF5Write(const std::string& file) const {

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
    const FloatColPtr fp = dynamic_pointer_cast<FloatCol>(ptr->second);
    for (size_t i = 0; i < fp->size(); i++) {
      flat_data_meta.push_back(fp->getData().at(i));
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
    const FloatColPtr fp = dynamic_pointer_cast<FloatCol>(ptr->second);    
    for (size_t i = 0; i < ptr->second->size(); i++) {
      flat_data.push_back(fp->getData().at(i));
    }
  }

  ///////
  // CELL ID
  ///////
  // setup for cellid
  std::vector<const char*> id_data;
  std::vector<std::unique_ptr<char[]>> unique_id_data;
  assert(m_table.find("id") != m_table.end());
  const IDColPtr ptr = dynamic_pointer_cast<IDCol>(m_table.at("id"));
  for (size_t i = 0; i < ptr->size(); i++) {
    std::string str = "cellid_" + std::to_string(static_cast<uint32_t>(ptr->getData().at(i) & 0xFFFFFFFF));
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
    assert(ptr != m_table.end());
    const FloatColPtr fp = dynamic_pointer_cast<FloatCol>(ptr->second);
    odata.reserve(CellCount());
    size_t nn = CellCount();
    for (size_t i = 0; i < nn; i++)
      odata.push_back(fp->getData().at(i));
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
  
  // hsize_t obs_dims[1] = {obs.size()};
  // H5::DataSpace obs_dataspace(1, obs_dims);
  // H5::DataSet obs_dataset = h5file.createDataSet("/obs", strdatatype, obs_dataspace);  
  // obs_dataset.write(obs_data.data(), strdatatype);

  // // write obs (numeric)
  // hsize_t obs_data_dims[2] = {CellCount(), obs.size()};
  // H5::DataSpace obs_data_dataspace(1, obs_data_dims);
  // H5::DataSet obs_data_dataset = h5file.createDataSet("/obs_data", H5::PredType::NATIVE_FLOAT, obs_dataspace);
  // obs_data_dataset.write(flat_data_meta.data(), H5::PredType::NATIVE_FLOAT);
  
  // write X (numeric)
  H5::DataSpace dataspace(2, x_dims);
  H5::DataSet dataset = h5file.createDataSet("/X", H5::PredType::NATIVE_FLOAT, dataspace);
  dataset.write(flat_data.data(), H5::PredType::NATIVE_FLOAT);

  h5file.close();
  
#else
  std::cerr << "Warning: Not able to read/write HDF5 file without including / linking HD5 library with build (and need -DHAVE_HDF5)" << std::endl;
#endif
  
}*/

std::ostream& operator<<(std::ostream& os, const CellTable& table) {

  // make the cell and sample ids
  //IntColPtr id_ptr = dynamic_pointer_cast<NumericColumn<uint64_t>>(table.m_table.at("id"));
  std::shared_ptr<NumericColumn<uint32_t>> sid = make_shared<NumericColumn<uint32_t>>();
  std::shared_ptr<NumericColumn<uint32_t>> cid = make_shared<NumericColumn<uint32_t>>();
  cid->reserve(table.size());
  sid->reserve(table.size());
  for (const auto& c : table.m_id_ptr->getData()) {
    sid->PushElem(static_cast<uint32_t>(c >> 32));
    cid->PushElem(static_cast<uint32_t>(c & 0xFFFFFFFF));
  }

  os << "CellID -- "   << cid->toString() << std::endl;
  os << "SampleID -- " << sid->toString() << std::endl;
  os << "RawID -- " << table.m_id_ptr->toString() << std::endl;
  os << "Pheno Flag -- "   << table.m_pflag_ptr->toString() << std::endl;
  os << "Cell Flag -- "    << table.m_cflag_ptr->toString() << std::endl;
  os << "X -- "      << table.m_x_ptr->toString() << std::endl;
  os << "Y -- "      << table.m_y_ptr->toString() << std::endl;
  
  for (const auto& t : table.m_header.GetDataTags()) {
    
    const FloatColPtr col_ptr = table.m_table.at(t.id);
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
  m_id_ptr->PushElem(cell.id);
  m_pflag_ptr->PushElem(cell.pflag);
  m_cflag_ptr->PushElem(cell.cflag);
  m_x_ptr->PushElem(cell.x);
  m_y_ptr->PushElem(cell.y);

  // add the info data
  if (!nodata) {
    size_t i = 0;
    for (const auto& t : m_header.GetDataTags()) {
      m_table[t.id]->PushElem(cell.cols.at(i));
      i++;
    }
  }
  
}

void CellTable::clusterDBSCAN(CellSelector select,
			      float epsilon,
			      size_t min_size,
			      size_t min_cluster_size
			      ) {

  validate();
  
#ifdef HAVE_MLPACK

  std::vector<double> xvec(m_x_ptr->getData().begin(),m_x_ptr->getData().end());
  std::vector<double> yvec(m_y_ptr->getData().begin(),m_y_ptr->getData().end());
  
  // Convert std::vector to Armadillo row vectors.
  arma::rowvec x_row(CellCount());
  arma::rowvec y_row(CellCount());

  size_t count = 0;
  std::vector<size_t> idx(CellCount());
  for (size_t i = 0; i < m_pflag_ptr->size(); i++) {
    // only select on certain cells
    if (select.TestFlags(m_pflag_ptr->getData().at(i),
			 m_cflag_ptr->getData().at(i))) {
      x_row(count) = m_x_ptr->getData().at(i);
      y_row(count) = m_y_ptr->getData().at(i);
      idx[count] = i;
      count++;
    }
  }
  x_row.resize(count);
  y_row.resize(count);
  idx.resize(count);
  
  // Combine the two row vectors into a matrix.
  arma::mat dataset(2, x_row.size());
  dataset.row(0) = x_row;
  dataset.row(1) = y_row;

  if (m_verbose)
    std::cerr << "...running dbscan (mlpack) on " << AddCommas(x_row.size()) <<
      " cells. Epsilon: " << epsilon << " Min points: " <<
      min_size << " min cluster size: " << min_cluster_size << std::endl;
  
  // Perform DBSCAN clustering.
  arma::Row<size_t> assignments(idx.size(), arma::fill::zeros); // to store cluster assignments
  
  // The parameters are: epsilon (radius of neighborhood), minimum points in neighborhood
  const bool batchmode = false;
  mlpack::DBSCAN<> dbscan(epsilon, min_size, batchmode);
  dbscan.Cluster(dataset, assignments);
  
  // histogram the data
  std::unordered_map<size_t, size_t> cluster_hist;
  for (const auto& a : assignments) {
    cluster_hist[a]++;
  }
  
  // store the clusters data
  FloatColPtr fc = std::make_shared<FloatCol>();
  fc->resize(CellCount()); // resize zeros as well
  for (const auto& cc : fc->getData())
    assert(cc == 0);
  
  // make the assignments
  std::unordered_set<size_t> assignment_set;
  assert(idx.size() == assignments.size());
  for (size_t i = 0; i < idx.size(); i++) {

    // don't add if not a big cluster
    if (cluster_hist[assignments[i]] < min_cluster_size) {
      fc->SetValueAt(idx[i], 0);
    // otherwise add the cluster to m_table
    } else {
      if (assignments[i] == SIZE_MAX) {
	fc->SetValueAt(idx[i], -1);
      } else {
	fc->SetValueAt(idx[i], assignments[i] + 1);	
	assignment_set.insert(assignments[i]); // just to keep track of number of assignments
      }
    }
    
  }
  
  if (m_verbose)
    std::cerr << "...produced " << AddCommas(assignment_set.size()) << " clusters" << std::endl;

  Tag ttag1(Tag::CA_TAG, "dbscan_cluster","");
  assert(fc->size() == CellCount());
  AddColumn(ttag1, fc);
  
#else  
  std::cerr << "...mlpack header library not available. Need to include this during compilation" << std::endl;
#endif  
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

  validate();
  size_t prev_size = m_x_ptr->size();
  size_t n = 0;

  // loop the table and find the maximum length
  for (const auto &c : m_table) {
    size_t current_size = c.second->size();

    if (current_size != prev_size) {
      std::cerr << "Warning: Column sizes do not match. Column: " <<
	c.first << " prev size " << prev_size <<
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

  validate();

  int xwidth, ywidth;
  assert(otif);
  TIFFGetField(otif->get(), TIFFTAG_IMAGEWIDTH, &xwidth);
  TIFFGetField(otif->get(), TIFFTAG_IMAGELENGTH, &ywidth);  
  
  TiffImage tif(xwidth, ywidth, 16);
  
  if (m_verbose)
    std::cerr << "...generating the integral sum image. Cell table has " << AddCommas(m_x_ptr->size()) << " cells" << std::endl;

  // Store points in an unordered set for fast lookup.  
  vector<vector<int>> dp(xwidth+1, vector<int>(ywidth+1, 0));  
  for (size_t i = 0; i < m_x_ptr->size(); i++) {
    int xp = static_cast<int>(m_x_ptr->GetNumericElem(i) / microns_per_pixel);
    int yp = static_cast<int>(m_y_ptr->GetNumericElem(i) / microns_per_pixel);
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

  validate();
  
  std::vector<size_t> indices(m_x_ptr->size()); // index vector
  std::iota(indices.begin(), indices.end(), 0); // fill with 0, 1, ..., n-1

  if (m_verbose)
    std::cerr << "...sorting " << AddCommas(CellCount()) << " cells" << std::endl;
  
  // sort indices based on comparing values in (x,y)
  if (reverse) {
    std::sort(indices.begin(), indices.end(),
	      [this](size_t i1, size_t i2) {
		return std::hypot(m_x_ptr->getData().at(i1), m_y_ptr->getData().at(i1)) > 
		       std::hypot(m_x_ptr->getData().at(i2), m_y_ptr->getData().at(i2));
	      });
  } else {
    std::sort(indices.begin(), indices.end(),
	      [this](size_t i1, size_t i2) {
		return std::hypot(m_x_ptr->getData().at(i1), m_y_ptr->getData().at(i1)) <
		       std::hypot(m_x_ptr->getData().at(i2), m_y_ptr->getData().at(i2));
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
  FloatColPtr c_ptr = it->second; //dynamic_pointer_cast<FloatCol>(it->second);
  
  std::vector<size_t> indices(c_ptr->size()); // index vector
  std::iota(indices.begin(), indices.end(), 0); // fill with 0, 1, ..., n-1

  if (m_verbose)
    std::cerr << "...sorting " << AddCommas(CellCount()) << " cells" << std::endl;
  
  // sort indices based on comparing values in (x,y)
  if (reverse) {
    std::sort(indices.begin(), indices.end(),
	      [&c_ptr](size_t i1, size_t i2) {
		return c_ptr->getData().at(i1) > c_ptr->getData().at(i2);
	      });
  } else {
    std::sort(indices.begin(), indices.end(),
	      [&c_ptr](size_t i1, size_t i2) {
		return c_ptr->getData().at(i1) < c_ptr->getData().at(i2);
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

  validate();

  float xmax = 0;
  float ymax = 0;  

  size_t ncells = CellCount();
  std::vector<double> coords(ncells * 2);
  
  for (size_t i = 0; i < ncells; i++) {
    coords[i*2  ] = m_x_ptr->getData().at(i);
    coords[i*2+1] = m_y_ptr->getData().at(i);
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
  FloatColPtr d_label = std::make_shared<FloatCol>();
  FloatColPtr d_size  = std::make_shared<FloatCol>();
  
  // count the number of nodes for each component
  std::unordered_map<int, size_t> dcount;
  for (const auto& c : pointToComponentId) {
    dcount[c.second]++;
  }
  
  // fill the data into the columns
  for (size_t i = 0; i < m_x_ptr->size(); i++) {
    JPoint p = {m_x_ptr->getData().at(i),m_y_ptr->getData().at(i)};    
    
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
      points[i] = Point(m_x_ptr->getData().at(i), m_y_ptr->getData().at(i));
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

  validate();

  std::vector<FloatColPtr> col_ptr;
  for (const auto& t : m_header.GetDataTags()) {
    assert(m_table.find(t.id) != m_table.end());
    col_ptr.push_back(m_table.at(t.id));
  }

  for (size_t i = 0; i < numRows; i++) {

    Cell cell;
    
    cell.id         = m_id_ptr->getData().at(i);
    cell.cflag      = m_cflag_ptr->getData().at(i); 
    cell.pflag      = m_pflag_ptr->getData().at(i);
    cell.x          = m_x_ptr->getData().at(i); 
    cell.y          = m_y_ptr->getData().at(i); 

    // if there are any filters, pass if not in the set
    if (m_cells_to_write.size())  {
      if (!m_cells_to_write.count(cell.id))
	continue;
    }

    for (const auto& c : col_ptr) {
      cell.cols.push_back(c->getData().at(i));
    }

    // write it
    (*m_archive)(cell);
  }
}

void CellTable::Crop(float xlo, float xhi, float ylo, float yhi) {
  
  validate();
  
  // total number of rows of table
  size_t nc = CellCount();
  
  std::vector<size_t> valid_indices;
  for (size_t i = 0; i < nc; ++i) {
    
    float x_value = m_x_ptr->getData().at(i);
    float y_value = m_y_ptr->getData().at(i);    
    
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
  std::vector<std::string> labels;
  std::vector<std::vector<float>> data;
  for (const auto& t : m_header.GetDataTags()) {
    if (t.type == Tag::MA_TAG) {
      auto ptr = m_table.find(t.id);
      assert(ptr != m_table.end());
      data.push_back(ptr->second->getData());
      labels.push_back(ptr->first);
    }
  }    
  
  // get the correlation matrix
  size_t n = data.size();
  std::vector<std::vector<float>> correlation_matrix(n, std::vector<float>(n, 0));
  for (size_t i = 0; i < n; i++) {
    for (size_t j = i + 1; j < n; j++) {
      correlation_matrix[i][j] = pearsonCorrelation(data[i], data[j]); 
    }
  }
  
  
  // if csv, print to csv instead of triangle table
  if (csv) {
    print_correlation_matrix(labels, correlation_matrix, sort);
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
    std::cout << std::setw(spacing) << labels.at(i);
  }
  std::cout << std::endl;

  // print the matrix
  for (size_t i = 0; i < n; i++) {
    
    // Print y-axis marker
    std::cout << std::setw(spacing) << labels.at(i);
    
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
  validate();
  
  float x_min = m_x_ptr->Min();
  float x_max = m_x_ptr->Max();
  float y_min = m_y_ptr->Min();
  float y_max = m_y_ptr->Max();
  
  // scale it to fit on the plot
  size_t nc = CellCount();
  std::vector<std::pair<int, int>> scaled_coords;
  for (size_t i = 0; i < nc; i++) {
    float x_ = m_x_ptr->getData().at(i);
    float y_ = m_y_ptr->getData().at(i);    
    
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

void CellTable::DeleteColumn(const std::string& col) {

  if (m_table.find(col) == m_table.end()) {
    throw std::runtime_error("Trying to delete column that doesn't exists: " + col);
  }
  
  m_table.erase(col);
  
}

void CellTable::AddICPXYColumn(ColPtr value,
			       const std::string& type) {

  if (value->size() != CellCount()) {
    throw std::runtime_error("Adding column of incorrect size");
  }
  
  // add to the table
  if (type == "x") {
    m_x_ptr = std::dynamic_pointer_cast<FloatCol>(value);
  } else if (type == "y") {
    m_y_ptr = std::dynamic_pointer_cast<FloatCol>(value);
  } else if (type == "id") {
    m_id_ptr = std::dynamic_pointer_cast<IDCol>(value);
  } else if (type == "pflag") {
    m_pflag_ptr = std::dynamic_pointer_cast<IntCol>(value);
  } else if (type == "cflag") {
    m_cflag_ptr = std::dynamic_pointer_cast<IntCol>(value);
  } else {
    assert(false);
  }
  
  return;
}

void CellTable::AddColumn(const Tag& tag,
			  FloatColPtr value) {
  
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
  validate();
  
  // Initialize the "sparoi" column
  FloatColPtr new_data = std::make_shared<FloatCol>();
  new_data->reserve(nc);

  // Loop the table and check if the cell is in the ROI
  for (size_t i = 0; i < nc; i++) {
    float x_ = m_x_ptr->getData().at(i);
    float y_ = m_y_ptr->getData().at(i);    

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

void CellTable::print_correlation_matrix(const std::vector<std::string>& labels,
					 const std::vector<std::vector<float>>& correlation_matrix,
					 bool sort) const {

  assert(labels.size() == correlation_matrix.size());
  std::vector<std::tuple<int, int, double>> sorted_correlations;
  
  for (int i = 0; i < labels.size(); i++) {
    for (int j = i + 1; j < labels.size(); j++) {
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
    
    std::cout << labels.at(i) << "," << labels.at(j) << "," << corr << std::endl;
  }
}

void CellTable::initialize_cols() {
  
  // initialize cell id
  m_id_ptr = std::make_shared<IDCol>();
  m_cflag_ptr = std::make_shared<IntCol>();
  m_pflag_ptr = std::make_shared<IntCol>();
  m_x_ptr = std::make_shared<FloatCol>();
  m_y_ptr = std::make_shared<FloatCol>();  
  
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

void CellTable::PrintJaccardSimilarity(bool csv, bool sort) const {

  // collect the marker data in one structure
  std::vector<Tag> markers = m_header.GetMarkerTags();

  validate();
  
  std::vector<std::vector<bool>> data;
  std::vector<std::string> labels;
  
  // loop each marker (each bit in the pflags)
  for (size_t bitPos = 0; bitPos < markers.size(); bitPos++) {

    // Extract bits from the current position from all elements
    std::vector<bool> bits;
    for (size_t j = 0; j < m_pflag_ptr->size(); j++) {
      cy_uint val = m_pflag_ptr->getData().at(j);
      bits.push_back((val >> bitPos) & 1);
    }

    labels.push_back(markers.at(bitPos).id);
    data.push_back(bits);
  }    
  
  // get the correlation matrix
  size_t n = data.size();
  std::vector<std::vector<float>> correlation_matrix(n, std::vector<float>(n, 0));
  for (size_t i = 0; i < n; i++) {
    for (size_t j = i + 1; j < n; j++) {
      correlation_matrix[i][j] = jaccardSimilarity(data.at(i), data.at(j));
    }
  }
  
  // if csv, print to csv instead of triangle table
  if (csv) {
    print_correlation_matrix(labels, correlation_matrix, sort);
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
  for (const auto& t : markers)
    if (t.id.length() > spacing)
      spacing = t.id.length() + 2;
  
  // Print x-axis markers
  std::cout << std::setw(spacing) << " ";
  for (size_t i = 0; i < n; i++) {
    std::cout << std::setw(spacing) << labels.at(i);
  }
  std::cout << std::endl;

  // print the matrix
  for (size_t i = 0; i < n; i++) {
    
    // Print y-axis marker
    std::cout << std::setw(spacing) << labels.at(i);
    
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

int CellTable::PlotPNG(const std::string& file,
		       float scale_factor,
		       const std::string& colname
		       ) const {
  
#ifdef HAVE_CAIRO

  std::random_device rd; 
  std::mt19937 gen(rd()); 
  std::uniform_int_distribution<> dis(0,1);
  
  const float radius_size = 3.0f;
  const float alpha_val = 0.5f;
  constexpr float TWO_PI = 2.0 * M_PI;
  
  // get the x y coordinates of the cells
  validate();
  
  // Open the PNG for drawing
  const int width  = m_x_ptr->Max();
  const int height = m_y_ptr->Max();    
  
  // open PNG for drawing
  cairo_surface_t *surfacep = cairo_image_surface_create (CAIRO_FORMAT_ARGB32, width*scale_factor, height*scale_factor);
  cairo_t *crp = cairo_create (surfacep);
  cairo_set_source_rgb(crp, 0, 0, 0); // background color
  
  // loop and draw
  size_t count = 0;
  for (size_t j = 0; j < CellCount(); j++) { // loop the cells
    
    // Draw the arc segment
    const float x = m_x_ptr->getData().at(j);
    const float y = m_y_ptr->getData().at(j);

    const cy_uint pf = m_pflag_ptr->getData().at(j);
    const cy_uint cf = m_cflag_ptr->getData().at(j);
    
    float start_angle = 0.0;
    
    CellFlag pflag(pf);
    CellFlag cflag(cf);
    
    Color c;
    
    if (pflag.testAndOr(4416,0) && pflag.testAndOr(32768,0)) // T-cell - PD-1+ 
      c = color_red;
    else if (pflag.testAndOr(4416,0) && !pflag.testAndOr(32768,0)) // T-cell - PD-1-
      c = color_light_red;
    else if (pflag.testAndOr(1024,0)) // B-cell
      c = color_purple;
    else if (pflag.testAndOr(0,18432) || pflag.testAndOr(0,133120)) // PD-L1 POS tumor cell
      c = color_dark_green;
    else if (pflag.testAndOr(147456,0)) // PD-L1 NEG tumor cell
      c = color_light_green;
    else if (pflag.testAndOr(2048,0)) // PD-L1 any cell
      c = color_cyan;      
    else
      c = color_gray;
    
    cairo_set_source_rgba(crp, c.redf(), c.greenf(), c.bluef(), alpha_val);
    cairo_arc(crp, x*scale_factor, y*scale_factor, radius_size, 0, TWO_PI);
    cairo_line_to(crp, x*scale_factor, y*scale_factor);
    cairo_fill(crp);
    
    // red radius
    /*    if (pflag.testAndOr(147456,0) && j % 100000 == 0) { //dis(gen) < 0.00002)  {
      cairo_set_source_rgb(crp, 1, 0, 0);
      cairo_set_line_width(crp, 4);
      cairo_arc(crp, x*scale_factor, y*scale_factor, 200.0f/0.325f*scale_factor, 0, 2*M_PI);
      cairo_stroke(crp);
      }*/
  }

  int legend_width = 3500*scale_factor;
  int legend_height = 400*scale_factor; // Height of each color box
  int font_size = 70;
  int legend_padding = 20;
  int legend_x = width*scale_factor - legend_width - legend_padding;
  int legend_y = legend_padding;

  ColorMap cm = {color_red, color_light_red,
		 color_purple, color_dark_green, color_light_green,
		 color_cyan, color_gray};
  std::vector<std::string> labels = {
    "T-cell PD-1 pos",
    "T-cell PD-1 neg",    
    "B-cell",
    "Tumor-cell - PD-L1 pos",
    "Tumor-cell - PD-L1 neg",
    "Other PD-L1 pos",
    "Stromal"
  };
  
  add_legend_cairo(crp, font_size, legend_width, legend_height, legend_x, legend_y,
		   cm, labels);
  
  cairo_destroy (crp);
  cairo_surface_write_to_png (surfacep, file.c_str());
  cairo_surface_destroy (surfacep);
  
#else
  std::cerr << "Unable to make PNG -- need to include / link Cairo library" << std::endl;
#endif
  
  return 0;
}

void CellTable::ScramblePflag(int seed, bool lock_flags) {

  validate();
  
  // flags are locked, so we are permuting flags
  if (lock_flags) {
    m_pflag_ptr->Scramble(seed);
    return;
  }

  // store new flags
  std::vector<cy_uint> new_flags(m_pflag_ptr->size());

  // make a random number generator
  std::default_random_engine rng(seed);
  
  // flags are not locks, so we are permutting individual bits
#ifdef USE_64_BIT  
  for (int bitPos = 0; bitPos < 64; ++bitPos) {
#else
  for (int bitPos = 0; bitPos < 32; ++bitPos) {
#endif
    
    // Extract bits from the current position from all elements
    std::vector<bool> bits;
    for (size_t i = 0; i < new_flags.size(); i++) {
      cy_uint val = m_pflag_ptr->getData().at(i);
      bits.push_back((val >> bitPos) & 1);
    }
    
    // Shuffle these bits
    std::shuffle(bits.begin(), bits.end(), rng);
    
    // Set these shuffled bits back to the vector
    for (size_t i = 0; i < new_flags.size(); ++i) {
      new_flags[i] &= ~(1u << bitPos);           // Clear the bit at bitPos (should already be cleared)
      new_flags[i] |= (bits[i] << bitPos);       // Set the bit at bitPos if bits[i] is 1
    }
  }

  // assign the new flags
  m_pflag_ptr = std::make_shared<NumericColumn<cy_uint>>(new_flags);
  validate();
  
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
 
 void CellTable::validate() const {
   
   assert(m_x_ptr);
   assert(m_y_ptr);
   assert(m_pflag_ptr);
   assert(m_cflag_ptr);
   assert(m_id_ptr);

   assert(m_x_ptr->size() == m_y_ptr->size());
   assert(m_x_ptr->size() == m_pflag_ptr->size());
   assert(m_x_ptr->size() == m_cflag_ptr->size());
   assert(m_x_ptr->size() == m_id_ptr->size());         
   
 }
