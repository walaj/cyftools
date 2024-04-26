#include "cell_table.h"
#include "polygon.h"
#include "cysift.h"
#include "csv.h"
#include "color_map.h"
#include "cell_selector.h"
#include "cell_utils.h"
#include "cell_graph.h"

#include <random>
#include <cstdlib>

#ifdef HAVE_TIFFLIB
#include "tiff_writer.h"
#endif

#ifdef HAVE_CAIRO
#include "cairo/cairo.h"
#include "cairo/cairo-pdf.h"
#endif

//#ifdef HAVE_HDF5
//#include <H5Cpp.h>
//#endif

bool CellTable::HasColumn(const std::string& col) const {
  return m_table.find(col) != m_table.end();
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

#include <iomanip>  // for std::setw and std::left

std::ostream& operator<<(std::ostream& os, const CellTable& table) {

    // make the cell and sample ids
    std::shared_ptr<NumericColumn<uint32_t>> sid = make_shared<NumericColumn<uint32_t>>();
    std::shared_ptr<NumericColumn<uint32_t>> cid = make_shared<NumericColumn<uint32_t>>();
    cid->reserve(table.size());
    sid->reserve(table.size());
    for (const auto& c : table.m_id_ptr->getData()) {
        sid->PushElem(static_cast<uint32_t>(c >> 32));
        cid->PushElem(static_cast<uint32_t>(c & 0xFFFFFFFF));
    }

    // Decide on a field width for the descriptions; 
    // this can be adjusted based on your expected max width or based on the actual data.
    const int descWidth = 12;
    const int dataWidth = 50;  // Adjust this based on the width of your data.

    os << std::left << std::setw(descWidth) << "CellID"    << " -- " << cid->toString()           << std::endl;
    os << std::left << std::setw(descWidth) << "SampleID"  << " -- " << sid->toString()           << std::endl;
    os << std::left << std::setw(descWidth) << "RawID"     << " -- " << table.m_id_ptr->toString() << std::endl;
    os << std::left << std::setw(descWidth) << "Pheno Flag"   << " -- " << table.m_pflag_ptr->toString()  << std::endl;
    os << std::left << std::setw(descWidth) << "Cell Flag"    << " -- " << table.m_cflag_ptr->toString()   << std::endl;
    os << std::left << std::setw(descWidth) << "X"      << " -- " << table.m_x_ptr->toString() << std::endl;
    os << std::left << std::setw(descWidth) << "Y"      << " -- " << table.m_y_ptr->toString() << std::endl;
    
    for (const auto& t : table.m_header.GetDataTags()) {
        
        const FloatColPtr col_ptr = table.m_table.at(t.id);
        std::string ctype;

        switch (t.type) {
            case Tag::MA_TAG: ctype = "Marker"; break;
            case Tag::GA_TAG: ctype = "Graph"; break;
            case Tag::CA_TAG: ctype = "Meta"; break;
            default: ctype = "UNKNOWN"; break;
        }
        
        os << std::left << std::setw(descWidth) << t.id << " -- " << std::setw(dataWidth) << ctype << " -- " << col_ptr->toString() << std::endl;
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
      assert(false);
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
#ifdef HAVE_OMP  
#pragma omp parallel for num_threads(m_threads)
#else
  std::cerr << "OMP not included, no support for multithreading. Compile with to support" << std::endl;
#endif  
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

void CellTable::PrintJaccardSimilarity(bool csv, bool sort, bool subset_score) const {

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
      correlation_matrix[i][j] = subset_score ?
	jaccardSubsetSimilarity(data.at(i), data.at(j)) :
        jaccardSimilarity(data.at(i), data.at(j));      
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
		       const std::string& module,
		       const std::string& roifile,
		       const std::string& title,
		       ColorLabelVec palette
		       ) const {
  
#ifdef HAVE_CAIRO

  validate();
  
  std::random_device rd; 
  std::mt19937 gen(rd()); 
  std::uniform_int_distribution<> dis(0,1);

  float micron_per_pixel = 0.325f;
  const float radius_size = 6.0f * scale_factor;
  const float ALPHA_VAL = 0.7f; // alpha for cell circles
  constexpr float TWO_PI = 2.0 * M_PI;
  
  // Original dimensions
  const int original_width = m_x_ptr->Max();
  const int original_height = m_y_ptr->Max();    

  // Additional height to accommodate the legend at the top
  const int legend_total_height = original_height * 0.05; // Adjust as needed

  // New dimensions with extra space for the legend
  const int width = original_width;
  const int height = original_height + legend_total_height; // Increase height for legend

  // open PNG for drawing
  cairo_surface_t *surfacep = cairo_image_surface_create(CAIRO_FORMAT_ARGB32, width*scale_factor, height*scale_factor);
  cairo_t *crp = cairo_create(surfacep);
  
  // Fill background with black for the legend area
  cairo_set_source_rgb(crp, 0, 0, 0);
  cairo_rectangle(crp, 0, 0, width*scale_factor, legend_total_height*scale_factor);
  cairo_fill(crp);
  
  // Set background color for the rest
  cairo_set_source_rgb(crp, 1, 1, 1); // Assuming white is the desired background color
  cairo_rectangle(crp, 0, legend_total_height*scale_factor, width*scale_factor, height*scale_factor);
  cairo_fill(crp);

  // get a pre-selected module if the map is zero
  if (palette.empty()) {
    palette = ColorLabelVecForModule(module);
  } 

  if (m_verbose)
    for (const auto& pp : palette)
      std::cerr << pp << std::endl;
  
  // defined below in color_map.cpp
  //ColorLabelMap cm = ColorLabelMapForModule(module);

  // loop and draw
  size_t count = 0;
  for (size_t j = 0; j < CellCount(); j++) { // loop the cells
    
    // Draw the arc segment
    const float x = m_x_ptr->at(j);
    const float y = m_y_ptr->at(j) + legend_total_height;
    
    const cy_uint pf = m_pflag_ptr->at(j);
    const cy_uint cf = m_cflag_ptr->at(j);
    
    float start_angle = 0.0;

    // get the P and C-flags
    CellFlag pflag(pf);
    CellFlag cflag(cf);

    // select the color for this cell from the palette
    Color c = select_color(cf, pf, palette);
    
    cairo_set_source_rgba(crp, c.redf(), c.greenf(), c.bluef(), ALPHA_VAL);
    cairo_arc(crp, x*scale_factor, y*scale_factor, radius_size, 0, TWO_PI);
    cairo_fill(crp); 

  }

  // red radius
  if (false)
  for (size_t j = 0; j < CellCount(); j++) { // loop the cells
    const float x = m_x_ptr->at(j);
    const float y = m_y_ptr->at(j) + legend_total_height;
    if (j % 10000 == 0) { 
      cairo_set_source_rgb(crp, 1, 0, 0);
      cairo_set_line_width(crp, 16);
      cairo_arc(crp, x*scale_factor, y*scale_factor, 50.0f/0.325f*scale_factor, 0, 2*M_PI);
      cairo_stroke(crp);
    }
  }

  
  ////////
  // ROI
  // read in the roi file
  if (!roifile.empty()) {
    std::vector<Polygon> rois = read_polygons_from_file(roifile);
    
    for (const auto& polygon : rois) {

      if (polygon.Text.find("tumor") == std::string::npos && polygon.Name.find("tumor") == std::string::npos &&
	  polygon.Text.find("normal") == std::string::npos && polygon.Name.find("normal") == std::string::npos &&
	  polygon.Text.find("blacklist") == std::string::npos && polygon.Name.find("blacklist") == std::string::npos &&
	  polygon.Text.find("rtifact") == std::string::npos && polygon.Name.find("rtifact") == std::string::npos) {
	std::cerr << "unknown roi name: " << polygon.Text << " -- " << polygon.Name << std::endl;
	continue;
      }
      
      if (polygon.size() == 0)
	continue;
      
      // Move to the first vertex
      auto firstVertex = *polygon.begin();
      cairo_move_to(crp, firstVertex.x*scale_factor*micron_per_pixel,
    (firstVertex.y*micron_per_pixel+legend_total_height)*scale_factor);
      
      // Draw lines to each subsequent vertex
      for (const auto& v : polygon) {
	
	// draw the points
	//cairo_set_source_rgba(crp, 1, 0, 0, 1);
	//cairo_arc(crp, v.first*scale_factor*micron_per_pixel, v.second*scale_factor*micron_per_pixel, 10, 0, TWO_PI);
	//cairo_fill(crp);
	
	cairo_line_to(crp, v.x*scale_factor*micron_per_pixel, (v.y*micron_per_pixel+legend_total_height)*scale_factor);
      }
      
      // Close the polygon
      cairo_close_path(crp);

      const float alphaval = 0.2;
      // Set the source color for fill and line (red with high alpha for transparency)
      if (polygon.Text.find("tumor") != std::string::npos || polygon.Name.find("tumor") != std::string::npos)
	cairo_set_source_rgba(crp, 1, 0, 0, alphaval); // Red with alpha = alphaval
      else if (polygon.Text.find("normal") != std::string::npos || polygon.Name.find("normal") != std::string::npos) 
	cairo_set_source_rgba(crp, 1, 1, 0, alphaval); // Yelslow with alpha = alphaval
      else if (polygon.Text.find("eminal") != std::string::npos || polygon.Name.find("normal") != std::string::npos) 
	cairo_set_source_rgba(crp, 0, 1, 0, alphaval); // Grexen with alpha = alphaval
      else if (polygon.Text.find("blacklist") != std::string::npos || polygon.Name.find("blacklist") != std::string::npos ||
	       polygon.Text.find("rtifact") != std::string::npos || polygon.Name.find("rtifact") != std::string::npos)
	cairo_set_source_rgba(crp, 0.5, 0.5, 0.5, alphaval); // Gray with alpha = 0.5	
      else
	assert(false);
      
      // Fill the polygon
      cairo_fill_preserve(crp);
      
      // Set the source color for the boundary (solid red)
      cairo_set_source_rgb(crp, 1, 0, 0); // Solid red
      
      // Draw the boundary
      cairo_stroke(crp);
    }
  }

  ///////
  // CONVEX HULL
  // compute and display the convex hull
  std::unordered_map<std::string, FloatColPtr>::const_iterator it = m_table.find("tls_id");
  std::unordered_map<float, std::vector<JPoint>> hull_map;
  int n = CellCount();
  if (it != m_table.end()) {
    FloatColPtr cptr = it->second; // This is your reference/pointer to the FloatColPtr
    // Use a set to find unique elements, since sets automatically
    // remove duplicates and store elements in sorted order
    std::set<float> unique_clusters(cptr->begin(), cptr->end());

    // loop the clusters and make the convex hull
    for (const auto& cl : unique_clusters) {
      
      // cluster 0 is holder for not a cluster
      if (cl == 0 || true) // debug - true means dont' plot
	continue;

      // fill polygon with the points that will need to have hull around them
      std::vector<JPoint> polygon;
      polygon.reserve(n);
      for (size_t i = 0; i < n; i++) {
	if (cptr->at(i) == cl)
	  polygon.push_back(JPoint(m_x_ptr->at(i), m_y_ptr->at(i)));
      }
      
      // get the convex hull
      hull_map[cl] = convexHull(polygon);

      if (m_verbose)
	std::cerr << "...cyftools png - hull map for " << cl << " has " << hull_map[cl].size() << " points " << std::endl;
    }

    // draw the hulls
    ///////
    cairo_set_source_rgb(crp, 1.0, 0.0, 0.0);
    cairo_set_line_width(crp, 60.0*scale_factor); // Set the line width to 5.0
    // loop through invidual hulls
    for (auto& hullm : hull_map) {
      auto& hull = hullm.second;
      
      // loop though individual hull points
      for (size_t i = 0; i < hull.size(); ++i) {
	const auto& start = hull.at(i);
	const auto& end = hull.at((i + 1) % hull.size()); // Wrap around to first point
	
	cairo_move_to(crp, start.x*scale_factor, (start.y+legend_total_height)*scale_factor);
	cairo_line_to(crp, end.x*scale_factor, (end.y+legend_total_height)*scale_factor);
      }
      cairo_stroke(crp); // Actually draw the lines
    }
    
  } else {
    if (m_verbose)
      std::cerr << "...cyftools png - no tls clusters to convex hull" << std::endl;
  }

  
  ///////
  // LEGEND
  int font_size = 200*scale_factor;
  
  add_legend_cairo_top(crp, font_size,
		       legend_total_height*scale_factor,
		       width*scale_factor,
		       palette); 

  // this is where the y is for bottom of text
  float y_base = legend_total_height*scale_factor/2 + font_size/2;
  
  // add the title to the legend
  if (!title.empty()) {
    cairo_select_font_face(crp, "Arial",
			   CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_BOLD);
    cairo_set_font_size(crp, font_size); // Adjust font size to your needs

    // Draw the title
    cairo_set_source_rgb(crp, 255, 255, 255); // Set color to white for the text
    cairo_move_to(crp, 100*scale_factor, y_base);
    cairo_show_text(crp, title.c_str());
  }

  ///////
  // scale bar
  ///////
  //std::cerr << "...cyftools png - drawing scale-bar assuming " << micron_per_pixel << " microns / pixel" << std::endl;
  draw_scale_bar(crp, width*scale_factor*0.85, y_base, 1000*scale_factor, 50*scale_factor, "1 mm");
  
  // clean up the PNG
  cairo_destroy (crp);
  cairo_surface_write_to_png (surfacep, file.c_str());
  cairo_surface_destroy (surfacep);

#else
  std::cerr << "Unable to make PNG -- need to include / link Cairo library" << std::endl;
#endif
  
  return 0;
}

size_t CellTable::CountCFlag(int flag) const {

  assert(flag > 0);
  
  size_t count = 0;
  const size_t n = CellCount();
  for (int i = 0; i < n; i++)
    if (IS_FLAG_SET(m_cflag_ptr->at(i), flag))
      count++;
  return count;
}


void CellTable::ClearCFlag(int flag) {
  const size_t n = CellCount();
  validate();
  for (int i = 0; i < n; i++)
    CLEAR_FLAG((*m_cflag_ptr)[i], flag);
}

void CellTable::CopyCFlag(cy_uint flag_from, cy_uint flag_to) {

  assert(flag_from > 0);
  assert(flag_to > 0);
  
  const size_t n = CellCount();

  for (size_t i = 0; i < n; i++) {
    cy_uint& cf  = (*m_cflag_ptr)[i];
    if (IS_FLAG_SET(cf, flag_from))
      SET_FLAG(cf, flag_to);
  }
  
}

void CellTable::FlagToFlag(const bool clear_flag_to,
			   const int flag_from, bool flag_from_negative, 
			   const int flag_to, bool flag_to_negative) {

  validate();
  const size_t n = CellCount();
  
  // set flags all on or all off
  if (clear_flag_to) {
    for (size_t i = 0; i < n; i++) {
      if (!flag_to_negative)    
	CLEAR_FLAG((*m_cflag_ptr)[i], flag_to);
      else
	SET_FLAG((*m_cflag_ptr)[i], flag_to);
    }
  }
  
  // transfer
  for (size_t i = 0; i < n; i++) {
    
    if (!flag_from_negative &&
	!flag_to_negative && 
	IS_FLAG_SET(m_cflag_ptr->at(i), flag_from))
      SET_FLAG((*m_cflag_ptr)[i], flag_to);
    else if (flag_from_negative &&
	     !flag_to_negative &&
	     !IS_FLAG_SET(m_cflag_ptr->at(i), flag_from))
      SET_FLAG((*m_cflag_ptr)[i], flag_to);
    else if (!flag_from_negative &&
	     flag_to_negative &&
	     IS_FLAG_SET(m_cflag_ptr->at(i), flag_from))
      CLEAR_FLAG((*m_cflag_ptr)[i], flag_to);
    else if (flag_from_negative &&
	     flag_to_negative &&
	     !IS_FLAG_SET(m_cflag_ptr->at(i), flag_from))
      CLEAR_FLAG((*m_cflag_ptr)[i], flag_to);
  }
}

void CellTable::ScramblePflag(int seed, bool lock_flags, bool phenotype_only) {

  
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

#ifdef USE_64_BIT
#define BIT_LITERAL 1ULL // Unsigned long long literal for 64-bit
#else
#define BIT_LITERAL 1u // Unsigned integer literal for 32-bit
#endif
  
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

    // skip if nothing to scramble
    int sum = std::accumulate(bits.begin(), bits.end(), 0);
    if (!sum) {
      continue;
    }
    
    // Shuffle these bits across cells
    std::shuffle(bits.begin(), bits.end(), rng);
    
    // Set these shuffled bits back to the vector
    for (size_t i = 0; i < new_flags.size(); ++i) {
      new_flags[i] &= ~(BIT_LITERAL << bitPos);
      new_flags[i] |= (static_cast<cy_uint>(bits[i]) << bitPos);      
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
  Cell cell;
  while (true) {
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
  
void CellTable::StreamTableCSV(CerealProcessor& proc, const std::string& file) {

  // Assume x and y indiceis are the first and second
  int x_index = 0;
  int y_index = 1;

  // assume starts at 2 and ends and end of file
  int start_index = 0; // start and end indicies for markers
  int end_index = 0;
  
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

      // assumng "jeremiah" format of X,Y,markers...
      proc.SetXInd(0);
      proc.SetYInd(1);
      proc.SetStartIndex(2);
      proc.SetEndIndex(1000000); // just go to end of line they're all markers
      
      // add the tag to the header
      Tag tag(line);
      m_header.addTag(tag);

    }

    // line starts with CellID (assuming...)
    else if (line[0] == 'C') {

      // make sure there are no header lines in the middle of file
      if (header_read)
	throw std::runtime_error("Misformed file: header lines should all be at top of file");

      std::vector<std::string> header_lines = tokenize_comma_delimited(line);

      size_t ind = 0;
      for (const auto& s : header_lines) {
	
	if (is_mcmicro_meta(s)) {
	  if (s == "X_centroid") 
	    x_index = ind;
	  else if (s == "Y_centroid")
	    y_index = ind;
	  // right now, do nothing with meta data
	  
	  // if we're already past markers
	  if (start_index > 0 && end_index == 0)
	    end_index = ind; // this is the last index seen
	}

	// this is a marker, add the tag
	else {

	  // clean out the ARgo etc
	  std::string s2 = clean_marker_string(s);
	  
	  Tag tag(Tag::MA_TAG, s2, "");
	  m_header.addTag(tag);
	  // if first marker we've seen, then set start
	  if (start_index == 0) {
	    start_index = ind;
	  }
	}

	// just a sanity check here
	if (ind == 0 && s != "CellID") {
	  std::cerr << "Error: cyftools convert -- saw header in csv as starting with C but its not CellID, double check?" << std::endl;
	  assert(false);
	} else if (ind == 1 && s != "Hoechst") {
	  std::cerr << "Warning: cyftools conveort -- usually assume mcmicro header is CellID,Hoechst,... but see here that 2nd elem is " << s << std::endl;
	  std::cerr << "         OK to proceed, but will assume " << s << " is first marker." << std::endl;
	} 
	
	ind++;
      }

      // we know X_centroid isn't at 0 slow, because we got into here with 'C' as first letter
      // so what happened is that it never found X_centroid
      if (x_index == 0) {
	std::cerr << "Error: cyftools convert -- X_centroid not found" << std::endl;
	assert(false);
      }
      if (y_index == 0) {
	std::cerr << "Error: cyftools convert -- Y_centroid not found" << std::endl;
	assert(false);
      }
      
      // some more sanity check
      assert(start_index > 0);
      assert(end_index > 0);
      
      proc.SetXInd(x_index);
      proc.SetYInd(y_index);
      proc.SetStartIndex(start_index); // assuming that 0 is CellID
      proc.SetEndIndex(end_index);

    }


    else {
      
      // Label that the header has already been read
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

  // validate all data cells are same length
  size_t prev_size = m_x_ptr->size();
  size_t n = 0;
  // loop the table and find the maximum length
  for (const auto &c : m_table) {
    size_t current_size = c.second->size();
    
    if (current_size != prev_size) {
      std::cerr << "Warning: Column sizes do not match. Column: " <<
	c.first << " prev size " << prev_size <<
	" current_size " << current_size << std::endl;
      assert(false);
    }
    
    prev_size = current_size;
    n = std::max(n, current_size);
  }
}
 
#ifdef HAVE_LDAPLUSPLUS
 const std::vector<std::vector<double>> CellTable::create_inverse_distance_weights() {
   
  validate();
  int n = m_x_ptr->size();
  
  std::vector<std::vector<double>> W(n, std::vector<double>(n, 0)); // Initialize the weights matrix
  
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      if (i != j) { // Avoid self-comparison
	// Calculate distance
	double dist = std::sqrt(
				std::pow(m_x_ptr->at(i) - m_x_ptr->at(j), 2) + 
				std::pow(m_y_ptr->at(i) - m_y_ptr->at(j), 2)
				);
	if (dist > 0) { // Avoid division by zero
	  W[i][j] = 1 / dist;
	}
      }
    }
  }
  
  // Set the diagonal to zero to avoid self-influence
  for (int i = 0; i < n; ++i) {
    W[i][i] = 0;
  }
  
  return W;
}

#endif

