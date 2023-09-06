#include "cell_processor.h"
#include "polygon.h"
#include "cell_utils.h"
#include "cell_flag.h"
#include "cell_selector.h"

#include "cell_row.h"

#include <limits>

int SelectProcessor::ProcessHeader(CellHeader& header) {
  
  m_header = header;

  // find which index the field to select is
  size_t i = 0;
  for (const auto& t : header.GetDataTags()) {

    auto it = m_criteria.find(t.id);

    // convert the string -> SelectOp map to a int -> SelectOp map
    if (it != m_criteria.end()) {
      m_criteria_int[i] = it->second;
    }
    i++;
  }

  if (m_criteria.size() != m_criteria_int.size())
    throw std::runtime_error("Unable to find all requested fields in the header");

  m_header.addTag(Tag(Tag::PG_TAG, "", m_cmd));

  m_header.SortTags();
  
  // just in time, make the output stream
  this->SetupOutputStream();
  
  // output the header
  assert(m_archive);
  (*m_archive)(m_header);
  
  return HEADER_NO_ACTION;
}

int AverageProcessor::ProcessHeader(CellHeader& header) {

  m_header = header;

  // which columns store data
  /*for (size_t i = 0; i < m_header.GetAllTags()) {
    const Tag& tag = m_header.GetAllTags().at(i);
    if (tag.type == Tag::CA_TAG || tag.type == Tag::MA_TAG))
    inds.push_back(i);
    }*/

  // initialize sums to zero
  sums.resize(m_header.GetDataTags().size());
  if (sums.size())
    assert(sums.at(0) == 0);
  
  m_header.addTag(Tag(Tag::PG_TAG, "", m_cmd));

  m_header.SortTags();
  
  // just in time, make the output stream
  this->SetupOutputStream();
  
  // output the header
  assert(m_archive);
  (*m_archive)(m_header);
  
  return HEADER_NO_ACTION;
}

int AverageProcessor::ProcessLine(Cell& cell) {

  for (size_t i = 0; i < cell.cols.size(); i++) {
    sums[i] += cell.cols.at(i);
  }
  n++;

  // don't emit anything in StreamTable
  return CellProcessor::NO_WRITE_CELL;
}

void AverageProcessor::EmitCell() const {

  Cell cell;

  std::vector<float> means(sums.size());
  if (n > 0) {
    for (size_t i = 0; i < means.size(); i++) {
      means[i] = sums.at(i) / n;
    }
  }

  cell.cols = means;

  // write the one cell
  OutputLine(cell);
}

int CellCountProcessor::ProcessHeader(CellHeader& header) {
  
  m_header = header;

  // initialize the count vector
  size_t num_marker_tags = header.GetMarkerTags().size();
  m_counts = std::vector<size_t>(num_marker_tags);

  m_header.ClearMeta();

  for (const auto& m : header.GetMarkerTags()) {
    m_header.addTag(Tag(Tag::CA_TAG, m.id + "_count", ""));
  }

  m_header.addTag(Tag(Tag::PG_TAG, "", m_cmd));

  m_header.SortTags();

  
  // just in time, make the output stream
  this->SetupOutputStream();
  
  // output the header
  assert(m_archive);
  (*m_archive)(m_header);
  
  return HEADER_NO_ACTION;
}

int CellCountProcessor::ProcessLine(Cell& cell) {

  for (size_t i = 0; i < m_counts.size(); i++) {
    CellFlag f(cell.pflag);
    cy_uint to_test = static_cast<cy_uint>(std::pow(2, i));
    //    std::cerr << " i " << i << " totest " << to_test << " f " << f << " test " << f.testAndOr(to_test,0) << std::endl;
    if (f.testAndOr(to_test, 0))
      m_counts[i]++;
  }

  return CellProcessor::NO_WRITE_CELL;
}

void CellCountProcessor::EmitCell() const {

  Cell cell;

  // add markers as dummies
  for (size_t i = 0; i < m_counts.size(); i++)
    cell.cols.push_back(0);

  // add the count data
  for (size_t i = 0; i < m_counts.size(); i++)
    cell.cols.push_back(m_counts.at(i));

  // zero the hard data since they are meaningless
  cell.x = 0;
  cell.y = 0;  
  cell.pflag = 0;
  cell.cflag = 0;
  cell.id = 0;

  // write the one cell
  OutputLine(cell);
  
}

int SummaryProcessor::ProcessHeader(CellHeader& header) {
  m_header = header;
  size_t num_data_tags = header.GetDataTags().size();

  // 5 is the number of pre-specified data colums. e.g
  // x, y, cflag, pflag, id
  m_min = std::vector<float>(num_data_tags + 5, std::numeric_limits<float>::max());
  m_max = std::vector<float>(num_data_tags + 5, std::numeric_limits<float>::min());
  m_mean = std::vector<float>(num_data_tags + 5, 0.0f); 

  return HEADER_NO_ACTION; // do nothing  
}

int SummaryProcessor::ProcessLine(Cell& cell) {

  m_count++;

  // fixed mins
  if (cell.x < m_min[0])
    m_min[0] = cell.x;
  if (cell.y < m_min[1])
    m_min[1] = cell.y;
  if (cell.pflag < m_min[2])
    m_min[2] = cell.pflag;
  if (cell.cflag < m_min[3])
    m_min[3] = cell.cflag;
  if (cell.id < m_min[4])
    m_min[4] = cell.id;

  // fixed maxs
  if (cell.x > m_max[0])
    m_max[0] = cell.x;
  if (cell.y > m_max[1])
    m_max[1] = cell.y;
  if (cell.pflag > m_max[2])
    m_max[2] = cell.pflag;
  if (cell.cflag > m_max[3])
    m_max[3] = cell.cflag;
  if (cell.id > m_max[4])
    m_max[4] = cell.id;

  // fixed means
  // assert is to make sure we don't have numeric overflow
  assert(m_mean[0] >= 0); // x and y should be > 0
  assert(m_mean[1] >= 0);
  assert(m_mean[4] >= 0);  
  assert(std::numeric_limits<float>::max() - m_mean[0] > cell.x);
  assert(std::numeric_limits<float>::max() - m_mean[1] > cell.y);
  assert(std::numeric_limits<float>::max() - m_mean[4] > cell.id);
  
  m_mean[0] += cell.x;
  m_mean[1] += cell.y;
  m_mean[4] += cell.id;

  // update remaining columns
  for (size_t i = 0; i < cell.cols.size(); i++) {
    if (m_mean.at(i+5) > 0)
      assert(std::numeric_limits<float>::max() - m_mean.at(i+5) > cell.cols.at(i));
    m_mean[i+5] += cell.cols.at(i);
  }

  return NO_WRITE_CELL;
}

void SummaryProcessor::Print() {

  if (m_count > 0)
    for (auto& m : m_mean)
      m /= static_cast<float>(m_count);
  
  // print x, y, id
  std::cout << "Cell count: " << AddCommas(m_count) << std::endl;
  std::cout << "X range (mean): [" << m_min[0] << "," << m_max[0] << "] (" << m_mean[0] << ")" << std::endl;
  std::cout << "Y range (mean): [" << m_min[1] << "," << m_max[1] << "] (" << m_mean[1] << ")"  << std::endl;
  std::cout << "ID range (mean): [" << m_min[4] << "," << m_max[4] << "] (" << m_mean[4] << ")"  << std::endl;
  std::cout << "pflag range: [" << m_min[2] << "," << m_max[2] << "]" << std::endl;
  std::cout << "cflag range: [" << m_min[3] << "," << m_max[3] << "]" << std::endl;      
  
  // print marker info
  std::vector<Tag> markers = m_header.GetMarkerTags();
  std::cout << "Number of markers: " << AddCommas(markers.size()) << std::endl;
  std::cout << "Markers: ";
  for (size_t i = 0; i < markers.size(); i++) {
    std::cout << markers.at(i).id;
    std::cout << ((i != markers.size() - 1) ? "," : "\n");
  }

  // print meat info
  std::vector<Tag> metas = m_header.GetMetaTags();
  std::cout << "Number of metas: " << AddCommas(metas.size()) << std::endl;
  std::cout << "Metas: ";
  for (size_t i = 0; i < metas.size(); i++) {
    std::cout << metas.at(i).id;
    std::cout << ((i != metas.size() - 1) ? "," : "\n");    
  }
}

int SubsampleProcessor::ProcessHeader(CellHeader& header) {

  m_header = header;

  if (m_rate <= 0 || m_rate > 1)
    throw std::runtime_error("Subsample rate should be between (0,1]). Rate currently: " + std::to_string(m_rate));

  // set the random generator
  std::mt19937 m_gen(m_seed); // Standard mersenne_twister_engine seeded with rd()
  std::uniform_real_distribution<> m_dis(0.0, 1.0); 
  
  m_header.addTag(Tag(Tag::PG_TAG, "", m_cmd));

  m_header.SortTags();
  
  // just in time, make the output stream
  this->SetupOutputStream();
  
  // output the header
  assert(m_archive);
  (*m_archive)(m_header);

  return HEADER_NO_ACTION;
  
}

int SubsampleProcessor::ProcessLine(Cell& line) {

  m_count++;
  
  if (m_count % 500000 == 0 && m_verbose) {
    std::cerr << "...sampled cell " << AddCommas(m_count) << std::endl;
  }

  // kludgy, but this processor will EITHER sample select OR
  // do a rate subsampling, but not both...
  if (m_samplenum.size() > 0) {
    if (m_samplenum.count(line.get_sample_id())) {
      m_kept++;

      if (m_kept % 100000 == 0 && m_verbose) {
	std::cerr << "...kept cell " << AddCommas(m_kept) << std::endl;
      }
      
    }
  }
  
  else if (m_dis(m_gen) <= m_rate) {
    m_kept++;

    if (m_kept % 100000 == 0 && m_verbose) {
      std::cerr << "...kept cell " << AddCommas(m_kept) << " giving rate of " <<
	static_cast<float>(m_kept)/m_count << std::endl;
    }
    
    return WRITE_CELL;
    
  }
  return NO_WRITE_CELL;
  
}

int SelectProcessor::ProcessLine(Cell& cell) {

  bool write_cell = false;

  ///////
  // FLAGS
  ///////
  // NB: even if flags are all empty, default should be to trigger a "write_cell = true"

  //std::cerr << cell.cflag << " - " << cell.pflag << std::endl << m_select << std::endl;
  
  // if flags met, print the cell
  if ( m_select.TestFlags(cell.pflag, cell.cflag)) { 
    write_cell = true;
  }

  //std::cerr << " WRITE CELL " << write_cell << std::endl;
  
  ///////
  // FIELD
  ///////
  bool flag_write_cell = write_cell;
  write_cell = m_or_toggle ? false : write_cell; // if or toggle is on, then ignore flag criteria and set true
  
  for (const auto& c : m_criteria_int) {
    
    assert(c.first < cell.cols.size());
    float value = cell.cols.at(c.first);
    
    for (const auto& cc : c.second) { // loop the vector of {optype, float}

      /*      std::cerr << "m_or_toggle " << m_or_toggle <<
	" criteria int " << c.first << " optype equal? " <<
	(cc.first == optype::EQUAL_TO) << " numeric " << cc.second <<
	" value " << value << " write_cell before " << write_cell << std::endl;
      */
	
      switch (cc.first) { // switching on the optype
      case optype::GREATER_THAN: write_cell          = write_cell = m_or_toggle ? (write_cell || (value >  cc.second)) : (write_cell && (value >  cc.second));  break;
      case optype::LESS_THAN: write_cell             = write_cell = m_or_toggle ? (write_cell || (value <  cc.second)) : (write_cell && (value <  cc.second));  break;
      case optype::GREATER_THAN_OR_EQUAL: write_cell = write_cell = m_or_toggle ? (write_cell || (value >= cc.second)) : (write_cell && (value >= cc.second));  break;
      case optype::LESS_THAN_OR_EQUAL: write_cell    = write_cell = m_or_toggle ? (write_cell || (value <= cc.second)) : (write_cell && (value <= cc.second));  break;
      case optype::EQUAL_TO: write_cell              = write_cell = m_or_toggle ? (write_cell || (value == cc.second)) : (write_cell && (value == cc.second));  break;
      default: assert(false);
      }

      //std::cerr << " write_cell after " << write_cell << std::endl;
    }
    
  }
  
  if (flag_write_cell && write_cell)
    return CellProcessor::WRITE_CELL;
  
  return CellProcessor::NO_WRITE_CELL; // don't write if not selected
}


int HallucinateProcessor::ProcessHeader(CellHeader& header) {
  m_header = header;

  // setup generator
  m_gen = std::mt19937(m_seed); // Seed the generator
  m_distrib = std::uniform_int_distribution<>(0, m_nphenotypes - 1); // Define the range
  
  // remove old markers
  size_t i = 0;
  for (const auto& t : header.GetAllTags()) {
    if (t.type == Tag::MA_TAG)
      m_to_remove.insert(i);
    i++;
  }
  m_header.Cut(m_to_remove);

  // add the new fake tags
  for (size_t i = 0; i < m_nphenotypes; i++) 
    m_header.addTag(Tag(Tag::MA_TAG, "FakeMarker" + std::to_string(i + 1), ""));
  
  m_header.addTag(Tag(Tag::PG_TAG, "", m_cmd));
  m_header.SortTags();

  // just in time output
  this->SetupOutputStream();
  
  // output the header
  assert(m_archive);
  (*m_archive)(m_header);

  return HEADER_NO_ACTION; 
}

int HallucinateProcessor::ProcessLine(Cell& cell) {

  // just transfer to a new set of columns,
  // not including what is to be cut
  std::vector<float> cols_new;
  assert(cell.cols.size() >= m_to_remove.size());
  cols_new.reserve(cell.cols.size() - m_to_remove.size() + m_nphenotypes);

  // add back the old non-marker data
  for (size_t i = 0; i < cell.cols.size(); i++) {
    if (!m_to_remove.count(i)) {
      cols_new.push_back(cell.cols.at(i));
    }
  }

  // the hallucinated data is just zeros in the markers
  for (size_t i = 0; i < cell.cols.size(); i++) {
    cols_new.push_back(0);
  }
  cell.cols = cols_new;

  
  assert(m_nphenotypes > 0);

  // Generate random flag to turn on
  CellFlag flag;
  flag.setFlagOn(m_distrib(m_gen));
  cell.pflag = flag.toBase10();

  return WRITE_CELL;
}

int CountProcessor::ProcessHeader(CellHeader& header) {
  return HEADER_NO_ACTION; // do nothing
}

int CountProcessor::ProcessLine(Cell& cell) {
  m_count++;
  return NO_WRITE_CELL; // do nothing
}

void CountProcessor::PrintCount() {
  std::cout << m_count << std::endl;
}

int ScatterProcessor::ProcessHeader(CellHeader& header) {
  m_header = header;
  
  m_header.addTag(Tag(Tag::PG_TAG, "", m_cmd));
  m_header.SortTags();

  // setup random number generator
  m_gen = std::mt19937(m_seed);
  m_wdistrib = std::uniform_int_distribution<>(0, m_width);
  m_hdistrib = std::uniform_int_distribution<>(0, m_height);  

  assert(m_width > 0);
  assert(m_height > 0);
  
  // just in time output
  this->SetupOutputStream();
  
  // output the header
  assert(m_archive);
  (*m_archive)(m_header);

  return HEADER_NO_ACTION; 

}

int ScatterProcessor::ProcessLine(Cell& cell) {

  cell.x = m_wdistrib(m_gen);
  cell.y = m_hdistrib(m_gen);  

  return WRITE_CELL;
}

int HeadProcessor::ProcessHeader(CellHeader& header) {
  m_header = header;

  m_header.addTag(Tag(Tag::PG_TAG, "", m_cmd));
  m_header.SortTags();
  
  // just in time output
  this->SetupOutputStream();
  
  // output the header
  assert(m_archive);
  (*m_archive)(m_header);

  return HEADER_NO_ACTION; 
}

int CutProcessor::ProcessHeader(CellHeader& header) {

  m_header = header;

  // find indicies of columns to remove
  size_t i = 0;
  for (const auto& t : header.GetAllTags()) {

    // don't cut non-data tags
    if (t.type != Tag::MA_TAG && t.type != Tag::CA_TAG)
      continue;
    
    if (!m_include.count(t.id)) {
      m_to_remove.insert(i);
    }
    i++;
  }
  
  // cut down the header
  m_header.Cut(m_to_remove);

  m_header.addTag(Tag(Tag::PG_TAG, "", m_cmd));
  m_header.SortTags();

  // just in time output
  this->SetupOutputStream();

  // output the header
  assert(m_archive);
  (*m_archive)(m_header);
  
  return HEADER_NO_ACTION; 
}

int CutProcessor::ProcessLine(Cell& cell) {

  // just transfer to a new set of columns,
  // not including what is to be cut
  std::vector<float> cols_new;
  assert(cell.cols.size() >= m_to_remove.size());
  cols_new.reserve(cell.cols.size() - m_to_remove.size());
  
  for (size_t i = 0; i < cell.cols.size(); i++) {
    if (!m_to_remove.count(i)) {
      cols_new.push_back(cell.cols.at(i));
    }
  }
  cell.cols = cols_new;
  
  return WRITE_CELL;
}

int HeadProcessor::ProcessLine(Cell& cell) {

  m_current_n++;
  if (m_current_n <= m_n) {
    return WRITE_CELL;
  }
  return NO_WRITE_CELL;
}

int DivideProcessor::ProcessHeader(CellHeader& header) {
  
  m_header = header;

  const std::string div_ans = m_numer_string + "div" + m_denom_string;
  
  // error handling
  if (m_numer_string == m_denom_string) {
    throw std::runtime_error("Error: Numerator and denominator must be different");
  }

  // find the indicies
  size_t i = 0;
  for (const auto& t : header.GetAllTags()) { 
    if (t.id == m_numer_string) {
      m_numer = i;
    } else if (t.id == m_denom_string) {
      m_denom = i;
    } else if (t.id == div_ans) {
      std::cerr << "Warning: Overwriting existing column " << div_ans << std::endl;
      m_existing_column = i;
    }
    i++;
  }

  // make sure we found something
  if (m_numer < 0) 
    throw std::runtime_error("Error: Numerator column " + m_numer_string + " not in cell table");
  if (m_denom < 0) 
    throw std::runtime_error("Error: Denominator column " + m_denom_string + " not in cell table");

  // add cmd and new tag
  m_header.addTag(Tag(Tag::PG_TAG, "", m_cmd));
  
  m_header.addTag(Tag(Tag::CA_TAG, div_ans, ""));
  
  m_header.SortTags();
  
  // just in time output
  this->SetupOutputStream();
  
  // output the header
  assert(m_archive);
  (*m_archive)(m_header);

  return HEADER_NO_ACTION;
  
}

int DivideProcessor::ProcessLine(Cell& cell) {

  m_count++;

  assert(m_numer >= 0 && m_denom >= 0);
  assert(m_existing_column == -1 || m_existing_column < cell.cols.size());
  assert(m_numer < cell.cols.size() && m_denom < cell.cols.size());

  // divide by zero
  if (cell.cols.at(m_denom) == 0) {
    if (m_existing_column >= 0)
      cell.cols[m_existing_column] = m_div_zero;
    else
      cell.cols.push_back(m_div_zero);
    return WRITE_CELL;
  }
  
  float val = cell.cols.at(m_numer) / cell.cols.at(m_denom);

  // store the value
  if (m_existing_column >= 0)
    cell.cols[m_existing_column] = val;
  else
    cell.cols.push_back(val);

  return WRITE_CELL;
  
}

int LogProcessor::ProcessHeader(CellHeader& header) {

  m_header = header;

  // setup which are marker indicies
  size_t i = 0;

  for (const auto& t : header.GetAllTags()) { 
    if (t.type == Tag::MA_TAG) {
      m_to_log.insert(i);
    }
    i++;
  }
  
  m_header.addTag(Tag(Tag::PG_TAG, "", m_cmd));
  m_header.SortTags();
  
  // just in time output
  this->SetupOutputStream();
  
  // output the header
  assert(m_archive);
  (*m_archive)(m_header);

  return HEADER_NO_ACTION;
}

int LogProcessor::ProcessLine(Cell& cell) {

  m_count++;
  for (size_t i = 0; i < cell.cols.size(); i++) {
    if (m_to_log.count(i)) {
      if (cell.cols[i] > 0) {
	cell.cols[i] = std::log10(cell.cols.at(i));
      } else {
	
	if (!m_bool_warning_emitted) {
	  std::vector<Tag> tags = m_header.GetDataTags();
	  m_bool_warning_emitted = true;
	  std::cerr << "Warning: encountered zero or negative number to log on " <<
	    "line " << m_count << " with value " << cell.cols[i] << " on column " <<
	    tags.at(i).id << " - removing this line in output. This warning " <<
	    "will emit only once per file" << std::endl;
	}
	
	return NO_WRITE_CELL;
      }
    }
  }
    
  return WRITE_CELL;
  
}

int ROIProcessor::ProcessHeader(CellHeader& header) {

  m_header = header;

  m_header.addTag(Tag(Tag::PG_TAG, "", m_cmd));

  Tag roi_tag(Tag::CA_TAG,"roi", "");
  m_header.addTag(roi_tag);
  
  m_header.SortTags();

  // just in time, make the output stream
  this->SetupOutputStream();
  
  // output the header
  assert(m_archive);
  (*m_archive)(m_header);

  return HEADER_NO_ACTION;
}

int CleanProcessor::ProcessHeader(CellHeader& header) {

  m_header = header;

  // cut the header
  std::unordered_set<size_t> to_remove;
  for (size_t i = 0; i < m_header.size(); i++) {
    //std::cerr << " i " << i << " tag " << m_header.at(i).id << " typ " <<
    //  static_cast<unsigned int>(m_header.at(i).type) << std::endl;
    if ( (m_clean_meta   && m_header.at(i).type == Tag::CA_TAG) ||
         (m_clean_marker && m_header.at(i).type == Tag::MA_TAG) ||
	 (m_clean_graph  && m_header.at(i).type == Tag::GA_TAG))
      to_remove.insert(i);
  }
  m_header.Cut(to_remove);

  // which cells to cut
  const std::vector<Tag> orig_data_tags = header.GetDataTags();
  for (size_t i = 0; i < orig_data_tags.size(); i++) {
    if ( (m_clean_meta   && orig_data_tags.at(i).type == Tag::CA_TAG) ||
         (m_clean_marker && orig_data_tags.at(i).type == Tag::MA_TAG))
      m_to_remove.insert(i);
  }

  // just in time, make the output stream
  this->SetupOutputStream();

  m_header.addTag(Tag(Tag::PG_TAG, "", m_cmd));
  m_header.SortTags();
  
  // output the header
  assert(m_archive);
  (*m_archive)(m_header);

  return HEADER_NO_ACTION;

}

int CleanProcessor::ProcessLine(Cell& cell) {

  // just transfer to a new set of columns,
  // not including what is to be cut
  std::vector<float> cols_new;
  assert(cell.cols.size() >= m_to_remove.size());
  cols_new.reserve(cell.cols.size() - m_to_remove.size());
  
  for (size_t i = 0; i < cell.cols.size(); i++) {
    if (m_to_remove.count(i) == 0) {
      cols_new.push_back(cell.cols.at(i));
    }
  }
  cell.cols = cols_new;

  if (m_clean_cflags)
    cell.cflag = 0;
  if (m_clean_pflags)
    cell.pflag = 0;
  
  // clean the graph
  /*  if (m_clean_graph) {
    cell.m_spatial_ids.clear();
    cell.m_spatial_dist.clear();
    cell.m_spatial_flags.clear();
    }*/

  return 1;
  
}

int ROIProcessor::ProcessLine(Cell& cell) {
  //int ROIProcessor::ProcessLine(const std::string& line) {

  /*
  // convert the line to CellRow
  CellRow values(m_col_count); 
  int num_elems = read_one_line_to_cellrow(line, values, m_header);
  if (num_elems != m_col_count) {
    throw std::runtime_error("Error on line: " + line + "\n\tread " +
			     std::to_string(num_elems) + " expected " +
			     std::to_string(m_col_count));
  }
  
  // get x and y from the string
  float x_, y_;
  get_two_elements_as_floats(line, x_i, y_i, x_, y_);
  
  // Loop the table and check if the cell is in the ROI
  bool print_line = false;
  
  // Loop through all polygons and check if the point is inside any of them
  for (const auto &polygon : m_rois) {
    
    // if point is in this polygon, add the polygon id number to the roi
    if (polygon.PointIn(x_,y_)) {
      
      print_line = true;
      
      //std::cerr << " ADDING POINT " << x_ << "," << y_ << std::endl;
      //new_data->SetNumericElem(polygon.Id, i);
      
      // uncomment below if want to prevent over-writing existing
      // break;
    }
  }
  
  if (print_line)
    std::cout << line << std::endl;

  */
  return 1;
  
}

int PhenoProcessor::ProcessHeader(CellHeader& header) {

  m_header = header;

  m_header.addTag(Tag(Tag::PG_TAG, "", m_cmd));
  
  // build up a map of the indices of the markers
  // in the Cell
  size_t i = 0;
  for (const auto& t : header.GetDataTags()) {
    if (t.type == Tag::MA_TAG)
      m_marker_map[t.id] = i;
    i++;
  }

  // loop through once and just warn about missing data
  for (const auto& b : m_p) {
    if (m_marker_map.find(b.first) == m_marker_map.end()) {
      std::cerr << "Warning: Marker in phenotype file " <<
	b.first << " is not in the header" << std::endl;
    }
  }

  for (const auto& m : m_marker_map) {
    if (m_p.find(m.first) == m_p.end()) {
      std::cerr << "Warning: Marker in cell table " <<
	m.first << " is not in the phenotype file. Bit will be OFF" << std::endl;
    }
  }
  
  // just in time output, so as not to write an empty file if the input crashes
  // set the output to file or stdout
  this->SetupOutputStream(); 
  
  m_header.SortTags();

  // output the header
  assert(m_archive);
  (*m_archive)(m_header);
  
  return HEADER_NO_ACTION;
}

int PhenoProcessor::ProcessLine(Cell& cell) {

  // make a copy of the cell for output
  Cell cell_new = cell;
  
  // initialize an empty flag
  CellFlag flag;
  
  // loop through the gates
  for (const auto& b : m_p) {

    // skip if don't have the marker
    auto m = m_marker_map.find(b.first);
    if (m == m_marker_map.end()) {
      continue;
    }

    // set the index of where
    size_t marker_index = m->second;
    
    // set the flag on if it clears the gates
    if (cell.cols.at(m->second) >= b.second.first &&
	cell.cols.at(m->second) <= b.second.second) {
      //      std::cerr << "Marker: " << b.first << " cleared" << std::endl;
      //std::cerr << "...before " << flag << std::endl; 
      flag.setFlagOn(marker_index);
    }
    
  }

  // convert to cy_uint for storage
  cell.pflag = flag.toBase10();
  
  return 1;
}

int ViewProcessor::ProcessHeader(CellHeader& header) {

  m_header = header;

  std::unordered_set<size_t> to_remove;
  
  // find indicies of columns to display
  size_t i = 0;
  if (m_to_view.size()) {
    for (const auto& t : header.GetAllTags()) {
      
      // don't cut non-data tags
      if (t.type != Tag::MA_TAG && t.type != Tag::CA_TAG)
	continue;
      
      if (m_to_view.count(t.id)) {
	m_to_view_indicies.insert(i);
      } else {
	to_remove.insert(i);
      }
      i++;
    }
    
    // remove
    m_header.Cut(to_remove);
  }

  
  // print it, that's all
  if (m_print_header || m_header_only) 
    m_header.Print();
  
  return m_header_only ? ONLY_WRITE_HEADER : HEADER_NO_ACTION;
  
}

int ViewProcessor::ProcessLine(Cell& cell) {

  // classic view, no cut
  if (!m_to_view.size())  {
    cell.Print(m_round);
    return NO_WRITE_CELL; // don't output, since already printing it
  }
  
  // just transfer to a new set of columns,
  // not including what is to be cut
  std::vector<float> cols_new;
  assert(cell.cols.size() >= m_to_view.size());
  cols_new.reserve(m_to_view.size());

  for (size_t i = 0; i < cell.cols.size(); i++) {
    if (m_to_view_indicies.count(i)) {
      cols_new.push_back(cell.cols.at(i));
    }
  }
  cell.cols = cols_new;
  cell.Print(m_round);
    
  return NO_WRITE_CELL; // don't output, since already printing it
}

int CatProcessor::ProcessLine(Cell& line) {

  // update the cell id
  size_t n_cell_id = m_offset + line.id; 
  if (n_cell_id > m_max_cellid)
    m_max_cellid = n_cell_id;
  line.id = n_cell_id;
  
  /*  if (line.m_spatial_ids.size() && !m_error_emitted) {
    std::cerr << "Warning: Graph concatenation not supported" << std::endl;
    m_error_emitted = true;
    }*/
  
  // update the graph ids
  /*  for (const auto& i : m_graph_indicies) {
    const std::string& graph_line = tokens.at(i);

    // parse the node and rest the cell-ids with the offset
    CellNode node(graph_line);
    node.OffsetNodes(m_offset);

    tokens[i] = node.toString(false); // false is for "integerize"
    }*/
  
  /*
  // update the cell id
  tokens[m_cellid_index] = std::to_string(n_cell_id);

  // add the sample number
  tokens.emplace_back(std::to_string(m_sample));
  
  // output
  std::cout << tokens_to_comma_string(tokens) << std::endl;

  */
  return WRITE_CELL;
}

int CatProcessor::ProcessHeader(CellHeader& header) {

  // if we already have the master header, just compare for error checking
  if (m_master_set) {

    // check compatibilty. The checker will emit the error message
    if (!m_header.isConcatenatable(header)) {
      assert(false);
    }
    
  } else {
    
    m_header = header;
    m_master_set = true;
    m_header.addTag(Tag(Tag::PG_TAG, "", m_cmd));

    // just in time output, so as not to write an empty file if the input crashes
    // set the output to file or stdout
    this->SetupOutputStream();
    
    m_header.SortTags();
    
    // write the header
    (*m_archive)(m_header);
    
  }
  
  // add the sample tag if not already there
  /*  if (!m_header.hasTag("sample")) {
    Tag sample_tag("CA","sample");
    m_header.addTag(sample_tag);
  } else {
    std::cerr << "Warning: Header already has sample tag, overwriting" << std::endl;
    }*/

  return HEADER_NO_ACTION;
}

int CerealProcessor::ProcessHeader(CellHeader& header) {
  m_header = header;

  m_header.addTag(Tag(Tag::PG_TAG, "", m_cmd));
  
  assert(!m_filename.empty());

  // set the output to file or stdout
  if (m_filename == "-") {
    m_archive = std::make_unique<cereal::PortableBinaryOutputArchive>(std::cout);
  } else {
    m_os = std::make_unique<std::ofstream>(m_filename, std::ios::binary);
    m_archive = std::make_unique<cereal::PortableBinaryOutputArchive>(*m_os);
  }

  m_header.SortTags();
  
  // archive the header
  (*m_archive)(m_header);
  
  return 0; // minor, but the HEADER_NO_ACTION is in CellProcessor abstract class, of which CerealProcessor is uniquely not a child
}

int CerealProcessor::ProcessLine(const std::string& line) {

  Cell row(line, m_header, m_cellid, m_sampleid);
  m_cellid++;

  // serialize it
  (*m_archive)(row);

  // serializing here, so don't need to write in the cell_table call
  // minor, but NO_WRITE_CELL is in CellProcessor abstract class, of which CerealProcessor is uniquely not a child  
  return 0; 
  
}

int DebugProcessor::ProcessHeader(CellHeader& header) {

  m_header = header;
  
  // just in time, make the output stream
  this->SetupOutputStream();

  //m_header.addTag(Tag(Tag::PG_TAG, "", m_cmd));
  m_header.SortTags();
  
  // output the header
  assert(m_archive);
  (*m_archive)(m_header);

  return HEADER_NO_ACTION;
  
}

int DebugProcessor::ProcessLine(Cell& cell) {

  cell.set_cell_id(m_cell_id);
  m_cell_id++;
  return WRITE_CELL;
  
}

int BuildProcessor::ProcessHeader(CellHeader& header) {
  m_header = header; // store but don't print

  return 2;
}

// just a pass through to place it in the table
int BuildProcessor::ProcessLine(Cell& cell) {

  return SAVE_CELL;
  
}
