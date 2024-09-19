#include "cell_processor.h"
#include "polygon.h"
#include "cell_utils.h"
#include "cell_flag.h"
#include "cell_selector.h"

#include "cell_row.h"

#include <regex>
#include <limits>

static inline bool roikey(const Polygon& poly, const std::string& keyword) {
    return poly.Text.find(keyword) != std::string::npos ||
           poly.Name.find(keyword) != std::string::npos;
}

int FilterProcessor::ProcessHeader(CellHeader& header) {
  
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

  if (m_verbose) 
  for (const auto& m : m_criteria_int) {
    std::cerr <<   "PH --- Marker/Meta field: " << m.first << std::endl;
    for (const auto& v : m.second) {
      std::cerr << "           - " << static_cast<int>(v.first) << " - value " << v.second << std::endl;     
    }
  }
  
  if (m_criteria.size() != m_criteria_int.size())
    throw std::runtime_error("cyftools filter -- Unable to find all requested fields in the header");

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

  // check that group by is in the header
  if (!m_group_by.empty()) {
    size_t i = 0;
    for (const auto& tag : m_header.GetDataTags()) {
      if (tag.id == m_group_by) {
	m_group_by_id = i;
	continue;
      }
      i++;
    }
    if (m_group_by_id < 0) {
      std::cerr << "Error: cyftools mean - group by of " << m_group_by << " not found in data header" << std::endl;
      assert(false);
    }
  }
  
  // which column (if any) stores tls-id
  /*size_t i = 0;
  for (const auto& tag : m_header.GetDataTags()) {
    if (tag.id == "tls_id") {
      m_tls_column = i;
      break;
    }
    i++;
    }*/
  
  // initialize sums to zero
  //sums.resize(m_header.GetDataTags().size());
  //if (sums.size())
  //  assert(sums.at(0) == 0);
  
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
    
    // get the group-by key, unless no group then all are ZERO key
    int key = m_group_by_id < 0 ? 0 : cell.cols.at(m_group_by_id);

    // add the data point
    allgroups[key].Add(i, cell.cols.at(i));
  }

  // if tls is a column, add ID to set
  //if (m_tls_column > 0)
  //  m_tls_set.insert(cell.cols.at(m_tls_column));
  
  // don't emit anything in StreamTable
  return CellProcessor::NO_WRITE_CELL;
}

void AverageProcessor::EmitCells() const {

  // loop each group and emit the cell with the means
  for (const auto& d : allgroups) {
    OutputLine(d.second.EmitCell());
  }
  
  // for tls id column, we want the number of unique TLS
  //cell.cols[m_tls_column] = m_tls_set.size();
  
  // write the one cell
  //OutputLine(cell);
}

int CellCountProcessor::ProcessHeader(CellHeader& header) {
  
  m_header = header;

  // initialize the count vector
  m_num_marker_tags = header.GetMarkerTags().size();
  m_counts = std::vector<size_t>(m_num_marker_tags + m_additional_flags.size());

  m_header.ClearMeta();

  // add single markers
  for (const auto& m : header.GetMarkerTags()) {
    m_header.addTag(Tag(Tag::CA_TAG, m.id + "_count", ""));
  }

  // add additional flag markers
  for (const auto& m : m_additional_flags) {
    m_header.addTag(Tag(Tag::CA_TAG, std::to_string(m) + "_count", ""));
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

  for (size_t i = 0; i < m_num_marker_tags; i++) {
    CellFlag f(cell.pflag);
    cy_uint to_test = static_cast<cy_uint>(std::pow(2, i));
    assert(to_test > 0);
    if (f.testAndOr(to_test, 0))
      m_counts[i]++;
  }

  // additional flags
  for (size_t i = m_num_marker_tags; i < (m_num_marker_tags + m_additional_flags.size()); i++) {
    if (m_additional_flags[i-m_num_marker_tags] == 0) {
      m_counts[i] += (cell.pflag == 0);
    } else if (IS_FLAG_SET(cell.pflag, m_additional_flags[i - m_num_marker_tags])) {
      m_counts[i]++;
    }
  }
  
  return CellProcessor::NO_WRITE_CELL;
}

void CellCountProcessor::EmitCell() const {

  Cell cell;

  // add markers as dummies
  for (size_t i = 0; i < m_num_marker_tags; i++)
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

int TrimProcessor::ProcessHeader(CellHeader& header) {
  m_header = header;

  m_header.addTag(Tag(Tag::PG_TAG, "", m_cmd));

  m_header.SortTags();
  
  // just in time, make the output stream
  this->SetupOutputStream();
  
  // output the header
  assert(m_archive);
  (*m_archive)(m_header);
  
  return HEADER_NO_ACTION;

}

int TrimProcessor::ProcessLine(Cell& cell) {
  
  if (cell.cflag & MARK_FLAG) 
    return WRITE_CELL;
  return NO_WRITE_CELL;
}

int FilterProcessor::ProcessLine(Cell& cell) {
  
  bool write_cell = false;

  ///////
  // FLAGS
  ///////
  // NB: even if flags are all empty, default should be to trigger a "write_cell = true"

  //std::cerr << cell.cflag << " - " << cell.pflag << std::endl << m_select << std::endl;
  
  // if C- or P-flags met, print the cell
  if ( m_select.TestFlags(cell.pflag, cell.cflag)) { 
    write_cell = true;
  }

  // if only selecting on field, start with default of "true"
  if (m_select.size() == 0) {
    write_cell = true;
  }

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


      //std::cerr << "B write cell " << write_cell << " c.first " << c.first << " cc.second (val) " << cc.second << " m or toggle " << m_or_toggle << std::endl;	 
      switch (cc.first) { // switching on the optype
      case optype::GREATER_THAN: write_cell          = write_cell = m_or_toggle ? (write_cell || (value >  cc.second)) : (write_cell && (value >  cc.second));  break;
      case optype::LESS_THAN: write_cell             = write_cell = m_or_toggle ? (write_cell || (value <  cc.second)) : (write_cell && (value <  cc.second));  break;
      case optype::GREATER_THAN_OR_EQUAL: write_cell = write_cell = m_or_toggle ? (write_cell || (value >= cc.second)) : (write_cell && (value >= cc.second));  break;
      case optype::LESS_THAN_OR_EQUAL: write_cell    = write_cell = m_or_toggle ? (write_cell || (value <= cc.second)) : (write_cell && (value <= cc.second));  break;
      case optype::EQUAL_TO: write_cell              = write_cell = m_or_toggle ? (write_cell || (value == cc.second)) : (write_cell && (value == cc.second));  break;
      default: assert(false);
      }

      //std::cerr << "A write cell " << write_cell << " c.first " << c.first << " cc.second (val) " << cc.second << " m or toggle " << m_or_toggle << std::endl;
    }
    
  }

  write_cell = write_cell && flag_write_cell;

  bool m_trim = !(m_mark1 || m_mark2);

  //  std::cerr << " m_mark1 " << m_mark1 << " m_mark2 " << m_mark2 << " mtrim " << m_trim << " write_cell " << write_cell << " cflag " <<
  //  cell.cflag << std::endl;
  
  // trim is set but cell passes
  if (m_trim && write_cell) {
    CLEAR_FLAG(cell.cflag, MARK_FLAG);
    return CellProcessor::WRITE_CELL;
    // trim is set, but cell doesn't pass       
  } else if (m_trim && !write_cell) {
    return CellProcessor::NO_WRITE_CELL;
    // cell passes, mark as 1
  } else if (m_mark1 && write_cell) {
    SET_FLAG(cell.cflag, MARK_FLAG);
    return CellProcessor::WRITE_CELL;
    // cell passes, mark as 2
  } else if (m_mark2 && write_cell) {
    SET_FLAG(cell.cflag, BUILD_GRAPH_FLAG);
    return CellProcessor::WRITE_CELL;
  // cell fails, but mode is mark not trim
  } else if (!m_trim && !write_cell) {
    if (m_mark1)
      CLEAR_FLAG(cell.cflag, MARK_FLAG);
    if (m_mark2)
      CLEAR_FLAG(cell.cflag, BUILD_GRAPH_FLAG);
    return CellProcessor::WRITE_CELL;
  // should never get here
  } else {
    assert(false);
  }

  // should never get here
  assert(false);
  return CellProcessor::WRITE_CELL; // don't write if not selected
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
  
  //  if (IS_FLAG_SET(cell.cflag, m_c_and_flags) && IS_FLAG_SET(cell.pflag, m_p_and_flags) &&
  //    ARE_FLAGS_OFF(cell.cflag, m_c_not_flags) && ARE_FLAGS_OFF(cell.pflag, m_p_not_flags))
    m_count++;
  
  return NO_WRITE_CELL; // do nothing
}

void CountProcessor::PrintCount() {
  std::cout << m_count << std::endl;
}

int MagnifyProcessor::ProcessHeader(CellHeader& header) {
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

int RescaleProcessor::ProcessHeader(CellHeader& header) {
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

int RescaleProcessor::ProcessLine(Cell& cell) {
  
  size_t i = 0;
  size_t marker_i = 0;
  for (auto& m : m_header.GetDataTags()) {
    if (m.type == Tag::MA_TAG) {
      if (!IS_FLAG_I_SET(cell.pflag, marker_i)) {
	cell.cols[i] = 0;
      }
      marker_i++;
    }
    i++;
  }
  
  return WRITE_CELL;
}


int MagnifyProcessor::ProcessLine(Cell& cell) {

  cell.x = cell.x * m_factor;
  cell.y = cell.y * m_factor;
  return WRITE_CELL;
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

  // find indicies of columns to keep
  size_t i = 0; // all tag index
  std::unordered_set<size_t> m_to_remove;  
  size_t d = 0; // data tag index

  for (const auto& t : header.GetAllTags()) {

    // skip non-data
    if (!t.isData()) {
      i++;
      continue;
    }

    // remove meta/marker not in list (for header Cut)
    if (!m_include.count(t.id)) 
      m_to_remove.insert(i);
    else
      // store the data index, for lookup in cell table
      m_data_to_keep.insert(d);

    // update iterators
    d++;
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
  assert(cell.cols.size() >= m_data_to_keep.size());
  cols_new.reserve(m_data_to_keep.size());
  
  for (size_t i = 0; i < cell.cols.size(); i++) {
    if (m_data_to_keep.count(i)) {
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

int FlagsetProcessor::ProcessHeader(CellHeader& header) {
  
  m_header = header;

  // add cmd and sort
  m_header.addTag(Tag(Tag::PG_TAG, "", m_cmd));
  m_header.SortTags();
  
  // just in time output
  this->SetupOutputStream();
  
  // output the header
  assert(m_archive);
  (*m_archive)(m_header);

  return HEADER_NO_ACTION;
  
}

int FlagsetProcessor::ProcessLine(Cell& cell) {

  if (IS_FLAG_SET(cell.cflag, MARK_FLAG)) {

    // set the pflag
    if (m_pflag_to_set) {
      SET_FLAG(cell.pflag, m_pflag_to_set);
    }

    // set the cflag
    if (m_cflag_to_set) {
      SET_FLAG(cell.cflag, m_cflag_to_set);      
    }

    // clear the pflag
    if (m_pflag_to_clear) {
      CLEAR_FLAG(cell.pflag, m_pflag_to_clear);
    }

    // clear the cflag
    if (m_cflag_to_clear) {
      CLEAR_FLAG(cell.cflag, m_cflag_to_clear);
    }

  }

  // no matter what, clear the mark
  CLEAR_FLAG(cell.cflag, MARK_FLAG);

  return CellProcessor::WRITE_CELL;
  
}

int MarkCheckProcessor::ProcessLine(Cell& cell) {

  if (IS_FLAG_SET(cell.cflag, MARK_FLAG)) {
    std::cerr << "MARK_FLAG is set on cell " << cell << std::endl;
    assert(false);    
  }

  return CellProcessor::WRITE_CELL;
}

// just a pass through, do absolutely nothing but stream out header
int MarkCheckProcessor::ProcessHeader(CellHeader& header) {

  m_header = header;
  
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

  // check if ROI regions are in the header
  /*for (const auto &polygon : m_rois) {
    if (roikey(polygon, "region")) {
      Tag roi_tag(Tag::CA_TAG,"roi_region", "");
      m_header.addTag(roi_tag);
      m_roi_region = true;
    }
    }*/
  
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
         (m_clean_marker && m_header.at(i).type == Tag::MA_TAG))
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

  if (m_clean_programs)
    m_header.CleanProgramTags();
  
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
  
  return 1;
  
}


int ROIProcessor::ProcessLine(Cell& cell) {

  // Loop the table and check if the cell is in the ROI
  bool print_line = true; //m_blacklist_remove;

  // Loop through all polygons and check if the point is inside any of them
  for (const auto &polygon : m_rois) {
    
    int number = 0; // region number
    
    // if point is in this polygon, add the polygon id number to the roi
    if (polygon.PointIn(cell.x,cell.y)) {
      
      if (roikey(polygon, "lacklist") || roikey(polygon, "rtifact")) {
	print_line = false;
      } else if (roikey(polygon, "ormal")) {
	CLEAR_FLAG(cell.cflag, TUMOR_FLAG);
	CLEAR_FLAG(cell.cflag, MARGIN_FLAG);	
	CLEAR_FLAG(cell.cflag, TUMOR_MANUAL_FLAG);
	print_line = true;
      } else if (roikey(polygon, "umor")) {
	SET_FLAG(cell.cflag, TUMOR_MANUAL_FLAG);
	print_line = true;
      } else if (roikey(polygon, "cd3panck_error")) {
	if (IS_FLAG_SET(cell.pflag, ORION_PANCK))
	  CLEAR_FLAG(cell.pflag, ORION_CD3);
	print_line = true;
      } else if (roikey(polygon, "region") && false) {
	std::regex pattern(R"(region(\d+))");
	std::smatch matches;
	// Check if the input string matches the pattern
	if (std::regex_match(polygon.Text, matches, pattern)) {
	  // Extract the trailing integer part and convert to an integer
	  number = std::stoi(matches[1].str());
	  assert(m_roi_region);
	} else if (std::regex_match(polygon.Name, matches, pattern)) {
	  number = std::stoi(matches[1].str());
	  assert(m_roi_region);	  
	} else {
	  std::cerr << "Warning: ROI with name " << polygon.Text << "-" <<
	    polygon.Name << " cannot be parsed" << std::endl;
	}
      } else {
	print_line = true;
      }

      // add the region number if indicated
      if (m_roi_region)
	cell.cols.push_back(number);
      
      // secondary annotations for PCa
      if (roikey(polygon, "3+3")) {
	SET_FLAG(cell.cflag, GLEASON_GRADE_GROUP_1);
	SET_FLAG(cell.cflag, TUMOR_MANUAL_FLAG);	
      } else if (roikey(polygon, "3+4")) {
	SET_FLAG(cell.cflag, GLEASON_GRADE_GROUP_2);
	SET_FLAG(cell.cflag, TUMOR_MANUAL_FLAG);	
      } else if (roikey(polygon, "4+3")) {
	SET_FLAG(cell.cflag, GLEASON_GRADE_GROUP_3);
	SET_FLAG(cell.cflag, TUMOR_MANUAL_FLAG);	
      } else if (roikey(polygon, "4+4")) {
	SET_FLAG(cell.cflag, GLEASON_GRADE_GROUP_4);
	SET_FLAG(cell.cflag, TUMOR_MANUAL_FLAG);	
      } else if (roikey(polygon, "5+4") || roikey(polygon, "4+5") || roikey(polygon, "5+5")) {
	SET_FLAG(cell.cflag, GLEASON_GRADE_GROUP_5);
	SET_FLAG(cell.cflag, TUMOR_MANUAL_FLAG);	
      } else if (roikey(polygon, "eural") || roikey(polygon, "PNI") || roikey(polygon, "pni")) {
	SET_FLAG(cell.cflag, PERINEURAL_INVASION);
      } else if (roikey(polygon, "SV") || roikey(polygon, "eminal")) {
	SET_FLAG(cell.cflag, SEMINAL_VESICLES);
      }
      // uncomment below if want to prevent over-writing existing
      // break;
    }
  }
  
  if (print_line)
    return CellProcessor::WRITE_CELL;
    
  return CellProcessor::NO_WRITE_CELL;
  
}

int PhenoProcessor::ProcessHeader(CellHeader& header) {

  m_header = header;

  m_header.addTag(Tag(Tag::PG_TAG, "", m_cmd));
  
  // if randomly scaling gates
  if (m_random_scale > 0) {

    // random number generator
    std::random_device rd;
    std::mt19937 m_gen(rd());
    std::uniform_real_distribution<> dis(1 - m_random_scale, 1 + m_random_scale);

    // adjust the new gate
    std::string pheno_tag;    
    for (auto& p : m_p) {
      float random_gate_scale = dis(m_gen);
      int new_low_gate = static_cast<int>(p.second.first * random_gate_scale);
      pheno_tag += p.first + ":" + std::to_string(new_low_gate) + ",";
      p.second.first  = new_low_gate;
    }
    if (!pheno_tag.empty()) //remove last comma
      pheno_tag.pop_back();

    // add the tag
    m_header.addTag(Tag(Tag::PG_TAG, "", pheno_tag));
  }
  
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
      std::cerr << "Warning: Marker " <<
	std::left << std::setw(12) << std::setfill(' ') <<
	m.first << " not in phenotype file. Bit to OFF" << std::endl;
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
    if (cell.cols.at(m->second) >= b.second.first*m_scale &&
	cell.cols.at(m->second) <= b.second.second*m_scale) {
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
  
  // find indices of columns to display
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
  if (m_print_header || m_header_only) {
    m_header.Print();
  // or print as csv version 
  } else if (m_csv_header) {
    std::cout << "sid,cid,cflag,pflag,x,y,"; 
    const auto& tags = m_header.GetDataTags();
    for (auto it = tags.begin(); it != tags.end(); ++it) {
      std::cout << it->id;
      if (std::next(it) != tags.end()) {
        std::cout << ",";
      }
    }
    std::cout << std::endl;
  } else if (m_crevasse) { // print as crevasse format (csv header)
    std::cout << "CellID,X,Y,";
    const auto& tags = m_header.GetDataTags();
    for (auto it = tags.begin(); it != tags.end(); ++it) {
      if (it->type == Tag::MA_TAG) {
	std::cout << it->id;
	if (std::next(it) != tags.end() && std::next(it)->type == Tag::MA_TAG) {
	  std::cout << ","; 
	}
      }
    }
    std::cout << std::endl;
  } else if (m_adjacent) {
    ; // no print
  }
  
  return m_header_only ? ONLY_WRITE_HEADER : HEADER_NO_ACTION;
  
}

int OffsetProcessor::ProcessLine(Cell& cell) {
  
  cell.x = cell.x + m_x;
  cell.y = cell.y + m_y;

  return WRITE_CELL;
}

int OffsetProcessor::ProcessHeader(CellHeader& header) {

  m_header = header;
  
  m_header.addTag(Tag(Tag::PG_TAG, "", m_cmd));
  
  this->SetupOutputStream(); 
  
    m_header.SortTags();
  
  // output the header
  assert(m_archive);
  (*m_archive)(m_header);
  
  return HEADER_NO_ACTION;
  
}

int CheckProcessor::ProcessHeader(CellHeader& header) {

  m_header = header;

  this->SetupOutputStream();
  
  // output the header
  assert(m_archive);
  (*m_archive)(m_header);
  
  return HEADER_NO_ACTION;

}

int CheckProcessor::ProcessLine(Cell& cell) {
  return WRITE_CELL;
}

int FlipProcessor::ProcessHeader(CellHeader& header) {

  m_header = header;

  // adds tag to hold the command input in the header to keep track of what was done 
  m_header.addTag(Tag(Tag::PG_TAG, "", m_cmd));
  
  this->SetupOutputStream(); 
  
  m_header.SortTags();
  
  // output the header
  assert(m_archive);
  (*m_archive)(m_header);
  
  return HEADER_NO_ACTION;
  
}

int FlipProcessor::ProcessLine(Cell& cell) {

  // write the flip logic
  if (m_x != -100000)
    cell.x = 2 * m_x - cell.x + m_xmax;
  if (m_y != -100000)  
    cell.y = 2 * m_y - cell.y + m_ymax;  

  return WRITE_CELL;
}


int ViewProcessor::ProcessLine(Cell& cell) {

  // classic view, no cut
  if (!m_to_view.size())  {
    if (m_crevasse) {
      cell.PrintForCrevasse(m_header);
      return NO_WRITE_CELL;
    }
    
    cell.PrintWithHeader(m_round, m_tabprint, m_adjacent, m_header, false);

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

  cell.PrintWithHeader(m_round, m_tabprint, m_adjacent, m_header, m_strict_cut);
    
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
  
  Cell row(line,
	   m_id_index,
	   m_x_index,
	   m_y_index,
	   m_header,
	   m_cellid,
	   m_sampleid);

  // used here only if m_id_index < 0
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
  
  return CellProcessor::SAVE_HEADER;
}

// just a pass through to place it in the table
int BuildProcessor::ProcessLine(Cell& cell) {

  return SAVE_CELL;
  
}
