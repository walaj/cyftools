#include "cell_processor.h"
#include "polygon.h"
#include "cell_utils.h"
#include "cell_flag.h"
#include "cell_selector.h"
#include "cyf_flags.h"   // standard @FL flag-bit declarations

#include "cell_row.h"

#include <regex>
#include <limits>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <array>
#include <cstdio>
#include <cstdlib>

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

  // clean cflags
  if (m_clean_cflags) {
    cell.cflag = 0;
  } else if (m_cflag > 0) {
    cell.cflag &= ~m_cflag;
  }

  // clean pflags
  if (m_clean_pflags) {
    cell.pflag = 0;
  } else if (m_pflag > 0) {
    cell.pflag &= ~m_pflag;
  }

  return 1;

}

// The header edit, shared by ProcessHeader (full rewrite) and the fast block-copy
// reheader. Returns false if the ROI add policy refuses the operation.
bool AddTagProcessor::computeNewHeader(CellHeader& h) {

  // ROI add/overwrite policy (used by addroi): refuse to silently add onto a file
  // that already has ROIs unless the user picked --add or --overwrite.
  if (m_roi_mode != ROI_NONE) {
    size_t existing = 0;
    for (const auto& t : h) if (t.type == Tag::RO_TAG) ++existing;
    if (existing > 0 && m_roi_mode == ROI_REQUIRE) {
      std::cerr << "cyftools addroi: file already has " << existing
                << " ROI(s). Pass --add to append or --overwrite to replace." << std::endl;
      return false;
    }
    if (existing > 0 && m_roi_mode == ROI_OVERWRITE)
      h.RemoveRoiTags("", -1);
  }

  // ROIs are appended (each polygon is its own @RO; an --add of the same file
  // really appends, even when ids collide). Non-ROI tags merge their fields into
  // an existing same-class/same-id tag (the addtag behavior).
  for (const auto& t : m_tags) {
    if (m_roi_mode != ROI_NONE) h.appendRawTag(t);
    else                        h.UpsertTag(t);
  }

  // optional sample group/category -> @SA GP field. Merge into the file's
  // existing @SA tag (keep its id/name); create one keyed by the default id
  // (the input filename stem) only if the file has no @SA tag yet.
  if (!m_group.empty()) {
    std::string sid = m_group_default_id;
    for (const auto& t : h) if (t.type == Tag::SA_TAG) { sid = t.id; break; }
    h.UpsertTag(Tag(Tag::SA_TAG, sid, std::string("GP:") + m_group));
  }

  h.addTag(Tag(Tag::PG_TAG, "", m_cmd));   // provenance
  h.SortTags();
  return true;
}

int AddTagProcessor::ProcessHeader(CellHeader& header) {
  m_header = header;
  if (!computeNewHeader(m_header)) {
    m_aborted = true;
    return ONLY_WRITE_HEADER;   // stop before opening output -> no file written
  }
  this->SetupOutputStream();    // just in time, make the output stream
  assert(m_archive);
  (*m_archive)(m_header);
  return HEADER_NO_ACTION;
}

int AddTagProcessor::ProcessLine(Cell& cell) {
  return WRITE_CELL;   // pass every cell through unchanged
}

bool ScaleRoiProcessor::computeNewHeader(CellHeader& h) {
  // clean number formatting (no scientific notation for typical coords)
  auto fmtnum = [](double v) {
    char buf[32];
    std::snprintf(buf, sizeof buf, "%.7g", v);
    return std::string(buf);
  };

  // scale every @RO polygon's PT coordinates; collect first, then upsert (so we
  // don't mutate the tag vector while iterating it).
  std::vector<Tag> updates;
  for (const auto& t : h) {
    if (t.type != Tag::RO_TAG) continue;
    std::istringstream ps(t.GetField("PT"));
    std::string pair, out;
    while (ps >> pair) {
      const size_t c = pair.find(',');
      if (c == std::string::npos) continue;
      try {
        double x = std::stod(pair.substr(0, c))  * m_factor;
        double y = std::stod(pair.substr(c + 1)) * m_factor;
        if (m_flipx) x = -x;
        if (m_flipy) y = -y;
        x += m_xoff; y += m_yoff;
        if (!out.empty()) out += " ";
        out += fmtnum(x) + "," + fmtnum(y);
      } catch (...) { /* skip a malformed coordinate pair */ }
    }
    updates.push_back(Tag(Tag::RO_TAG, t.id, std::string("PT:") + out));   // by-value ctor: no ODR-use of RO_TAG
  }
  for (auto& u : updates) h.UpsertTag(u);   // PT field overrides the old one

  h.addTag(Tag(Tag::PG_TAG, "", m_cmd));     // provenance
  h.SortTags();
  return true;
}

int ScaleRoiProcessor::ProcessHeader(CellHeader& header) {
  m_header = header;
  computeNewHeader(m_header);
  this->SetupOutputStream();
  (*m_archive)(m_header);
  return HEADER_NO_ACTION;
}

int ScaleRoiProcessor::ProcessLine(Cell& cell) {
  return WRITE_CELL;   // pass every cell through unchanged
}

bool ClearRoiProcessor::computeNewHeader(CellHeader& h) {
  const size_t n = h.RemoveRoiTags(m_name, m_sample, m_id);
  if (m_verbose)
    std::cerr << "cyftools clearroi: removed " << n << " @RO tag(s)" << std::endl;
  h.addTag(Tag(Tag::PG_TAG, "", m_cmd));     // provenance
  h.SortTags();
  return true;
}

int ClearRoiProcessor::ProcessHeader(CellHeader& header) {
  m_header = header;
  computeNewHeader(m_header);
  this->SetupOutputStream();
  (*m_archive)(m_header);
  return HEADER_NO_ACTION;
}

int ClearRoiProcessor::ProcessLine(Cell& cell) {
  return WRITE_CELL;   // pass every cell through unchanged
}

int ValidateProcessor::ProcessHeader(CellHeader& header) {
  version    = header.GetHeaderField("VN");
  sort_order = header.GetHeaderField("SO");
  mpp        = header.GetHeaderField("MP");
  units      = header.GetHeaderField("UN");
  return ONLY_WRITE_HEADER;   // header is all we need; stop before the cells
}

int RoiInfoProcessor::ProcessHeader(CellHeader& header) {
  mpp   = header.GetHeaderField("MP");
  units = header.GetHeaderField("UN");

  for (const auto& t : header.GetAllTags()) {
    if (t.type != Tag::RO_TAG) continue;

    std::vector<JPoint> verts;
    std::istringstream ps(t.GetField("PT"));
    std::string pair;
    while (ps >> pair) {
      const size_t comma = pair.find(',');
      if (comma == std::string::npos) continue;
      try {
        verts.emplace_back(std::stof(pair.substr(0, comma)),
                           std::stof(pair.substr(comma + 1)));
      } catch (...) { /* skip malformed coordinate */ }
    }

    RoiEntry e;
    e.id     = t.id;
    e.name   = t.GetField("NM");
    e.sample = t.GetField("SA");
    e.nverts = verts.size();
    e.area   = Polygon(verts).Area();
    rois.push_back(std::move(e));
  }

  return ONLY_WRITE_HEADER;   // header is all we need; stop before the cells
}

int FlagRoiProcessor::ProcessHeader(CellHeader& header) {

  m_header = header;

  // build polygons from the @RO header tags (PT = "x,y x,y ..." pairs)
  for (const auto& t : m_header.GetAllTags()) {
    if (t.type != Tag::RO_TAG) continue;

    if (!m_name_filter.empty() && t.GetField("NM").find(m_name_filter) == std::string::npos)
      continue;
    if (m_sample_filter >= 0) {
      const std::string sa = t.GetField("SA");
      if (!sa.empty() && std::stol(sa) != m_sample_filter) continue;
    }

    std::vector<JPoint> verts;
    std::istringstream ps(t.GetField("PT"));
    std::string pair;
    while (ps >> pair) {
      const size_t comma = pair.find(',');
      if (comma == std::string::npos) continue;
      try {
        verts.emplace_back(std::stof(pair.substr(0, comma)),
                           std::stof(pair.substr(comma + 1)));
      } catch (...) { /* skip malformed coordinate */ }
    }
    if (verts.size() >= 3)
      m_polys.emplace_back(verts);
  }

  if (m_polys.empty())
    std::cerr << "Warning: flagroi found no matching @RO polygons in the header" << std::endl;

  this->SetupOutputStream();
  m_header.addTag(Tag(Tag::PG_TAG, "", m_cmd));
  m_header.SortTags();
  assert(m_archive);
  (*m_archive)(m_header);
  return HEADER_NO_ACTION;
}

int FlagRoiProcessor::ProcessLine(Cell& cell) {
  uint64_t& reg = (m_reg == 'c') ? cell.cflag : cell.pflag;

  // overwrite mode: clear the bit on every cell first, so that after the full
  // pass the bit is set on exactly the cells inside the ROIs (and nowhere else).
  if (m_overwrite)
    reg &= ~m_mask;

  for (const auto& poly : m_polys) {
    if (poly.PointIn(cell.x, cell.y)) {
      reg |= m_mask;
      break;   // one containing region is enough to set the bit
    }
  }
  return WRITE_CELL;
}

int ExportProcessor::ProcessHeader(CellHeader& header) {
  m_header = header;
  for (const auto& t : header.GetDataTags()) {
    m_markers.push_back(t.id);
    m_kinds.push_back(t.type == Tag::MA_TAG ? 'M' : 'C');   // marker vs calculated
  }
  m_cols.resize(m_markers.size());
  return HEADER_NO_ACTION;   // no .byf output stream; we write our own pack
}

int ExportProcessor::ProcessLine(Cell& cell) {
  m_id.push_back(cell.id);
  m_x.push_back(cell.x);
  m_y.push_back(cell.y);
  m_cflag.push_back(cell.cflag);
  m_pflag.push_back(cell.pflag);
  const size_t n = std::min(m_cols.size(), cell.cols.size());
  for (size_t i = 0; i < n; ++i)            m_cols[i].push_back(cell.cols[i]);
  for (size_t i = n; i < m_cols.size(); ++i) m_cols[i].push_back(0.0f);
  return NO_WRITE_CELL;   // we don't re-emit cells
}

void ExportProcessor::finalize() {

  std::ofstream os(m_outpath, std::ios::binary);
  if (!os.good()) {
    std::cerr << "cyftools export: cannot open output: " << m_outpath << std::endl;
    return;
  }

  auto w    = [&](const void* p, std::size_t n) { os.write((const char*)p, (std::streamsize)n); };
  auto wu16 = [&](uint16_t v) { unsigned char b[2]={(unsigned char)v,(unsigned char)(v>>8)}; w(b,2); };
  auto wu32 = [&](uint32_t v) { unsigned char b[4]; for(int i=0;i<4;++i) b[i]=(unsigned char)(v>>(8*i)); w(b,4); };
  auto wu64 = [&](uint64_t v) { unsigned char b[8]; for(int i=0;i<8;++i) b[i]=(unsigned char)(v>>(8*i)); w(b,8); };

  const uint64_t n = m_id.size();

  // preamble: magic, version, n_cells, marker count + names
  w("CYFV", 4);
  wu16(5); wu16(0);          // version 5: trailing @FL flag maps (pflag then cflag)
  wu64(n);
  wu32(static_cast<uint32_t>(m_markers.size()));
  for (const auto& m : m_markers) { wu16(static_cast<uint16_t>(m.size())); w(m.data(), m.size()); }
  for (char k : m_kinds) os.put(k);   // 'M' = marker (MA), 'C' = calculated (CA)

  // pad to an 8-byte boundary so the column arrays below are aligned: a JS reader
  // can then make zero-copy BigUint64Array/Float32Array views straight over them.
  while ((static_cast<long long>(os.tellp()) & 7) != 0) os.put('\0');

  // column-major arrays (little-endian). u64 columns written element-wise to
  // guarantee byte order; f32 columns are contiguous IEEE-754 little-endian.
  for (uint64_t v : m_id)    wu64(v);
  w(m_x.data(), m_x.size() * sizeof(float));
  w(m_y.data(), m_y.size() * sizeof(float));
  for (uint64_t v : m_cflag) wu64(v);
  for (uint64_t v : m_pflag) wu64(v);
  for (const auto& col : m_cols)
    w(col.data(), col.size() * sizeof(float));

  // --- @RO polygons, so the pack is a self-contained viewer file ---
  const std::vector<Tag> tags = m_header.GetAllTags();
  uint32_t n_roi = 0;
  for (const auto& t : tags) if (t.type == Tag::RO_TAG) ++n_roi;
  wu32(n_roi);
  for (const auto& t : tags) {
    if (t.type != Tag::RO_TAG) continue;
    const std::string id = t.id;
    const std::string nm = t.GetField("NM");
    std::vector<float> pts;                       // x,y interleaved
    std::istringstream ps(t.GetField("PT"));
    std::string pair;
    while (ps >> pair) {
      const size_t c = pair.find(',');
      if (c == std::string::npos) continue;
      try { pts.push_back(std::stof(pair.substr(0, c)));
            pts.push_back(std::stof(pair.substr(c + 1))); }
      catch (...) { /* skip bad coord */ }
    }
    wu16(static_cast<uint16_t>(id.size())); w(id.data(), id.size());
    wu16(static_cast<uint16_t>(nm.size())); w(nm.data(), nm.size());
    wu32(static_cast<uint32_t>(pts.size() / 2));
    w(pts.data(), pts.size() * sizeof(float));
  }

  // --- @FL flag-bit maps (name -> bit), pflag block then cflag block ---
  // pflag lets a viewer derive a "phenotype gate" (lowest marker value among cells
  // whose bit is set); cflag carries the structural-flag names for labeling.
  auto writeFlags = [&](const char* reg){
    std::vector<std::pair<std::string,int>> fl;
    for (const auto& t : tags){
      if (t.type != Tag::FL_TAG || t.GetField("RG") != reg) continue;
      try { fl.emplace_back(t.id, std::stoi(t.GetField("BI"))); } catch (...) {}
    }
    wu32(static_cast<uint32_t>(fl.size()));
    for (const auto& f : fl){
      wu16(static_cast<uint16_t>(f.first.size())); w(f.first.data(), f.first.size());
      os.put(static_cast<char>(static_cast<unsigned char>(f.second & 0xff)));
    }
  };
  writeFlags("pflag");
  writeFlags("cflag");
}

// ---------------------------------------------------------------------------
// cohort: per-sample region densities via the "painting" area estimate
// ---------------------------------------------------------------------------
namespace {

// A rasterized union-of-discs accumulator. Each cell paints a filled disc of a
// fixed radius onto a grid of `pixel`-micron squares; the painted-pixel count ×
// pixel^2 is the union area, so two overlapping discs are counted once. Rows are
// 64-bit-word aligned (stride) so disjoint row-bands can be painted by separate
// threads with no shared word — i.e. the OpenMP path below is race-free.
struct PaintGrid {
  double x0 = 0, y0 = 0, pixel = 1.0;
  long   W = 0, H = 0, stride = 0;        // stride = uint64 words per row
  std::vector<uint64_t> bits;

  PaintGrid(double xmin, double ymin, double xmax, double ymax,
            double radius, double pixel_) : pixel(pixel_) {
    x0 = xmin - radius;                   // grow the box so edge discs fit
    y0 = ymin - radius;
    W = static_cast<long>(std::floor((xmax + radius - x0) / pixel)) + 1;
    H = static_cast<long>(std::floor((ymax + radius - y0) / pixel)) + 1;
    if (W < 1) W = 1;
    if (H < 1) H = 1;
    stride = (W + 63) / 64;
    bits.assign(static_cast<size_t>(stride) * static_cast<size_t>(H), 0ull);
  }

  inline void set(long ix, long iy) {     // ix,iy assumed in-bounds
    bits[static_cast<size_t>(iy) * static_cast<size_t>(stride) +
         static_cast<size_t>(ix >> 6)] |= (1ull << (ix & 63));
  }

  size_t painted() const {
    size_t c = 0;
    for (uint64_t w : bits) c += static_cast<size_t>(__builtin_popcountll(w));
    return c;
  }
};

// Paint discs of `radius` microns for the cells in `idx` (indices into xs/ys)
// into `g`. Threads own disjoint horizontal bands of rows, so writes never
// collide on a 64-bit word (rows are word-aligned in PaintGrid).
void paint_discs(PaintGrid& g, const std::vector<size_t>& idx,
                 const std::vector<float>& xs, const std::vector<float>& ys,
                 double radius, size_t nthreads) {
  const double Rpix  = radius / g.pixel;
  const long   Rceil = static_cast<long>(std::ceil(Rpix));
  const double R2    = Rpix * Rpix;

#ifdef HAVE_OMP
#pragma omp parallel num_threads(nthreads)
#endif
  {
    long y_lo = 0, y_hi = g.H;
#ifdef HAVE_OMP
    const int  tid  = omp_get_thread_num();
    const int  nth  = omp_get_num_threads();
    const long band = (g.H + nth - 1) / nth;
    y_lo = static_cast<long>(tid) * band;
    y_hi = std::min<long>(g.H, y_lo + band);
#endif
    for (size_t k = 0; k < idx.size(); ++k) {
      const size_t i  = idx[k];
      const long   cx = std::lround((xs[i] - g.x0) / g.pixel);
      const long   cy = std::lround((ys[i] - g.y0) / g.pixel);
      // only the rows this thread owns; rows outside [0,H) are skipped too
      long dlo = -Rceil, dhi = Rceil;
      if (cy + dlo < y_lo)     dlo = y_lo - cy;
      if (cy + dhi > y_hi - 1) dhi = y_hi - 1 - cy;
      for (long dy = dlo; dy <= dhi; ++dy) {
        const long iy = cy + dy;
        if (iy < 0 || iy >= g.H) continue;
        const double rem = R2 - static_cast<double>(dy) * static_cast<double>(dy);
        if (rem < 0) continue;
        const long hw   = static_cast<long>(std::floor(std::sqrt(rem)));
        long ix_lo = cx - hw; if (ix_lo < 0)      ix_lo = 0;
        long ix_hi = cx + hw; if (ix_hi >= g.W)   ix_hi = g.W - 1;
        for (long ix = ix_lo; ix <= ix_hi; ++ix)
          g.set(ix, iy);
      }
    }
  }
  (void)nthreads;
}

} // anonymous namespace

int CohortProcessor::ProcessHeader(CellHeader& header) {
  m_header = header;

  // Resolve this file's coordinate scale. Cohort paints in microns (the radius and
  // grid are micron-valued), so cells stored in pixels must be converted with the
  // file's own @HD MP (microns/pixel). This is per-file by design: each table
  // carries its own calibration; the cohort never assumes a shared pixel size.
  {
    const std::string un = header.GetHeaderField("UN");
    const std::string mp = header.GetHeaderField("MP");
    if (un == "micron") {
      m_scale = 1.0; m_scale_ok = true;
      m_scale_note = "coords already in microns (@HD UN:micron)";
    } else if (un == "pixel") {
      double mpv = 0.0;
      try { mpv = std::stod(mp); } catch (...) { mpv = 0.0; }
      if (mpv > 0.0) {
        m_scale = mpv; m_scale_ok = true;
        m_scale_note = "converted pixels -> microns via @HD MP:" + mp;
      } else {
        m_scale = 1.0; m_scale_ok = false;
        m_scale_note = "@HD UN:pixel but MP missing/invalid - cannot convert to microns";
      }
    } else if (un.empty()) {
      // legacy file with no declared units: assume microns (prior behavior)
      m_scale = 1.0; m_scale_ok = true;
      m_scale_note = "no @HD UN tag - assuming coords are already in microns";
    } else {
      m_scale = 1.0; m_scale_ok = false;
      m_scale_note = "unrecognized @HD UN:" + un + " - expected micron or pixel";
    }
  }

  // @FL bit -> name maps (the self-describing flag vocabulary, when present)
  for (const auto& t : header.GetFlagTags()) {
    const std::string reg = t.GetField("RG");
    int bit;
    try { bit = std::stoi(t.GetField("BI")); } catch (...) { continue; }
    if      (reg == "cflag") m_cflag_names[bit] = t.id;
    else if (reg == "pflag") m_pflag_names[bit] = t.id;
  }

  // marker-order pflag names (1st @MA marker -> bit 0, ...), as a fallback for
  // files that carry no @FL pflag tags. Matches how cyfview maps pflag bits.
  int midx = 0;
  for (const auto& t : header.GetDataTags())
    if (t.type == Tag::MA_TAG) m_pflag_auto[midx++] = t.id;

  // sample label from @SA, if any (else derived from the filename at emit time);
  // the GP field carries an optional group/category (set via `addtag --group`).
  const auto sa = header.GetSampleTags();
  if (!sa.empty()) { m_sample_name = sa.front().id; m_sample_group = sa.front().GetField("GP"); }

  return HEADER_NO_ACTION;   // cohort reads only; no .byf output stream
}

int CohortProcessor::ProcessLine(Cell& cell) {
  if (m_x.empty())
    m_sample_id = static_cast<uint64_t>(cell.id) >> 32;   // (sample_id<<32)|cell_id
  // store coordinates in microns (m_scale == 1 when already micron-valued)
  m_x.push_back(static_cast<float>(cell.x * m_scale));
  m_y.push_back(static_cast<float>(cell.y * m_scale));
  m_cflag.push_back(static_cast<uint64_t>(cell.cflag));
  m_pflag.push_back(static_cast<uint64_t>(cell.pflag));
  return NO_WRITE_CELL;       // cohort never re-emits cells
}

void CohortProcessor::WriteSampleJSON(std::ostream& os, const std::string& file,
                                      const std::string& ind) {
  const size_t N = m_x.size();

  // minimal JSON string escaper
  auto esc = [](const std::string& s) {
    std::string o; o.reserve(s.size() + 2);
    for (char c : s) {
      switch (c) {
        case '"':  o += "\\\""; break;
        case '\\': o += "\\\\"; break;
        case '\n': o += "\\n";  break;
        case '\t': o += "\\t";  break;
        case '\r': o += "\\r";  break;
        default:
          if (static_cast<unsigned char>(c) < 0x20) {
            char b[8]; std::snprintf(b, sizeof b, "\\u%04x", static_cast<unsigned char>(c)); o += b;
          } else o += c;
      }
    }
    return o;
  };

  // sample label: @SA, else the file's basename stem
  std::string sample = m_sample_name;
  if (sample.empty()) {
    std::string base = file;
    const size_t slash = base.find_last_of("/\\");
    if (slash != std::string::npos) base = base.substr(slash + 1);
    const size_t dot = base.find_last_of('.');
    if (dot != std::string::npos && dot != 0) base = base.substr(0, dot);
    sample = base;
  }

  // which cflag / pflag bits ever appear (OR over all cells), for the legends
  uint64_t cflag_any = 0, pflag_any = 0;
  for (size_t i = 0; i < N; ++i) { cflag_any |= m_cflag[i]; pflag_any |= m_pflag[i]; }

  std::vector<int> pbits;
  for (int b = 0; b < 64; ++b) if (pflag_any & (1ull << b)) pbits.push_back(b);

  // compartments -> synthetic "All" (bit -1) then each present cflag bit
  std::vector<int> cbits;
  cbits.push_back(-1);
  for (int b = 0; b < 64; ++b) if (cflag_any & (1ull << b)) cbits.push_back(b);

  // sample-level bbox
  double gxmin = 0, gymin = 0, gxmax = 0, gymax = 0;
  if (N) {
    gxmin = gxmax = m_x[0]; gymin = gymax = m_y[0];
    for (size_t i = 1; i < N; ++i) {
      gxmin = std::min<double>(gxmin, m_x[i]); gxmax = std::max<double>(gxmax, m_x[i]);
      gymin = std::min<double>(gymin, m_y[i]); gymax = std::max<double>(gymax, m_y[i]);
    }
  }

  // name resolvers: header @FL first, then built-in cflag vocabulary, then generic
  auto cname = [&](int bit) -> std::string {
    if (bit < 0) return "All";
    auto it = m_cflag_names.find(bit);
    if (it != m_cflag_names.end()) return it->second;
    for (const auto& f : cyf::standardFlags())
      if (std::string(f.reg) == "cflag" && f.bit == bit) return f.name;
    return "cflag_bit_" + std::to_string(bit);
  };
  auto pname = [&](int bit) -> std::string {
    auto it = m_pflag_names.find(bit);
    if (it != m_pflag_names.end()) return it->second;       // @FL declaration wins
    auto ja = m_pflag_auto.find(bit);
    if (ja != m_pflag_auto.end()) return ja->second;        // else @MA marker order
    return "pflag_bit_" + std::to_string(bit);
  };

  std::cerr << "cyftools cohort:   " << sample << ": " << AddCommas(N)
            << " cells, " << cbits.size() << " compartments, "
            << pbits.size() << " phenotype bits" << std::endl;

  // --- precompute each compartment's painted area (the density denominators) ---
  // A compartment is one cflag bit's region, or "All" (every cell, bit -1).
  struct Comp { int bit; size_t n; double area_mm2; };
  std::vector<Comp> comps;
  comps.reserve(cbits.size());
  for (int b : cbits) {
    std::vector<size_t> idx;
    if (b < 0) { idx.resize(N); for (size_t i = 0; i < N; ++i) idx[i] = i; }
    else { for (size_t i = 0; i < N; ++i) if (m_cflag[i] & (1ull << b)) idx.push_back(i); }

    double xmin = 0, ymin = 0, xmax = 0, ymax = 0;
    for (size_t k = 0; k < idx.size(); ++k) {
      const size_t i = idx[k];
      if (k == 0) { xmin = xmax = m_x[i]; ymin = ymax = m_y[i]; }
      else {
        xmin = std::min<double>(xmin, m_x[i]); xmax = std::max<double>(xmax, m_x[i]);
        ymin = std::min<double>(ymin, m_y[i]); ymax = std::max<double>(ymax, m_y[i]);
      }
    }

    std::cerr << "cyftools cohort:     compartment '" << cname(b) << "' ("
              << AddCommas(idx.size()) << " cells) painting..." << std::endl;

    double area_mm2 = 0.0;
    if (!idx.empty()) {
      PaintGrid g(xmin, ymin, xmax, ymax, m_radius, m_pixel);
      paint_discs(g, idx, m_x, m_y, m_radius, m_nthreads);
      area_mm2 = static_cast<double>(g.painted()) * m_pixel * m_pixel / 1.0e6;
    }
    comps.push_back({b, idx.size(), area_mm2});
  }

  // --- joint (cflag, pflag) histogram: the contingency table of the two flag ---
  // registers. The viewer counts any pflag combination in any compartment as a
  // masked sum over these rows, divided by that compartment's painted area:
  //   count = sum(row.count) over rows where (bit==-1 || row.cflag & (1<<bit))
  //                                       and  pflagPredicate(row.pflag)
  std::map<uint64_t, std::map<uint64_t, uint64_t>> joint;
  for (size_t i = 0; i < N; ++i) joint[m_cflag[i]][m_pflag[i]]++;
  size_t njoint = 0;
  for (const auto& cf : joint) njoint += cf.second.size();
  std::cerr << "cyftools cohort:     joint flag histogram: " << AddCommas(njoint)
            << " distinct (cflag,pflag) patterns" << std::endl;

  // --- emit the per-sample object ---
  os << ind << "{\n";
  os << ind << "  \"file\": \""   << esc(file)   << "\",\n";
  os << ind << "  \"sample\": \"" << esc(sample) << "\",\n";
  os << ind << "  \"group\": \""  << esc(m_sample_group) << "\",\n";   // @SA GP field, "" if untagged
  os << ind << "  \"sample_id\": " << m_sample_id << ",\n";
  os << ind << "  \"n_cells\": "   << N << ",\n";
  os << ind << "  \"bbox_um\": [" << gxmin << ", " << gymin << ", "
                                  << gxmax << ", " << gymax << "],\n";

  // legends: bit -> name for each register (drives the viewer's predicate UI)
  os << ind << "  \"cflag_bits\": [";
  bool cfirst = true;
  for (int b : cbits) {
    if (b < 0) continue;                       // skip the synthetic "All"
    os << (cfirst ? " " : ", "); cfirst = false;
    os << "{\"name\": \"" << esc(cname(b)) << "\", \"bit\": " << b << "}";
  }
  os << " ],\n";
  os << ind << "  \"pflag_bits\": [";
  for (size_t j = 0; j < pbits.size(); ++j)
    os << (j ? ", " : " ")
       << "{\"name\": \"" << esc(pname(pbits[j])) << "\", \"bit\": " << pbits[j] << "}";
  os << " ],\n";

  // compartments: precomputed painted areas (the density denominators)
  os << ind << "  \"compartments\": [\n";
  for (size_t k = 0; k < comps.size(); ++k) {
    const Comp&  c    = comps[k];
    const double dens = c.area_mm2 > 0 ? static_cast<double>(c.n) / c.area_mm2 : 0.0;
    os << ind << "    {\"cflag\": \"" << esc(cname(c.bit)) << "\", \"bit\": " << c.bit
       << ", \"n_cells\": " << c.n << ", \"area_mm2\": " << c.area_mm2
       << ", \"density_per_mm2\": " << dens << "}"
       << (k + 1 < comps.size() ? "," : "") << "\n";
  }
  os << ind << "  ],\n";

  // joint histogram: [cflag_value, pflag_value, count]. The flag *values* are
  // decimal strings, because a 64-bit register can exceed JS's 2^53 safe-integer
  // range -> a consumer parses them with BigInt and masks with BigInt.
  os << ind << "  \"flag_histogram\": [";
  bool first = true;
  for (const auto& cf : joint) {
    for (const auto& pf : cf.second) {
      os << (first ? "\n" : ",\n");
      first = false;
      os << ind << "    [\"" << cf.first << "\", \"" << pf.first << "\", " << pf.second << "]";
    }
  }
  os << (first ? "" : "\n") << ind << "  ]\n";
  os << ind << "}";
}


int ROIProcessor::ProcessLine(Cell& cell) {

  // Loop the table and check if the cell is in the ROI
  bool print_line = true; //m_blacklist_remove;

  // Loop through all polygons and check if the point is inside any of them
  for (const auto &polygon : m_rois) {
    
    int number = 0; // region number
    
    // if point is in this polygon, add the polygon id number to the roi
    if (polygon.PointIn(cell.x,cell.y)) {
      
      if (roikey(polygon, "lacklist") || roikey(polygon, "rtifact") || roikey(polygon, "error") || roikey(polygon, "artefact")) {
	SET_FLAG(cell.cflag, ARTIFACT_FLAG);		
	print_line = false; 
      } else if (roikey(polygon, "ormal")) {
	CLEAR_FLAG(cell.cflag, TUMOR_FLAG);
	CLEAR_FLAG(cell.cflag, MARGIN_FLAG);	
	CLEAR_FLAG(cell.cflag, TUMOR_MANUAL_FLAG);
	SET_FLAG(cell.cflag, NORMAL_FLAG);	
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
  // pflag bit = the marker's position among @MA tags (1st marker -> bit 0); the
  // column index tracks its position among all data columns (for reading values).
  size_t i = 0, mbit = 0;
  for (const auto& t : header.GetDataTags()) {
    if (t.type == Tag::MA_TAG)
      m_marker_map[t.id] = { i, mbit++ };
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

    const size_t col = m->second.first;    // column index (read the value)
    const size_t bit = m->second.second;   // pflag bit = position among @MA markers

    // set the flag on if it clears the gates
    if (cell.cols.at(col) >= b.second.first*m_scale &&
	cell.cols.at(col) <= b.second.second*m_scale) {
      flag.setFlagOn(bit);
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
  } else if (m_list_markers) {
    m_header.PrintMarkers();
    m_header_only = true; // to signify nothing more to print for this
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
    
    cell.PrintWithHeader(m_round, m_tabprint, m_adjacent, m_header, false, m_csv_header);

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

  cell.PrintWithHeader(m_round, m_tabprint, m_adjacent, m_header, m_strict_cut, m_csv_header);

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

  // declare the standard cflag bit meanings so the file is self-describing
  cyf::addStandardFlagTags(m_header);

  // ensure an @HD format/version line leads the header
  m_header.EnsureHeaderLine();

  // stamp the required coordinate-scale tags onto @HD (set by `convert`)
  if (!m_mpp.empty())   m_header.SetHeaderField("MP", m_mpp);
  if (!m_units.empty()) m_header.SetHeaderField("UN", m_units);

  assert(!m_filename.empty());

  // format follows the output extension (.cyf text, .byf binary); stdout "-" -> binary
  const cyf::OutFormat fmt = cyf::formatForPath(m_filename);
  if (m_filename == "-") {
    m_archive = std::make_unique<OutArchive>(std::cout, fmt);
  } else {
    m_os = std::make_unique<std::ofstream>(m_filename, std::ios::binary);
    if (!m_os->is_open())
      throw std::runtime_error("cannot open output file '" + m_filename + "' for writing");
    m_archive = std::make_unique<OutArchive>(*m_os, fmt);
  }

  m_header.SortTags();
  
  // archive the header
  (*m_archive)(m_header);
  
  return 0; // minor, but the HEADER_NO_ACTION is in CellProcessor abstract class, of which CerealProcessor is uniquely not a child
}

int CerealProcessor::ProcessLine(const std::string& line) {

  // Legacy positional parse when no column plan was supplied.
  if (m_col_kind.empty()) {
    Cell row(line, m_id_index, m_x_index, m_y_index, m_header, m_cellid, m_sampleid);
    if (m_id_index < 0) row.set_cell_id(m_cellid);
    m_cellid++;
    (*m_archive)(row);
    return 0;
  }

  const StringVec tokens = tokenize_comma_delimited<StringVec>(line);
  if (tokens.size() != m_col_kind.size())
    throw std::runtime_error("cyftools convert: line has " + std::to_string(tokens.size()) +
                             " columns but header has " + std::to_string(m_col_kind.size()) +
                             "; line: " + line);

  Cell cell;
  cell.id = 0; cell.x = 0; cell.y = 0; cell.cflag = 0; cell.pflag = 0;
  uint32_t cellid = m_cellid;

  // markers and metas collected separately then concatenated, so cols line up
  // with GetDataTags() after SortTags (all @MA in column order, then all @CA).
  std::vector<float> markers, metas;

  for (size_t i = 0; i < tokens.size(); ++i) {
    const char* t = tokens[i].c_str();
    switch (m_col_kind[i]) {
      case COL_ID:     cellid = static_cast<uint32_t>(std::strtol(t, nullptr, 10)); break;
      case COL_X:      cell.x = std::strtof(t, nullptr); break;
      case COL_Y:      cell.y = std::strtof(t, nullptr); break;
      case COL_MARKER: markers.push_back(std::strtof(t, nullptr)); break;
      case COL_CA:     metas.push_back(std::strtof(t, nullptr)); break;
      case COL_IC: {
        const float v = std::strtof(t, nullptr);
        metas.push_back(v);
        if (v != 0.0f) cell.cflag |= (static_cast<cy_uint>(1) << 5);   // bit 5 (IC)
        break;
      }
      case COL_GATE: {
        const int b = m_gate_bit[i];
        if (std::strtof(t, nullptr) != 0.0f && b >= 0 &&
            b < static_cast<int>(sizeof(cy_uint) * 8))
          cell.pflag |= (static_cast<cy_uint>(1) << b);                // gate -> pflag bit
        break;
      }
      case COL_REGION:
        if (static_cast<int>(std::strtof(t, nullptr)) == 1)
          cell.cflag |= (static_cast<cy_uint>(1) << 3);                // bit 3 (Region)
        break;
      default: break;   // COL_IGNORE
    }
  }

  cell.cols = std::move(markers);
  cell.cols.insert(cell.cols.end(), metas.begin(), metas.end());

  cell.set_sample_id(m_sampleid);
  cell.set_cell_id(cellid);
  m_cellid++;

  (*m_archive)(cell);
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
