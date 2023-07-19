#include "cell_processor.h"
#include "polygon.h"
#include "cell_utils.h"
#include "cell_flag.h"

#include "cell_row.h"

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

  for (size_t i = 0; i < cell.m_cols.size(); i++) {
    sums[i] += cell.m_cols.at(i);
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

  cell.m_cols = means;

  // write the one cell
  OutputLine(cell);
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

  if (m_dis(m_gen) <= m_rate) {
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
  
  // get the flag value from the line
  CellFlag pflag(cell.m_pheno_flag);
  CellFlag cflag(cell.m_cell_flag);  
  
  // test it and print line if so
  bool pflags_met = pflag.testAndOr(m_por, m_pand);
  bool cflags_met = cflag.testAndOr(m_cor, m_cand);  

  //std::cerr << "m_por " << m_por << " m_pand " << m_pand << " pflag " << pflag << " test " << pflags_met << " m_pnot " << m_pnot << std::endl;
  //std::cerr << "m_cor " << m_cor << " m_cand " << m_cand << " cflag " << cflag << " test " << cflags_met << " m_cnot " << m_cnot << std::endl;  
  
  // if flags met, print the cell
  if ( (pflags_met != m_pnot) && (cflags_met != m_cnot)) {
    write_cell = true;
    //std::cerr << " writing cell " << std::endl;
  }

  ///////
  // FIELD
  ///////
  bool flag_write_cell = write_cell;
  write_cell = m_or_toggle ? false : write_cell; // if or toggle is on, then ignore flag criteria and set true
  
  for (const auto& c : m_criteria_int) {
    
    assert(c.first < cell.m_cols.size());
    float value = cell.m_cols.at(c.first);
    
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

/*int RadialProcessor::ProcessHeader(CellHeader& header) {

  m_header = header;

  m_header.addTag(Tag(Tag::PG_TAG, "", m_cmd));
  
  // build a set (just for this method) to compare existing tags with
  std::unordered_set<std::string> tag_set; 
  for (const auto& t : m_header.GetDataTags())
    tag_set.insert(t.id);
  
  // assure tags aren't already there
  for (const auto& l : m_label) {
    
    // warn if already there
    if (tag_set.count(l)) {
      std::cerr << "Warning: header already contains column: " <<
	l << std::endl;
    } else {
      // otherwise add the new tags
      Tag dtag(Tag::CA_TAG, l, "");
      m_header.addTag(dtag);
    }
  }

  m_header.SortTags();
  
  // just in time output, so as not to write an empty file if the input crashes
  // set the output to file or stdout
  this->SetupOutputStream(); 

  // output the header
  assert(m_archive);
  (*m_archive)(m_header);
  
  return HEADER_NO_ACTION;
  
}

int RadialProcessor::ProcessLine(Cell& cell) {

  std::vector<float> cell_count(m_inner.size());

  assert(cell.m_spatial_ids.size() == cell.m_spatial_flags.size());
  assert(cell.m_spatial_ids.size() == cell.m_spatial_dist.size());
  
  // loop the nodes connected to each cell
  for (size_t i = 0; i < cell.m_spatial_ids.size(); i++) { 

    // test if the connected cell meets the flag criteria
    // n.first is cell_id of connected cell to this cell
    for (size_t j = 0; j < m_inner.size(); j++) {
      
      CellFlag tflag(cell.m_spatial_flags.at(i));
      
      // both are 0, so take all cells OR it meets flag criteria
      if ( (!m_logor[j] && !m_logand[j]) ||
	   tflag.testAndOr(m_logor[j], m_logand[j])) {
	
	// then increment cell count if cell in bounds
	cell_count[j] += cell.m_spatial_dist.at(i) >= m_inner[j] &&
    	              cell.m_spatial_dist.at(i) <= m_outer[j];

      }
    }
  }
  
  // calculate the density
  std::vector<float> area(m_inner.size());
  for (size_t j = 0; j < m_inner.size(); ++j) {
    float outerArea = static_cast<float>(m_outer[j]) * static_cast<float>(m_outer[j]) * 3.1415926535f;
    float innerArea = static_cast<float>(m_inner[j]) * static_cast<float>(m_inner[j]) * 3.1415926535f;
    area[j] = outerArea - innerArea;
  }
  
  // do the density calculation for each condition
  // remember, i is iterator over cells, j is over conditions
  for (size_t j = 0; j < area.size(); ++j) {
    float value = cell.m_spatial_ids.empty() ? 0 : cell_count[j] * 1000000 / area[j]; // density per 1000 square pixels
    cell.m_cols.push_back(value);
  }

  return WRITE_CELL;
}
*/

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
  std::vector<float> m_cols_new;
  assert(cell.m_cols.size() >= m_to_remove.size());
  m_cols_new.reserve(cell.m_cols.size() - m_to_remove.size());
  
  for (size_t i = 0; i < cell.m_cols.size(); i++) {
    if (!m_to_remove.count(i)) {
      m_cols_new.push_back(cell.m_cols.at(i));
    }
  }
  cell.m_cols = m_cols_new;
  
  return WRITE_CELL;
}

int HeadProcessor::ProcessLine(Cell& cell) {

  m_current_n++;
  if (m_current_n <= m_n) {
    return WRITE_CELL;
  }
  return NO_WRITE_CELL;
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
  for (size_t i = 0; i < cell.m_cols.size(); i++) {
    if (m_to_log.count(i)) {
      if (cell.m_cols[i] > 0) {
	cell.m_cols[i] = std::log10(cell.m_cols.at(i));
      } else {
	
	if (!m_bool_warning_emitted) {
	  std::vector<Tag> tags = m_header.GetDataTags();
	  m_bool_warning_emitted = true;
	  std::cerr << "Warning: encountered zero or negative number to log on " <<
	    "line " << m_count << " with value " << cell.m_cols[i] << " on column " <<
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
  std::vector<float> m_cols_new;
  assert(cell.m_cols.size() >= m_to_remove.size());
  m_cols_new.reserve(cell.m_cols.size() - m_to_remove.size());
  
  for (size_t i = 0; i < cell.m_cols.size(); i++) {
    if (m_to_remove.count(i) == 0) {
      m_cols_new.push_back(cell.m_cols.at(i));
    }
  }
  cell.m_cols = m_cols_new;

  // clean the graph
  if (m_clean_graph) {
    cell.m_spatial_ids.clear();
    cell.m_spatial_dist.clear();
    cell.m_spatial_flags.clear();
  }

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

/*int TumorProcessor::ProcessHeader(CellHeader& header) {
  
  m_header = header;
  m_header.addTag(Tag(Tag::PG_TAG, "", m_cmd));

  // just in time output, so as not to write an empty file if the input crashes
  // set the output to file or stdout
  this->SetupOutputStream(); 
  
  m_header.SortTags();

  // output the header
  assert(m_archive);
  (*m_archive)(m_header);

  return HEADER_NO_ACTION;
}

int TumorProcessor::ProcessLine(Cell& cell) {

  // make a copy of the cell for output
  Cell cell_new = cell;

  assert(cell.m_spatial_ids.size() == cell.m_spatial_flags.size());
  assert(cell.m_spatial_ids.size() == cell.m_spatial_dist.size());

  }*/

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
    if (cell.m_cols.at(m->second) >= b.second.first &&
	cell.m_cols.at(m->second) <= b.second.second) {
      //      std::cerr << "Marker: " << b.first << " cleared" << std::endl;
      //std::cerr << "...before " << flag << std::endl; 
      flag.setFlagOn(marker_index);
    }
    
  }

  // convert to cy_uint for storage
  cell.m_pheno_flag = flag.toBase10();
  
  return 1;
}

int ViewProcessor::ProcessHeader(CellHeader& header) {

  // print it, that's all
  if (m_print_header || m_header_only) 
    header.Print();
  
  return m_header_only ? ONLY_WRITE_HEADER : HEADER_NO_ACTION;
  
}

int ViewProcessor::ProcessLine(Cell& cell) {

  cell.Print(m_round);
    
  return NO_WRITE_CELL; // don't output, since already printing it
}

int CatProcessor::ProcessLine(Cell& line) {

  // update the cell id
  size_t n_cell_id = m_offset + line.m_id; 
  if (n_cell_id > m_max_cellid)
    m_max_cellid = n_cell_id;
  line.m_id = n_cell_id;
  
  if (line.m_spatial_ids.size() && !m_error_emitted) {
    std::cerr << "Warning: Graph concatenation not supported" << std::endl;
    m_error_emitted = true;
  }
  
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

  Cell row(line, m_header);

  // serialize it
  (*m_archive)(row);

  // serializing here, so don't need to write in the cell_table call
  // minor, but NO_WRITE_CELL is in CellProcessor abstract class, of which CerealProcessor is uniquely not a child  
  return 0; 
  
}

int BuildProcessor::ProcessHeader(CellHeader& header) {
  m_header = header; // store but don't print

  return 2;
}

// just a pass through to place it in the table
int BuildProcessor::ProcessLine(Cell& cell) {

  return 2;
  
}
