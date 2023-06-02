#include "cell_processor.h"
#include "cell_utils.h"
#include "cell_flag.h"

#include "cell_row.h"

int SelectProcessor::ProcessHeader(CellHeader& header) {

  m_header = header;

  m_header.addTag(Tag(Tag::PG_TAG, "", m_cmd));

  m_header.SortTags();
  
  // just in time, make the output stream
  this->SetupOutputStream();
  
  // output the header
  assert(m_archive);
  (*m_archive)(m_header);
  
  return 0;
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
  
  return 0;
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

int SelectProcessor::ProcessLine(Cell& cell) {

  // get the flag value from the line
  CellFlag flag(cell.m_flag);
  
  // test it and print line if so
  bool flags_met = flag.testAndOr(m_or, m_and);

  // if flags met, print the cell
  if (flags_met != m_not)
    return CellProcessor::WRITE_CELL;

  return CellProcessor::NO_WRITE_CELL; // don't write if not selected
}

int RadialProcessor::ProcessHeader(CellHeader& header) {

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
  
  return 0;
  
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

  return 1;
}


int CountProcessor::ProcessHeader(CellHeader& header) {
  return 0; // do nothing
}

int CountProcessor::ProcessLine(Cell& cell) {
  m_count++;
  return 0; // do nothing
}

void CountProcessor::PrintCount() {
  std::cout << m_count << std::endl;
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
  
  return 0;
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
  
  return 1;
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

  return 0;
}

int LogProcessor::ProcessLine(Cell& cell) {

  m_count++;
  for (size_t i = 0; i < cell.m_cols.size(); i++) {
    if (m_to_log.count(i)) {
      if (cell.m_cols[i] > 0) {
	cell.m_cols[i] = std::log10(cell.m_cols.at(i));
      } else {
	
	if (!m_bool_warning_emitted) {
	  m_bool_warning_emitted = true;
	  std::cerr << "Warning: encountered zero or negative number to log on " <<
	    "line " << m_count << " - removing this line in output. This warning " <<
	    "will emit only once per file";
	}
	
	return 0;
      }
    }
  }
    
  return 1;
  
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

  return 0;
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

  return 0;

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
  
  return 0;
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
  cell.m_flag = flag.toBase10();
  
  return 1;
}

int ViewProcessor::ProcessHeader(CellHeader& header) {

  // print it, that's all
  if (m_print_header || m_header_only) 
    header.Print();
  
  return m_header_only ? 1 : 0;
  
}

int ViewProcessor::ProcessLine(Cell& cell) {

  cell.Print(m_round);
    
  return 0; // don't output, since already printing it
}

int CatProcessor::ProcessLine(Cell& line) {
  //int CatProcessor::ProcessLine(const std::string& line) {

  //std::cerr << " processingl ine " << line << std::endl;
  //std::cerr << " graph indicies " << m_graph_indicies.size() << std::endl;
  /*  
  // get the cell id
  size_t o_cell_id = get_nth_element_as_integer(line, m_cellid_index);

  // update the new cell id and max
  size_t n_cell_id = m_offset + o_cell_id;
  if (n_cell_id > m_max_cellid)
    m_max_cellid = n_cell_id;

  
  //  std::cerr << "old cell id: " << o_cell_id << " new " << n_cell_id <<
  // " offset " << m_offset << " max " << this->GetMaxCellID() << std::endl;

  // get the tokens
  std::vector<std::string> tokens = tokenize_comma_delimited(line);

  // update the graph ids
  for (const auto& i : m_graph_indicies) {
    const std::string& graph_line = tokens.at(i);

    // parse the node and rest the cell-ids with the offset
    CellNode node(graph_line);
    node.OffsetNodes(m_offset);

    tokens[i] = node.toString(false); // false is for "integerize"
  }
  */
  
  /*
  // update the cell id
  tokens[m_cellid_index] = std::to_string(n_cell_id);

  // add the sample number
  tokens.emplace_back(std::to_string(m_sample));
  
  // output
  std::cout << tokens_to_comma_string(tokens) << std::endl;

  */
  return 1;
}

int CatProcessor::ProcessHeader(CellHeader& header) {
  
  /*
  // if we already have the master header, just compare for error checking
  if (m_master_set) {

    // loop both master header and new header
    // but allow for master to have extra sample column
    std::vector<std::string> header_cols = header.GetColOrder();
    std::vector<std::string> master_cols = m_master_header.GetColOrder();
    
    // remove "sample" from master_cols
    master_cols.erase(std::remove_if(master_cols.begin(), master_cols.end(),
				     [](const std::string& s) { return s == "sample"; }), master_cols.end());
    
    if (!std::equal(header_cols.begin(), header_cols.end(), master_cols.begin(), master_cols.end())) {
      throw std::runtime_error("Error: All headers have to have the same number of columns in same order");
    }

    return 0;
  }
  
  // find which one the cell id is
  size_t ii = 0;
  for (const auto& t : header.GetColTags()) {

    // get the ID tag
    if (t.isIDTag()) {

      // check that there is only one ID tag
      if (m_cellid_index != static_cast<size_t>(-1))
	throw std::runtime_error("Error: Cell ID already found - does your header have two?");
      m_cellid_index = ii;
    }

    // set graph tags
    if (t.isGraphTag())
      m_graph_indicies.push_back(ii);
    
    ii++;
  }
  
  // throw error if not found
  if (m_cellid_index == static_cast<size_t>(-1)) {
    throw std::runtime_error("Error: no cell id found. ID tag in header required");
  }

  // save a copy of the first header to be the master header
  if (!m_master_header.size()) {
    m_master_header = header;
    m_master_set = true;
  }

  // add the sample tag if not already there
  if (!m_master_header.hasTag("sample")) {
    Tag sample_tag("CA","sample");
    m_master_header.addTag(sample_tag);

  }

  // add the command tag
  Tag cmd_tag("CN", "cysift cat");
  m_master_header.addTag(cmd_tag);
  
  // print the header
  if (m_print_header && m_master_set)
    m_master_header.Print();
  */
  
  return 0;
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
  
  return 0;
}

int CerealProcessor::ProcessLine(const std::string& line) {

  Cell row(line, m_header);

  // serialize it
  (*m_archive)(row);
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
