#include "cell_processor.h"
#include "cell_utils.h"
#include "cell_flag.h"

int SelectProcessor::ProcessHeader(CellHeader& header) {

  // find which one the cellflag is
  size_t ii = 0;
  for (const auto& t : header.GetTags()) {
    if (t.isFlagTag()) {
      m_flag_index = ii;
    }
    ii++;
  }
  
  // must have a flag
  if (m_flag_index == static_cast<size_t>(-1)) {
    throw std::runtime_error("Header does not contain and flag tags");
  }
  
  // print the completed header
  if (m_print_header)
    header.Print();

  return 0;
}

int SelectProcessor::ProcessLine(const std::string& line) {

  // get the flag value from the line
  uint64_t flag_val = get_nth_element_as_integer(line, m_flag_index);
  CellFlag flag(flag_val);

  // test it and print line if so
  //if (flag.test(m_on, m_off))
  bool flags_met = flag.testAndOr(m_or, m_and);
  if (flags_met != m_not)
    std::cout << line << std::endl;

  return 0;
}

int CutProcessor::ProcessHeader(CellHeader& header) {

  // if not a strict cut, keep dims and id
  if (!m_strict_cut) {
    m_include.insert(header.GetX());
    m_include.insert(header.GetY());
    m_include.insert(header.GetZ());
    m_include.insert(header.GetID());
  }
  m_include.erase(std::string());

  // find indicies of columns to remove
  const std::vector<Tag>& tags = header.GetTags();
  size_t i = 0;
  for (const auto& t : tags) {

    // skip version tag
    if (t.isVersionTag())
      continue;
    
    // if not a strict cut, add dims and cellid
    if (!m_strict_cut &&
	(t.isDimTag() || t.isIDTag())) {
      i++;
      continue;
    } 
    
    if (!m_include.count(t.GetName())) {
      to_remove.insert(i);
    }

    i++;
    
  }

  // cut down the header
  header.Cut(m_include);

  if (m_header_print) {
    header.Print();
  }

  return 0;
}

int CutProcessor::ProcessLine(const std::string& line) {

  // tokenize the input line
  std::vector<std::string> tokens = tokenize_comma_delimited(line);

  // have a place to store the output tokens
  std::vector<std::string> output;
  output.reserve(tokens.size());
  
  // fill only tokens to output not labeled for removal
  for (size_t i = 0; i < tokens.size(); ++i) {
    if (!to_remove.count(i)) {
      output.push_back(tokens.at(i));
    }
  }
  
  // concatenate output tokens
  std::cout << tokens_to_comma_string(output) << std::endl;

  return 0;
}

int LogProcessor::ProcessHeader(CellHeader& header) {

  // setup which are marker indicies
  size_t i = 0;
  for (const auto& t : header.GetTags()) {

    if (t.isVersionTag())
      continue;
    
    if (t.isMarkerTag())
      m_to_log.insert(i);
    i++;
  }
  
  // print the header
  if (m_print_header)
    header.Print();

  return 0;
}

int LogProcessor::ProcessLine(const std::string& line) {

  m_line_number++;
  
  // tokenize the line
  std::vector<std::string> tokens = tokenize_comma_delimited(line);

  // log10 it
  std::vector<std::string> output;
  output.reserve(tokens.size());

  size_t i = 0;
  for (const auto& t : tokens) {
    
    // don't log it
    if (!m_to_log.count(i)) {
      output.push_back(t);
      // log it
    } else {
      try {
	float value = std::stof(t);

	float log_value;
	
	if (value <= 0) {
	  //throw std::invalid_argument("Cannot take the logarithm of a zero or negative number: line" + std::to_string(m_line_number) + " at token " + std::to_string(i+1));
	  std::cerr << "Warning: Cannot log <= 0 (setting output to size_t(-1)): line " + std::to_string(m_line_number) + " at token " + std::to_string(i+1) << std::endl;
	  log_value = static_cast<size_t>(-1);
	} else {
	  log_value = std::log10(value);
	}
	
	std::ostringstream oss;
	oss << log_value;
	output.push_back(oss.str());
      } catch (const std::invalid_argument& e) {
	throw std::invalid_argument("The input token '" + t + "' is not a valid number. Line: " + std::to_string(m_line_number) + " token " + std::to_string(i+1));
      } catch (const std::out_of_range& e) {
	throw std::out_of_range("The input token '" + t + "' is out of range for conversion to float. Line: " + std::to_string(m_line_number) + " token " + std::to_string(i+1));
      }
      
    }
    i++;
  }

  std::cout << tokens_to_comma_string(output) << std::endl;

  return 0;
  
}


int ROIProcessor::ProcessHeader(CellHeader& header) {

  // how many columns per line
  m_col_count = header.ColumnCount();

  // which column is x and y
  x_i = header.whichColumn(header.GetX());
  y_i = header.whichColumn(header.GetY());  
  
  // save a copy of the original header, for reading lines below
  m_header = header;
  
  // add the roi tag
  Tag roi_tag("CA","roi");
  header.addTag(roi_tag);

  // print the header
  if (m_print_header)
    header.Print();

  return 0;
}

int ROIProcessor::ProcessLine(const std::string& line) {

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

  return 0;
  
}

int ViewProcessor::ProcessHeader(CellHeader& header) {

  // print it, that's all
  if (m_print_header || m_header_only) 
    header.Print();
  
  return m_header_only ? 1 : 0;
  
}

int ViewProcessor::ProcessLine(const std::string& line) {

  // should never get here if header-only
  assert(!m_header_only);

  // no rounding, so just dump line
  if (m_round < 0) {
    std::cout << line << std::endl;
    return 0;
  }

  // rounding, so have to read the line and round
  std::vector<std::string> tokens = tokenize_comma_delimited(line);

  // store the modified strings
  std::vector<std::string> output;
  output.reserve(tokens.size());

  // round the string
  for (const auto& t : tokens)
    output.push_back(round_string(t, m_round));

  // concatenate output tokens
  std::cout << tokens_to_comma_string(output) << std::endl;
  
  return 0;
}
