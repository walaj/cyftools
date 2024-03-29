#include "cell_row.h"
#include <cassert>

#include "cell_utils.h"

static void printWithoutScientific(double number, int round) {

  // Convert to string with full precision
  std::string out = std::to_string(number);
  
  // Remove trailing zeros
  out.erase(out.find_last_not_of('0') + 1, std::string::npos);
  
  // Ensure only "round" characters after the decimal
  size_t decimalPos = out.find('.');
  if (decimalPos != std::string::npos && out.size() > decimalPos + round + 1) {
    out.erase(decimalPos + round + 1);
  }
  
  // Remove trailing decimal point if no fractional part
  if (out.back() == '.') {
    out.pop_back();
  }
  
  std::cout << out;
  
  /*
  if (std::abs(number - std::round(number)) < 1e-6) {
        // If the number is very close to an integer, print as an integer
        std::cout << static_cast<long long>(number);
    } else {
        // Else, print with full precision, without trailing zeros
        std::cout << std::fixed;
        std::string out = std::to_string(number);
        out.erase(out.find_last_not_of('0') + 1, std::string::npos);  // Remove trailing zeros
        if (out.back() == '.') {
            out.pop_back();  // Remove trailing decimal point if no fractional part
        }
        std::cout << out;
	}*/
}

void Cell::PrintForCrevasse(const CellHeader& header) const {

  char d = ',';

  uint32_t sampleID = static_cast<uint32_t>(id >> 32); // Extract the first 32 bits
  uint32_t cellID = static_cast<uint32_t>(id & 0xFFFFFFFF); // Extract the second 32 bits

  std::cout << cellID << d << x << d << y << d;

  size_t i = 0;
  const auto& tags = header.GetDataTags();
  for (auto it = tags.begin(); it != tags.end(); ++it) {
    if (it->type == Tag::MA_TAG) {
      std::cout << cols.at(i);
      if (std::next(it) != tags.end() && std::next(it)->type == Tag::MA_TAG) {
	std::cout << ","; 
      }
    }
    i++;
  }
  
  std::cout << std::endl;
}

void Cell::set_sample_id(uint32_t new_id) {
  uint32_t currentCellID = static_cast<uint32_t>(id & 0xFFFFFFFF);
  id = static_cast<uint64_t>(new_id) << 32 | currentCellID; 
}

void Cell::set_cell_id(uint32_t new_id) {
  uint32_t currentSampleID = static_cast<uint32_t>(id >> 32);
  id = static_cast<uint64_t>(currentSampleID) << 32 | new_id;
}

uint32_t Cell::get_sample_id() const {
    return static_cast<uint32_t>(id >> 32);
}

uint32_t Cell::get_cell_id() const {
    return static_cast<uint32_t>(id & 0xFFFFFFFF);
}

std::ostream& operator<<(std::ostream& os, const Cell& cell) {

  uint32_t sampleID = static_cast<uint32_t>(cell.id >> 32); // Extract the first 32 bits
  uint32_t cellID = static_cast<uint32_t>(cell.id & 0xFFFFFFFF); // Extract the second 32 bits
  
    os << sampleID << "\t"
       << cellID << "\t"
       << cell.cflag << "\t"
       << cell.pflag << "\t"      
       << cell.x << "\t"
       << cell.y << "\t";

    for (const auto &value : cell.cols) {
        os << value << "\t";
    }

    return os;
}

static void printValue(double value, int precision, int fixed_width) {

  double roundedValue = std::round(value * std::pow(10, precision)) / std::pow(10, precision);
  if (fixed_width > 0) {
    if (roundedValue == static_cast<int>(roundedValue)) {
      std::cout << std::setw(fixed_width) << static_cast<int>(roundedValue);
    } else {
      std::cout << std::setw(fixed_width) << std::fixed << std::setprecision(precision) << roundedValue;
    }
  } else {
    if (roundedValue == static_cast<int>(roundedValue)) {
      std::cout << static_cast<int>(roundedValue);
    } else {
      std::cout << std::fixed << std::setprecision(precision) << roundedValue;
    }
  }
}

/*void Cell::PrintWithHeader(int round, bool tabprint, const CellHeader& header) const {

  char d = tabprint ? '\n' : ',';

  uint32_t sampleID = static_cast<uint32_t>(id >> 32); // Extract the first 32 bits
  uint32_t cellID = static_cast<uint32_t>(id & 0xFFFFFFFF); // Extract the second 32 bits

  if (tabprint) {
    
    std::cout << std::fixed; // Set output to fixed. This applies to all floating-point output after this statement.
    std::cout << "sid:" << std::setw(8) << sampleID << d << "cid:" << std::setw(8) << cellID << d << "cflag:" << std::setw(8) << cflag <<
      d << "pflag:" << std::setw(8) << pflag << d << "x:";
    printValue(x, round, 8); // Assuming you modify printValue to use std::setw(8)
    std::cout << d << "y:";
    printValue(y, round, 8); // Assuming you modify printValue to use std::setw(8)
    
    // print cols
    size_t i = 0;
    for (const auto& t : header.GetDataTags()) {
      std::cout << d;
      std::cout << std::setw(8) << t.id << ":";
      printWithoutScientific(cols.at(i), round); // Assuming you modify printWithoutScientific to use std::setw(8)
      i++;
    }
  } else {

    std::cout << "sid:" << sampleID << d << "cid:" << cellID << d << "cflag:" <<cflag <<
      d << "pflag:" << pflag << d << "x:";
    printValue(x, round, -1);
    std::cout << d << "y:";
    printValue(y, round, -1);  
    
    // print cols
    size_t i = 0;
    for (const auto& t : header.GetDataTags()) {
      std::cout << d;
      std::cout << t.id << ":";
      printWithoutScientific(cols.at(i), round);
      i++;
    }

  }
  
    std::cout << std::endl;  
}

void Cell::Print(int round, bool tabprint) const {

  char d = tabprint ? '\t' : ',';

  uint32_t sampleID = static_cast<uint32_t>(id >> 32); // Extract the first 32 bits
  uint32_t cellID = static_cast<uint32_t>(id & 0xFFFFFFFF); // Extract the second 32 bits
  
  std::cout << sampleID << d << cellID << d << cflag << d << pflag << d;
  printValue(x, round);
  std::cout << d;
  printValue(y, round);  

  // print cols
  for (const auto& c : cols) {
    std::cout << d;
    printWithoutScientific(c, round);
  }

  std::cout << std::endl;
}
*/

void Cell::PrintWithHeader(int round,
			   bool tabprint,
			   bool header_print, 
			   const CellHeader& header) const {

  char d = tabprint ? '\t' : ',';
  uint32_t sampleID = static_cast<uint32_t>(id >> 32);
  uint32_t cellID = static_cast<uint32_t>(id & 0xFFFFFFFF);
  
  std::vector<std::string> headers(6,"");
  
  if (header_print)
    headers = {"sid:","cid:","cflag:","pflag:","x:","y:"};
  
  if (tabprint) std::cout << std::fixed; // Applies to all floating-point output after this statement
  
  int fixed_width = 6 + round;
  outputValue(headers[0], sampleID, tabprint, fixed_width, d);
  outputValue(headers[1], cellID, tabprint, fixed_width, d);
  outputValue(headers[2], cflag, tabprint, fixed_width, d);
  outputValue(headers[3], pflag, tabprint, fixed_width, d);
  std::cout << headers[4];
  printValue(x, round, tabprint ? fixed_width : -1);
  std::cout << d << headers[5];
  printValue(y, round, tabprint ? fixed_width : -1);

  // print cols
  size_t i = 0;
  for (const auto& t : header.GetDataTags()) {
    std::cout << (i > 0 || tabprint ? d : ','); // Conditional delimiter for the first column based on tabprint
    if (header_print)
      std::cout << t.id << ":";
    if (tabprint) std::cout << std::setw(fixed_width);
    printWithoutScientific(cols.at(i), round);
    i++;
  }
  
    std::cout << std::endl;
}

Cell::Cell(const std::string& row,
	   int x_index,// which column is X and Y
	   int y_index,
	   int start_index, // which columns start and end marker data
	   int end_index,
	   const CellHeader& header,
	   uint32_t cellid,
	   uint32_t sampleid) {
  
  const std::vector<std::string> tokens = tokenize_comma_delimited(row);  
  
  if (tokens.size() < 3) {
    throw std::runtime_error("CSV file should have at least three columns: id, x, y");
  }

  // check that we are in bounds
  if (x_index >= tokens.size() || y_index >= tokens.size() ||
      start_index >= tokens.size()) {
    std::cerr << "Error: cyftools convert - Not enough tokens for the header, line " << row << std::endl;
    assert(false);
  }
  
  pflag = 0;
  cflag = 0;
  id = 0;
  x = 0;
  y = 0;
  
  this->set_sample_id(sampleid);
  this->set_cell_id(cellid);
  
  // 1st and 2nd entry is x,y position
  x = std::strtof(tokens.at(x_index).c_str(), nullptr);
  y = std::strtof(tokens.at(y_index).c_str(), nullptr);

  // stop with as many columns as are in the header
  size_t num_cols = header.GetDataTags().size();

  // assume rest of the tags are column data
  size_t i = 0;
  while (i < tokens.size()) {

    // let user know if they are short on header specs
    /*    if (tokens.size() != m_last_line_count) {
      std::cerr << "Error: cyftools convert -- found "  << tokens.size() <<
	" on this line but " << m_last_line_count << " on last line" << std::endl;
      std::cerr << row << std::endl;

      for (const auto& t : header.GetDataTags())
	std::cerr << t.id << std::endl;
      assert(false);
      }*/

    // skip non-marker data and x and y
    if (i < start_index || i > end_index || i == x_index || i == y_index) {
      i++;
      continue;
    }

    cols.push_back(std::strtof(tokens.at(i).c_str(), nullptr));
    i++;
  }

  // let user know if they are short on data
  if (i < num_cols)
    std::cerr << "warning: only read in the available columns " << i << " but header specified " << num_cols << std::endl;
  
}
