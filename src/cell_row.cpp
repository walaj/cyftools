#include "cell_row.h"

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

static void printValue(double value, int precision) {
    double roundedValue = std::round(value * std::pow(10, precision)) / std::pow(10, precision);
    if (roundedValue == static_cast<int>(roundedValue)) {
        std::cout << static_cast<int>(roundedValue);
    } else {
        std::cout << std::fixed << std::setprecision(precision) << roundedValue;
    }
}
void Cell::Print(int round) const {

  char d = ',';

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

  // print graph
  /////////

  // print delimiter if graph non-empty
  /* if (m_spatial_ids.size())
    std::cout << d;

  // [print
  for (size_t i = 0; i < m_spatial_ids.size(); i++) {
    std::cout << m_spatial_ids.at(i) << "^" << m_spatial_dist.at(i) <<
      "&" << m_spatial_flags.at(i);
    if (i != m_spatial_ids.size() - 1)
      std::cout << ";";
      }
  */
  
  std::cout << std::endl;
}

Cell::Cell(const std::string& row, const CellHeader& header,
	   uint32_t cellid,
	   uint32_t sampleid) {

  
  
  const std::vector<std::string> tokens = tokenize_comma_delimited(row);  

  if (tokens.size() < 3) {
    throw std::runtime_error("CSV file should have at least three columns: id, x, y");
  }

  pflag = 0;
  cflag = 0;
  id = 0;
  x = 0;
  y = 0;

  this->set_sample_id(sampleid);
  this->set_cell_id(cellid);

  // 1st and 2nd entry is x,y position
  x = std::strtof(tokens.at(0).c_str(), nullptr);
  y = std::strtof(tokens.at(1).c_str(), nullptr);

  // stop with as many columns as are in the header
  size_t num_cols = header.GetDataTags().size();
  
  // assume rest of the tags are column data
  size_t i = 0;
  while (i < tokens.size() - 2) {

    // let user know if they are short on header specs
    if (i >= num_cols) {
      std::cerr << "warning: more data columns " <<
	tokens.size() << " in file than specified in header " <<
	num_cols << std::endl;
      break;
    }
    cols.push_back(std::strtof(tokens.at(i + 2).c_str(), nullptr));
    i++;
  }

  // let user know if they are short on data
  if (i < num_cols)
    std::cerr << "warning: only read in the available columns " << i << " but header specified " << num_cols << std::endl;
  
}
