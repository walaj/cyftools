#include "cell_row.h"

#include "cell_utils.h"

static void printWithoutScientific(double number) {
    if (std::abs(number - std::round(number)) < 1e-9) {
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
    }
}


std::ostream& operator<<(std::ostream& os, const Cell& cell) {
    os << cell.m_id << "\t"
       << cell.m_cell_flag << "\t"
       << cell.m_pheno_flag << "\t"      
       << cell.m_x << "\t"
       << cell.m_y << "\t";

    for (const auto &value : cell.m_cols) {
        os << value << "\t";
    }

    os << "Graph size: " <<
      cell.m_spatial_ids.size() << " " <<
      cell.m_spatial_dist.size();
    return os;
}

void Cell::Print(int round) const {

  char d = ',';

  std::cout << m_id << d << m_cell_flag << d << m_pheno_flag << d;
  std::cout << std::setprecision(round) << m_x << d <<
    m_y;

  // print cols
  for (const auto& c : m_cols) {
    std::cout << d;
    printWithoutScientific(c);
    //    std::cout << std::fixed << d << c;
  }

  // print graph
  /////////

  // print delimiter if graph non-empty
  if (m_spatial_ids.size())
    std::cout << d;

  // [print
  for (size_t i = 0; i < m_spatial_ids.size(); i++) {
    std::cout << m_spatial_ids.at(i) << "^" << m_spatial_dist.at(i) <<
      "&" << m_spatial_flags.at(i);
    if (i != m_spatial_ids.size() - 1)
      std::cout << ";";
  }
  
  std::cout << std::endl;
}

Cell::Cell(const std::string& row, const CellHeader& header) {

  const std::vector<std::string> tokens = tokenize_comma_delimited(row);  

  if (tokens.size() < 3) {
    throw std::runtime_error("CSV file should have at least three columns: id, x, y");
  }

  m_pheno_flag = 0;
  m_cell_flag = 0;  
  
  // assume the first entry is cellid
  m_id = std::stoi(tokens.at(0));

  // assume the 2nd,3rd entry is x,y position
  m_x = std::strtof(tokens.at(1).c_str(), nullptr);
  m_y = std::strtof(tokens.at(2).c_str(), nullptr);

  // stop with as many columns as are in the header
  size_t cols = header.GetDataTags().size();
  
  // assume rest of the tags are column data
  size_t i = 0;
  while (i < tokens.size() - 3) {

    // let user know if they are short on header specs
    if (i >= cols) {
      std::cerr << "warning: more data columns " << tokens.size() << " in file than specified in header " << cols << std::endl;
      break;
    }
    m_cols.push_back(std::strtof(tokens.at(i + 3).c_str(), nullptr));
    i++;
  }

  // let user know if they are short on data
  if (i < cols)
    std::cerr << "warning: only read in the available columns " << i << " but header specified " << cols << std::endl;
  
}
