#include "cell_utils.h"

#include <string_view>
#include <stdexcept>
#include <limits>

int get_nth_element_as_integer(const std::string_view& str, size_t n) {
    size_t start = 0;
    size_t end = str.find(',');

    for (size_t index = 1; index < n; ++index) {
      if (end == std::string_view::npos) {
	throw std::out_of_range("The nth element is out of range.");
      }
      
      start = end + 1;
      end = str.find(',', start);
    }
    
    std::string_view token = (end != std::string_view::npos) ? str.substr(start, end - start) : str.substr(start);
    
    // Check if the token is an integer
    char* end_ptr;
    long value = std::strtol(token.data(), &end_ptr, 10);
    
    if (end_ptr != token.data() && static_cast<size_t>(end_ptr - token.data()) == token.size() && value >= std::numeric_limits<int>::min() && value <= std::numeric_limits<int>::max()) {
        return static_cast<int>(value);
    } else {
        throw std::invalid_argument("The nth element is not an integer.");
    }
}
