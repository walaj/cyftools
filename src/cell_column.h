#ifndef CELL_COLUMN_H
#define CELL_COLUMN_H

#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include <typeinfo>
#include <type_traits>

#include "cy_utils.h"
#include "csv.h"

using namespace std;

/// @brief Alias for a row of cells, which can contain integers, floats, or strings.
using CellRow = vector<variant<int, float, string>>;

/// @brief Enum class representing the different types a Column can have.
enum class ColumnType {
    INT,
    FLOAT,
    STRING
};

/**
 * @class Column
 * @brief Base class for a column in a cell table, which can store elements of 
 *
 * This class provides an interface for columns of different types (integer, float, or string)
 * to be used in a table. Derived classes should implement the specific behavior for each type.
 */
class Column {
public:
    virtual ~Column() = default;

    virtual std::shared_ptr<Column> clone() const = 0;

    virtual std::shared_ptr<Column> CopyToFloat() const = 0;    
    
    virtual size_t size() const = 0;

    virtual std::string toString() const = 0;
    
    virtual ColumnType GetType() const = 0;

    virtual void Log10() = 0;

    virtual float GetNumericElem(size_t i) const = 0;

    virtual void SubsetColumn(const std::vector<size_t>& indices) = 0;

    virtual void PrintElem(size_t i) const = 0;

    virtual void SetPrecision(size_t n) = 0;

    virtual float Pearson(const Column& c) const = 0;

    virtual float Mean() const = 0;
    virtual float Min() const = 0;
    virtual float Max() const = 0;

    virtual void reserve(size_t n) = 0;
    
};

template <typename T>
class NumericColumn : public Column {
  static_assert(std::is_arithmetic<T>::value, "NumericColumn requires a numeric type");
  
 public:
  
  NumericColumn() {
    
    if (std::is_same_v<T, int>) {
      m_type = ColumnType::INT;
    } else if (std::is_same_v<T, float>) {
      m_type = ColumnType::FLOAT;
    } else {
      throw std::runtime_error("NumericColumn constructor: arg is neither an int nor a float");
    }
    
  }
  
  NumericColumn(const T& initial_elem) {
    
    if (std::is_same_v<T, int>) {
      m_type = ColumnType::INT;
    } else if (std::is_same_v<T, float>) {
      m_type = ColumnType::FLOAT;
    } else {
      throw std::runtime_error("NumericColumn constructor: arg is neither an int nor a float");      
    }
    
    m_vec.push_back(initial_elem);
  }
  
  std::shared_ptr<Column> clone() const override {
    return std::make_shared<NumericColumn<T>>(*this);
  }

  void reserve(size_t n) override {
    m_vec.clear();
    m_vec.reserve(n);
  }
  
  float GetNumericElem(size_t i) const override {
    if (i >= m_vec.size())
      throw std::out_of_range("Index out of range");
    return m_vec.at(i);
  }

  void SetNumericElem(T val, size_t i) {
    if (i >= m_vec.size())
      throw std::out_of_range("Index out of range");
    m_vec[i] = val;
  }

  float Pearson(const Column& c) const override {
    
    double mean_v1 = c.Mean();
    double mean_v2 = this->Mean();
    
    double num = 0.0, den_v1 = 0.0, den_v2 = 0.0;
    for (size_t i = 0; i < m_vec.size(); i++) {
      double val_v1 = c.GetNumericElem(i);
      double val_v2 = m_vec.at(i);
      
      num += (val_v1 - mean_v1) * (val_v2 - mean_v2);
      den_v1 += (val_v1 - mean_v1) * (val_v1 - mean_v1);
      den_v2 += (val_v2 - mean_v2) * (val_v2 - mean_v2);
    }
    
    return static_cast<float>(num / (std::sqrt(den_v1) * std::sqrt(den_v2)));
    
  }
  
  std::shared_ptr<Column> CopyToFloat() const override {
    
    auto currentColumn = static_cast<const NumericColumn<T>*>(this);
    std::shared_ptr<NumericColumn<float>> fcol = std::make_shared<NumericColumn<float>>(currentColumn->m_vec.size());

    for (size_t i = 0; i < currentColumn->m_vec.size(); ++i) {
      fcol->SetValueAt(i, static_cast<float>(currentColumn->m_vec[i]));
    }

    return fcol;
  }

  void SetValueAt(size_t index, const T& value) {
    m_vec[index] = value;
  }

  void PushElem(const T& elem) {
    m_vec.push_back(elem);
  }

  float Mean() const override {

    double nr = this->size();
    if (this->size() == 0) {
      throw std::runtime_error("Cannot compute the mean of an empty column.");
    }

    double sum = std::accumulate(m_vec.begin(), m_vec.end(), 0.0);
    double mean = sum / m_vec.size();
    
    return sum / nr;
  }
  
  float Min() const override {
    
    if (this->size() == 0) {
      throw std::runtime_error("Cannot compute the min of an empty column.");
    }
    return static_cast<float>(*std::min_element(m_vec.begin(), m_vec.end()));
  }
  
  float Max() const override {
    if (this->size() == 0) {
      throw std::runtime_error("Cannot compute the max of an empty column.");
    }
    
    return static_cast<float>(*std::max_element(m_vec.begin(), m_vec.end()));
  }
  
  size_t size() const override {
    return m_vec.size();
  }
  
  std::string toString() const override {
    const size_t print_lim = 3;
        std::stringstream ss;
        ss << "NumericColumn<" << typeid(T).name() << ">: Size: " <<
	  m_vec.size() << " [";
        for (size_t i = 0; i < std::min(size(), print_lim); i++) {
	  if (i > 0) ss << ", ";
	  ss << m_vec[i];
        }
        if (size() > print_lim) ss << ", ...";
        ss << "]";
        return ss.str();
    }

    void Log10() override {
      for (auto& elem : m_vec) {
	if (elem <= 0) {
	  throw std::invalid_argument("Log10 encountered a non-positive value");
	}
	elem = std::log10(elem);
      }
    }

    void SubsetColumn(const std::vector<size_t>& indices) override {
      std::vector<T> new_vec;
      new_vec.reserve(indices.size());
      for (const auto& index : indices) {
	if (index < 0 || index >= m_vec.size()) {
	  throw std::out_of_range("Index out of range");
	}
	new_vec.push_back(m_vec[index]);
      }
      m_vec = std::move(new_vec);
    }
    
    // Add this method to NumericColumn, StringColumn
    ColumnType GetType() const override {
      return m_type;
    }

    void PrintElem(size_t i) const override {
      
      if (i >= m_vec.size()) {
	throw std::runtime_error("Out of bounds in PrintElem");
	return;
      }
      
      switch (m_type) {
      case ColumnType::FLOAT:
        if (m_precision.has_value()) {
          std::cout << std::fixed << std::setprecision(m_precision.value()) <<
	    m_vec.at(i);
        } else {
          std::cout << m_vec.at(i);
        }
        break;
      case ColumnType::INT:
        std::cout << m_vec.at(i);
        break;
      default:
	throw std::runtime_error("Unexpected NumericColumn type");	
        break;
      }
    }
    
    void SetPrecision(size_t n) override {
      m_precision = n;
    }
    
protected:
    
    std::vector<T> m_vec;
    
    ColumnType m_type;

    std::optional<size_t> m_precision;
};

class StringColumn : public Column {
public:
    StringColumn(const std::string& initial_elem) {
        m_vec.push_back(initial_elem);
    }

    std::shared_ptr<Column> clone() const override {
      return std::make_shared<StringColumn>(*this);
    }
    
    ColumnType GetType() const override {
      return ColumnType::STRING; 
    }
    
    
    void PushElem(const std::string& elem) {
        m_vec.push_back(elem);
    }

    size_t size() const override {
        return m_vec.size();
    }

    std::string toString() const override {
        std::stringstream ss;
        ss << "StringColumn<string>: [";
        for (size_t i = 0; i < std::min(size(), (size_t)3); i++) {
            if (i > 0) ss << ", ";
            ss << m_vec[i];
        }
        if (size() > 3) ss << ", ...";
        ss << "]";
        return ss.str();
    }

    // Do nothing for StringColumn (or have dummies)
    void Log10() override {}
    float GetNumericElem(size_t i) const override { return 0; }
    float Pearson(const Column& c) const override { return 0; }
    float Mean() const override { return 0;   }
    float Min() const override { return 0;   }
    float Max() const override { return 0;   }
    void SetPrecision(size_t n) override {}
    
    std::shared_ptr<Column> CopyToFloat() const override {
      std::shared_ptr<NumericColumn<float>> fcol =
	std::make_shared<NumericColumn<float>>(1);
      return fcol;
    }

    void SubsetColumn(const std::vector<size_t>& indices) override {
      std::vector<std::string> new_vec;
      new_vec.reserve(indices.size());
      for (auto i : indices) {
	if (i >= 0 && i < m_vec.size()) {
	  new_vec.push_back(m_vec[i]);
	} else {
	  throw std::out_of_range("SubsetColumn: index out of range");
	}
      }
      m_vec = std::move(new_vec);
    }

    void PrintElem(size_t i) const override {

      if (i < m_vec.size()) {
	std::cout << m_vec.at(i);
      } else {
	throw std::runtime_error("Out of bounds in PrintElem");
      }
    }
    
    void reserve(size_t n) override {
      m_vec.clear();
      m_vec.reserve(n);
    }
    
 private:
    std::vector<std::string> m_vec;
};

#endif
