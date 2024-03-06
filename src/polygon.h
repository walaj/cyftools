#pragma once

// for JPoint
#include "cell_utils.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <utility>

class Polygon {
  
public:
  
  int Id;
  std::string Name;
  std::string Text;
  std::string type;
  std::vector<JPoint> vertices;

  Polygon(int Id, const std::string& Name,
	  const std::string& Text,
	  const std::string& type,
	  const std::vector<JPoint>& vertices)
    : Id(Id), Name(Name), Text(Text), type(type), vertices(vertices) {}
  
  Polygon(const std::vector<JPoint>& vec); 
  
  friend std::ostream& operator<<(std::ostream& os, const Polygon& polygon);
  
  bool PointIn(float x, float y) const;

  size_t size() const { return vertices.size(); }
  
  // Non-const iterator overloads
  using iterator = std::vector<JPoint>::iterator;
  iterator begin() { return vertices.begin(); }
  iterator end() { return vertices.end(); }
  
  // Const iterator overloads
  using const_iterator = std::vector<JPoint>::const_iterator;
  const_iterator begin() const { return vertices.begin(); }
  const_iterator end() const { return vertices.end(); }
  const_iterator cbegin() const { return vertices.cbegin(); }
  const_iterator cend() const { return vertices.cend(); }
  
};

std::vector<Polygon> read_polygons_from_file(const std::string& file_path);

std::vector<JPoint> parse_vertices(const std::string& vertex_str);
