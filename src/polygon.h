#ifndef POLYGON_H
#define POLYGON_H

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
  std::vector<std::pair<float, float>> vertices;

  Polygon(int Id, const std::string& Name, const std::string& Text, const std::string& type,
	  const std::vector<std::pair<float, float>>& vertices)
    : Id(Id), Name(Name), Text(Text), type(type), vertices(vertices) {}
  
  friend std::ostream& operator<<(std::ostream& os, const Polygon& polygon);
  
  bool PointIn(float x, float y) const;

  size_t size() const { return vertices.size(); }
  
  // Non-const iterator overloads
  using iterator = std::vector<std::pair<float, float>>::iterator;
  iterator begin() { return vertices.begin(); }
  iterator end() { return vertices.end(); }
  
  // Const iterator overloads
  using const_iterator = std::vector<std::pair<float, float>>::const_iterator;
  const_iterator begin() const { return vertices.begin(); }
  const_iterator end() const { return vertices.end(); }
  const_iterator cbegin() const { return vertices.cbegin(); }
  const_iterator cend() const { return vertices.cend(); }
  
};

std::vector<Polygon> read_polygons_from_file(const std::string& file_path);

std::vector<std::pair<float, float>> parse_vertices(const std::string& vertex_str);


/*int main() {
  std::string file_path = "roi_file.txt";
  std::vector<Polygon> polygons = read_polygons_from_file(file_path);

  for (const Polygon& polygon : polygons) {
  std::cout << "Polygon ID: " << polygon.Id << ", Name: " << polygon.Name << ", Text: " << polygon.Text
  << ", Type: " << polygon.type << std::endl;

  for (const auto& vertex : polygon.vertices) {
  std::cout << "(" << vertex.first << ", " << vertex.second << ") ";
  }
  std::cout << std::endl;
  }

  return 0;
  }*/


#endif
