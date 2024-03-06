#include "polygon.h"

std::vector<JPoint> parse_vertices(const std::string& vertex_str) {
    std::vector<JPoint> vertices;
    std::istringstream vertex_stream(vertex_str);
    std::string coordinate_pair;

    while (std::getline(vertex_stream, coordinate_pair, ' ')) {
        std::istringstream coordinate_stream(coordinate_pair);
        std::string x_str, y_str;

        std::getline(coordinate_stream, x_str, ',');
        std::getline(coordinate_stream, y_str, ',');

        float x = std::stof(x_str);
        float y = std::stof(y_str);

        vertices.emplace_back(x, y);
    }

    return vertices;
}

// ray casting method to see if point is in a polygon
bool Polygon::PointIn(float x, float y) const {

  size_t n = vertices.size();
  bool inside = false;
  
  for (size_t i = 0, j = n - 1; i < n; j = i++) {
    if (((vertices[i].y > y) != (vertices[j].y > y)) &&
	(x < (vertices[j].x - vertices[i].x) * (y - vertices[i].y) / (vertices[j].y - vertices[i].y) + vertices[i].x)) {
      inside = !inside;
    }
  }
  
  return inside;
}

std::vector<Polygon> read_polygons_from_file(const std::string& file_path) {
    std::vector<Polygon> polygons;
    std::ifstream infile(file_path);

    if (!infile) {
        std::cerr << "Error opening roi file: " << file_path << std::endl;
        return polygons;
    }

    std::string line;
    std::getline(infile, line); // Skip the header line

    while (std::getline(infile, line)) {
      std::istringstream ss(line);
      
      std::string Id_str, Name, Text, type, all_points;
      std::getline(ss, Id_str, ',');
      
      // Read Name within quotes
      std::getline(ss, Name, '"');
      std::getline(ss, Name, '"');
      ss.ignore();
      //std::getline(ss, Name, ',');
      
      // Read Text within quotes
      std::getline(ss, Text, '"');
      std::getline(ss, Text, '"');
      ss.ignore();
      //std::getline(ss, Text, ',');

      std::getline(ss, type, '"');
      std::getline(ss, type, '"');
      ss.ignore();
      
      // Read all_points within quotes
      std::getline(ss, all_points, '"');
      std::getline(ss, all_points, '"');
      
      int Id = std::stoi(Id_str);
      std::vector<JPoint> vertices = parse_vertices(all_points);
      
      polygons.emplace_back(Id, Name, Text, type, vertices);

    }
    return polygons;
}

std::ostream& operator<<(std::ostream& os, const Polygon& polygon) {
  os << "ID: " << polygon.Id << " Name: " << polygon.Name << ", Text: " << polygon.Text << ", Type: " << polygon.type
       << ", Number of vertices: " << polygon.vertices.size() << std::endl;

    int vertices_to_print = (polygon.vertices.size() > 3) ? 3 : polygon.vertices.size();

    for (int i = 0; i < vertices_to_print; ++i) {
        const auto& vertex = polygon.vertices[i];
        os << "(" << vertex.x << ", " << vertex.y << ") ";
    }

    if (polygon.vertices.size() > 3) {
        os << "...";
    }

    return os;
}
