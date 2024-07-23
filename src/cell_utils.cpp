#include "cell_utils.h"

#include <string_view>
#include <stdexcept>
#include <limits>
#include <cmath>
#include <fstream>
#include <unordered_set>

#ifdef HAVE_CAIRO
#include "cairo/cairo.h"
#include "cairo/cairo-pdf.h"
#endif

void column_to_row_major(std::vector<float>& data, int nobs, int ndim) {

  float* temp = new float[data.size()];
  
  for (int row = 0; row < nobs; ++row) {
    for (int col = 0; col < ndim; ++col) {
      int old_index = col * nobs + row;
      int new_index = row * ndim + col;
      temp[new_index] = data[old_index];
    }
  }
  
  for (size_t i = 0; i < data.size(); ++i) {
    data[i] = temp[i];
  }
  
  delete[] temp;
  
}

PhenoMap phenoread(const std::string& filename) {
  
  PhenoMap data;
  
  // open the phenotype file
  // of format: marker(string), low_bound(float), upper_bound(float)
  std::ifstream file(filename);
  if (!file.is_open()) {
    std::cerr << "Error opening file: " << filename << std::endl;
    return data;
  }

  bool jeremiah = true;
  bool crevasse = false;// if crevasse format, handle differently
  bool onetable = false;// if one table format, handle differently

  // if one table format, store the marker names and order
  StringVec onetable_header;
  
  // read and load the phenotype file
  std::string line;
  size_t linenum = 0;
  while (std::getline(file, line)) {

    // check if crevasse format
    if (line.find("Gate") != std::string::npos && linenum == 0) {
      crevasse = true;
      continue;
    }

    else if (line.find("Sample") != std::string::npos && linenum == 0) {
      onetable = true;
      onetable_header = tokenize_comma_delimited<StringVec>(line);
      continue;
    }

    else {
      jeremiah = true;
    }

    std::istringstream lineStream(line);
    std::string key;
    float value1, value2;

    // jeremiah format
    if (jeremiah) {
      if (std::getline(lineStream, key, ',')) {
	lineStream >> value1;
	if (lineStream.get() == ',') {
	  lineStream >> value2;
	  data[key] = std::make_pair(value1, value2);
	}
      }
    }
    // one table format
    //else if (onetable) {
    //  if (std::
    //}
    
    // crevasse format
    else if (crevasse) {
      std::string token;
      std::getline(lineStream, token, ','); // discard first token
      std::getline(lineStream, token, ','); // read second token
      key = token.substr(1, token.length() - 2); // Remove the quotes
      std::getline(lineStream, token, ','); // read third token
      value1 = std::stof(token);
      value1 = std::exp(value1); // crevasse stores as log2
      value2 = 100000000;
      data[key] = std::make_pair(value1, value2);      
    } 
      
    linenum++;
  }
  
  file.close();
  return data;
  
}


std::string columnTypeToString(ColumnType type) {
  switch (type) {
    case ColumnType::INT:
      return "INT";
    case ColumnType::FLOAT:
      return "FLOAT";
    case ColumnType::STRING:
      return "STRING";
    case ColumnType::FLAG:
      return "FLAG";
    case ColumnType::GRAPH:
      return "GRAPH";
    default:
      return "UNKNOWN";
  }
}

#ifdef HAVE_HDF5
void write_hdf5_dataframe_attributes(H5::Group& group) {

  H5::DataSpace scalar_dataspace(H5S_SCALAR);
  H5::StrType strdatatype(H5::PredType::C_S1, H5T_VARIABLE);
    
  // write the word "_index"
  H5::Attribute attr = group.createAttribute("_index", strdatatype, scalar_dataspace);
  const char* temp_cstr1 = "_index";
  attr.write(strdatatype, &temp_cstr1);

  // write the encoding-type
  attr = group.createAttribute("encoding-type", strdatatype, scalar_dataspace);
  const char* temp_cstr2 = "dataframe";
  attr.write(strdatatype, &temp_cstr2);

  // write the encoding-version
  attr = group.createAttribute("encoding-version", strdatatype, scalar_dataspace);
  const char* temp_cstr3 = "0.1.0";
  attr.write(strdatatype, &temp_cstr3);
  
}
#endif

float euclidean_distance(float x1, float y1, float x2, float y2) {
    float dx = x2 - x1;
    float dy = y2 - y1;
    return std::sqrt(dx*dx + dy*dy);
}

float euclidean_distance_squared(float x1, float y1, float x2, float y2) {
    float dx = x2 - x1;
    float dy = y2 - y1;
    return std::sqrt(dx*dx + dy*dy);
}

bool check_readable(const std::string& filename) {

  std::ifstream file(filename);

  bool answer = file.is_open();
  
  // Close the file if it was opened
  if (file) {
    file.close();
  }

  return answer;
}

std::pair<std::string, std::string> colon_parse(const std::string& str) {

    auto colonPos = str.find(':');
    if (colonPos == std::string::npos || colonPos == 0 || colonPos == str.length() - 1) {
        throw std::invalid_argument("Invalid string format");
    }

    auto colonCount = std::count(str.begin(), str.end(), ':');
    if (colonCount > 1) {
        throw std::invalid_argument("Multiple colons are not allowed");
    }

    return std::make_pair(str.substr(0, colonPos), str.substr(colonPos + 1));
}


void draw_scale_bar(cairo_t* cr, double x, double y, double bar_width, double bar_height, const std::string& text) {

#ifdef HAVE_CAIRO
    // Set color to white for the bar
    cairo_set_source_rgb(cr, 1.0, 1.0, 1.0); // White color
    cairo_set_line_width(cr, bar_height); // Set bar height as line width
    
    // Draw the bar
    cairo_move_to(cr, x, y);
    cairo_line_to(cr, x + bar_width, y);
    cairo_stroke(cr); // Apply the drawing

    // Prepare to draw text
    cairo_set_source_rgb(cr, 255, 255, 255); // Text color (black for contrast)
    cairo_select_font_face(cr, "Arial", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_BOLD);
    cairo_set_font_size(cr, bar_height*2); // Set font size relative to bar height

    // Calculate text width and height for centering
    cairo_text_extents_t extents;
    cairo_text_extents(cr, text.c_str(), &extents);
    double text_x = x + (bar_width - extents.width) / 2 - extents.x_bearing;
    double text_y = y - (bar_height + extents.height) / 2 - extents.y_bearing - bar_width * 0.1;

    // Draw the text
    cairo_move_to(cr, text_x, text_y);
    cairo_show_text(cr, text.c_str());
    cairo_stroke(cr); // Apply the drawing
#else
    std::cerr << "Warning - unable to draw scale bar, need to link Cairo library" << std::endl;
#endif    
}


/*void add_legend_cairo(cairo_t* crp, int font_size,
		      int legend_width, int legend_height,
		      int legend_x, int legend_y,
		      const ColorLabelMap& cm) {
#ifdef HAVE_CAIRO
  // Setting up font face
  cairo_select_font_face(crp, "Arial",
			 CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_BOLD);
  cairo_set_font_size(crp, font_size); // Adjust font size to your needs

  // Dimensions for the legend
  int legend_padding = 20; // Space between color boxes
  
  for (int i = 0; i < cm.size(); ++i) {
    
    // Set color for this entry
    Color c = cm[i].first;
    cairo_set_source_rgb(crp, c.redf(), c.greenf(), c.bluef());
    
    // Draw the color box
    cairo_rectangle(crp, legend_x, legend_y + i*(legend_height+legend_padding), legend_width, legend_height);
    cairo_fill(crp);
    
    // Draw the label
    cairo_set_source_rgb(crp, 0, 0, 0); // Set color to black for the text
    cairo_move_to(crp, legend_x, legend_y + i*(legend_height+legend_padding) + legend_height);
    cairo_show_text(crp, cm[i].second.c_str());
  }
#else
  std::cerr << "Warning: Can't plot, need to build with Cairo library" << std::endl;
#endif  
  
}
*/

void add_legend_cairo_top(cairo_t* crp, int font_size,
                      int legend_height,
                      int width,
			  const ColorLabelVec& cm) {
#ifdef HAVE_CAIRO
  
  // Black background for legend area
  cairo_set_source_rgb(crp, 0, 0, 0); // Set color to black
  cairo_rectangle(crp, 0, 0, width, legend_height); // Top strip
  cairo_fill(crp);

  // Measure the total width of the labels and spaces
  double total_width = 0;
  const std::string space = "  "; // Two spaces
  cairo_text_extents_t extents;
  cairo_set_font_size(crp, font_size); // Adjust font size to your needs  
  for (const auto& entry : cm) {
    if (entry.label == "nolabel") // keyword to skip the label
      continue;
    std::string label = entry.label + space; // Add spaces between labels
    cairo_text_extents(crp, label.c_str(), &extents);
    total_width += extents.x_advance; // Use x_advance for accurate spacing
  }
  
  // Subtract the width of the final two spaces, as they're not needed after the last label
  cairo_text_extents(crp, space.c_str(), &extents);
  total_width -= extents.x_advance;
  
  // Calculate starting x position to center the block of text
  int start_x = (width - total_width) / 2;
  int text_x = start_x;
  int text_y = legend_height / 2 + font_size / 2; // Adjust as needed for vertical centering
  
  // Draw the labels
  for (const auto& entry : cm) {
    if (entry.label == "nolabel") // keyword to skip the label
      continue;

    Color c = entry.c;
    std::string label = entry.label;
    
    // Set color for text
    cairo_set_source_rgb(crp, c.redf(), c.greenf(), c.bluef());
    
    // Draw the label
    cairo_move_to(crp, text_x, text_y);
    cairo_show_text(crp, label.c_str());
    
    // Measure the label width with the two spaces and update text_x for the next label
    cairo_text_extents(crp, (label + space).c_str(), &extents);
    text_x += extents.x_advance;
  }
#else
  std::cerr << "Warning: Can't plot, need to build with Cairo library" << std::endl;
#endif  
}
  
double jaccardSimilarity(const std::vector<bool>& v1, const std::vector<bool>& v2) {
    int intersectionCount = 0;
    int unionCount = 0;

    for (size_t i = 0; i < v1.size(); i++) {
        if (v1[i] && v2[i]) {
            intersectionCount++;
        }
        if (v1[i] || v2[i]) {
            unionCount++;
        }
    }

    if (unionCount == 0) return 0.0; // to avoid division by zero
    
    return static_cast<double>(intersectionCount) / unionCount;
}

double jaccardSubsetSimilarity(const std::vector<bool>& v1, const std::vector<bool>& v2) {
    int intersectionCount = 0;
    int minSizeCount = std::min(v1.size(), v2.size());

    for (size_t i = 0; i < v1.size(); i++) {
        if (v1[i] && v2[i]) {
            intersectionCount++;
        }
    }

    if (minSizeCount == 0) return 0.0; // to avoid division by zero
    
    return static_cast<double>(intersectionCount) / minSizeCount;
}

////////
/// adapted from: https://www.geeksforgeeks.org/convex-hull-using-graham-scan/
// Utility function to find the point with the lowest Y coordinate (or the leftmost in case of tie)
bool compareYThenX(const JPoint& p1, const JPoint& p2) {
    return (p1.y < p2.y) || (p1.y == p2.y && p1.x < p2.x);
}

// To find orientation of ordered triplet (p, q, r).
int orientation(const JPoint& p, const JPoint& q, const JPoint& r) {
    int val = (q.y - p.y) * (r.x - q.x) - (q.x - p.x) * (r.y - q.y);
    if (val == 0) return 0; // Collinear
    return (val > 0) ? 1 : 2; // Clockwise or Counterclockwise
}

// Comparator for sorting points with respect to the first point
bool compare(const JPoint& p1, const JPoint& p2, const JPoint& p0) {
    int o = orientation(p0, p1, p2);
    if (o == 0)
        return (p1.x - p0.x) * (p1.x - p0.x) + (p1.y - p0.y) * (p1.y - p0.y) < (p2.x - p0.x) * (p2.x - p0.x) + (p2.y - p0.y) * (p2.y - p0.y);
    return (o == 2);
}

std::vector<JPoint> convexHull(std::vector<JPoint>& points) {

  // returning since can't compute so small
  if (points.size() < 3)
    return std::vector<JPoint>();
  
  //std::cerr << "finding bottomost point" << std::endl;
  // Step 1: Find the bottommost point
  std::swap(points[0], *std::min_element(points.begin(), points.end(), compareYThenX));
  JPoint p0 = points[0];

  //std::cerr << "sorting polar" << std::endl;
  // Step 2: Sort points based on polar angle with p0
  std::sort(points.begin() + 1, points.end(), [p0](const JPoint& a, const JPoint& b) { return compare(a, b, p0); });

  //std::cerr << "removing colinear points" << std::endl;
  // Step 3: Remove points that are collinear with p0
  int m = 1; // Initialize size of modified array
  for (int i = 1; i < points.size(); ++i) {
    while (i < points.size() - 1 && orientation(p0, points[i], points[i+1]) == 0) i++;
    points[m++] = points[i];
  }
  
  // Convex hull is not possible if there are less than 3 unique points
  if (m < 3) return std::vector<JPoint>();
  
  // Step 4: Create an empty stack and push first three points to it
  std::stack<JPoint> S;
  S.push(points[0]);
  S.push(points[1]);
  S.push(points[2]);
  
  // Step 5: Process remaining points
  //std::cerr << "processing remaining points" << std::endl;
  for (int i = 3; i < m; i++) {
    while (S.size() > 1 && orientation(nextToTop(S), S.top(), points[i]) != 2) S.pop();
    S.push(points[i]);
  }
  
  std::vector<JPoint> hull;
  
  // Now stack has the output points, print contents of stack
  //std::cerr << "emptying stack" << std::endl;
  while (!S.empty()) {
    JPoint p = S.top();
    hull.push_back(p);
    //std::cout << "(" << p.x << ", " << p.y << ")" << std::endl;
    S.pop();
  }
  return hull;
}


// A utility function to find next to top in a stack of JPoints
JPoint nextToTop(std::stack<JPoint> &S) {
    JPoint p = S.top();
    S.pop();
    JPoint res = S.top();
    S.push(p);
    return res;
}
///////////////
//////////////

std::string clean_marker_string(const std::string& input) {
    // Copy the input string to work with
    std::string result = input;
    
    // Remove hyphens
    result.erase(std::remove(result.begin(), result.end(), '-'), result.end());
    
    // Find the position of the last underscore
    /*    size_t underscorePos = result.rfind('_');
    if (underscorePos != std::string::npos) {
      // Remove everything from the last underscore to the end
      result.erase(underscorePos);
    }
    */
    return result;
}
