#pragma once

#include <vector>
#include <string>
#include <sstream>
#include <iomanip>
#include <numeric>
#include <iostream>
#include "color_map.h"
#include <stack>

// forward declare
struct _cairo;
typedef _cairo cairo_t;

#include "cell_column.h"

#ifdef HAVE_HDF5
#include <H5Cpp.h>
#endif

struct JPoint {
  
  float x;
  float y;

  JPoint(float mx, float my) : x(mx), y(my) {}
  
  std::string print() const { return std::to_string(x) + "," + std::to_string(y); }

  bool operator==(const JPoint& other) const {
    return x == other.x && y == other.y;
  }

  friend std::ostream& operator<<(std::ostream& os, const JPoint& point) {
    os << point.x << "," << point.y;
    return os;
}

};


/** Format an integer to include commas
 * @param data Number to format
 * @return String with formatted number containing commas
 */
template <typename T> inline std::string AddCommas(T data) {
  std::stringstream ss; 
  ss << std::fixed << data; 
  std::string s = ss.str();
  if (s.length() > 3)
    for (int i = s.length()-3; i > 0; i -= 3)
      s.insert(i,",");
  return s;
}

template<typename T>
T mean(const std::vector<T>& v) {
  static_assert(std::is_arithmetic<T>::value, "Type must be numerical!");
  
  if (v.empty()) return T(0);
  return std::accumulate(v.begin(), v.end(), T(0)) / v.size();
}

template<typename T>
float pearsonCorrelation(const std::vector<T>& v1, const std::vector<T>& v2) {
  double mean_v1 = mean(v1);
  double mean_v2 = mean(v2); 
  
  double num = 0.0, den_v1 = 0.0, den_v2 = 0.0;
  for (size_t i = 0; i < v1.size(); i++) {
    double val_v1 = v1.at(i);
    double val_v2 = v2.at(i);
    
    num += (val_v1 - mean_v1) * (val_v2 - mean_v2);
    den_v1 += (val_v1 - mean_v1) * (val_v1 - mean_v1);
    den_v2 += (val_v2 - mean_v2) * (val_v2 - mean_v2);
    
  }
  
  return static_cast<float>(num / (std::sqrt(den_v1) * std::sqrt(den_v2))); 
}

///////
/////// CONVEX HULL
///////

bool compareYThenX(const JPoint& p1, const JPoint& p2);

int orientation(const JPoint& p, const JPoint& q, const JPoint& r);

bool compare(const JPoint& p1, const JPoint& p2, const JPoint& p0);

std::vector<JPoint> convexHull(std::vector<JPoint>& points);

JPoint nextToTop(std::stack<JPoint> &S);

///////

// Comparator function for sorting
bool ComparePoints(const JPoint& a, const JPoint& b, const JPoint& centroid);

double jaccardSubsetSimilarity(const std::vector<bool>& v1, const std::vector<bool>& v2);

double jaccardSimilarity(const std::vector<bool>& v1, const std::vector<bool>& v2);

//////
/// plotting
void add_legend_cairo_top(cairo_t* crp, int font_size,
			  int legend_height,
			  int width, const ColorLabelMap& cm);


void add_legend_cairo(cairo_t* crp, int font_size,
		      int legend_width, int legend_height,
		      int legend_x, int legend_y,
		      const ColorLabelMap& cm);

void draw_scale_bar(cairo_t* cr, double x, double y, double bar_width, double bar_height, const std::string& text);

void column_to_row_major(std::vector<float>& data, int nobs, int ndim);

PhenoMap phenoread(const std::string& filename);

std::vector<std::string> tokenize_comma_delimited(const std::string& str);

float euclidean_distance(float x1, float y1, float x2, float y2);

bool check_readable(const std::string& filename);

float euclidean_distance_squared(float x1, float y1, float x2, float y2);

std::pair<std::string, std::string> colon_parse(const std::string& str);

bool is_mcmicro_meta(const std::string& str);

std::string clean_marker_string(const std::string& input);

#ifdef HAVE_HDF5
void write_hdf5_dataframe_attributes(H5::Group& group);
#endif

enum class ColumnType;
std::string columnTypeToString(ColumnType type);
