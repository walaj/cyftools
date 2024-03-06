#include "cell_table.h"
#include <delaunator.hpp>

#ifdef HAVE_CAIRO
#include "cairo/cairo.h"
#include "cairo/cairo-pdf.h"
#endif

// cgal is needed only for Voronoi diagram output
#ifdef HAVE_CGAL
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/draw_triangulation_2.h>
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Delaunay_triangulation_2<K> DelaunayData;
typedef K::Point_2 Point;
#endif

#ifdef HAVE_BOOST
#include <boost/functional/hash.hpp>
#endif


namespace std {
    template<> struct hash<JPoint> {
        std::size_t operator()(const JPoint& p) const {
            std::size_t hx = std::hash<float>{}(p.x);
            std::size_t hy = std::hash<float>{}(p.y);

#ifdef HAVE_BOOST
            std::size_t seed = 0;
            boost::hash_combine(seed, hx);
            boost::hash_combine(seed, hy);
            return seed;
#else
            return hx ^ (hy << 1);  // Shift hy 1 bit to the left and XOR it with hx.
#endif
        }
    };
}

// hash table structures for Delaunay and Voronoi (to keep from duplicating lines)
struct pair_hash {
    template <class T1, class T2>
    std::size_t operator () (const std::pair<T1,T2> &pair) const {
        return std::hash<T1>()(pair.first) ^ std::hash<T2>()(pair.second);
    }
};

struct jline_eq {
    bool operator() (const std::pair<JPoint, JPoint>& lhs, const std::pair<JPoint, JPoint>& rhs) const {
        return (lhs.first == rhs.first && lhs.second == rhs.second) ||
               (lhs.first == rhs.second && lhs.second == rhs.first);
    }
};

void CellTable::Delaunay(const std::string& pdf_delaunay,
		const std::string& pdf_voronoi,
	        int limit) {

#ifndef HAVE_BOOST
  std::cerr << "Warning: Will run but may be signifiacntly slower without including Boost library (uses Boost hash function for speed)" << std::endl;
#endif

  validate();

  float xmax = 0;
  float ymax = 0;  

  size_t ncells = CellCount();
  std::vector<double> coords(ncells * 2);
  
  for (size_t i = 0; i < ncells; i++) {
    coords[i*2  ] = m_x_ptr->getData().at(i);
    coords[i*2+1] = m_y_ptr->getData().at(i);
    if (coords[i*2] > xmax)
      xmax = coords[i*2];
    if (coords[i*2+1] > ymax)
      ymax = coords[i*2+1];
  }

#ifdef HAVE_CAIRO
  int width = xmax;
  int height = ymax;
#endif
  
  // construct the graph
  if (m_verbose) 
    std::cerr << "...constructing Delaunay triangulation on " << AddCommas(ncells) << " cells" << std::endl;
  delaunator::Delaunator d(coords);

  float micron_per_pixel = 0.325f;
  int limit_sq = limit <= 0 ? INT_MAX : limit * limit;
  
  // hash set to check if point already made
  std::unordered_set<std::pair<JPoint, JPoint>, pair_hash, jline_eq> lines;
  
  // reserve memory to avoid dynamic reallocation
  lines.reserve(d.triangles.size() / 3);
  
  size_t skip_count = 0;
  size_t draw_count = 0;
  
  // Find which lines to keep, by de-duplicating and removing lines > size limit
  for(size_t i = 0; i < d.triangles.size(); i+=3) {
    
    if (i/3 % 100000 == 0 && m_verbose)
      std::cerr << "...drawing triangle " << AddCommas(i/3) << " of " << AddCommas(d.triangles.size() / 3) << " for " <<
	AddCommas(ncells) << " cells" << std::endl;
    
    float x0 = static_cast<float>(d.coords[2 * d.triangles[i  ]]);
    float y0 = static_cast<float>(d.coords[2 * d.triangles[i  ] + 1]);
    float x1 = static_cast<float>(d.coords[2 * d.triangles[i+1]]);
    float y1 = static_cast<float>(d.coords[2 * d.triangles[i+1] + 1]);
    float x2 = static_cast<float>(d.coords[2 * d.triangles[i+2]]);
    float y2 = static_cast<float>(d.coords[2 * d.triangles[i+2] + 1]);

    std::array<std::pair<JPoint, JPoint>, 3> line_array = {      
      {std::make_pair(JPoint(x0, y0), JPoint(x1, y1)), 
       std::make_pair(JPoint(x1, y1), JPoint(x2, y2)), 
       std::make_pair(JPoint(x2, y2), JPoint(x0, y0))}
    };

    // deduplicate and remove Delaunay connections that are too short
    for (auto& line : line_array) {
      float dx = line.first.x - line.second.x;
      float dy = line.first.y - line.second.y;
      if (lines.find(line) == lines.end() && dx*dx + dy*dy <= limit_sq) {
        lines.insert(line);
        ++draw_count;
      } else {
        ++skip_count;
      }
    }
  }
  
  // store the adjaceny list
  std::unordered_map<JPoint, std::vector<JPoint>> adjList;  
  for (const auto& line : lines) {
    adjList[line.first].push_back(line.second);
    adjList[line.second].push_back(line.first); // assuming undirected graph
  }
  
  // build the adjaceny map
  std::unordered_set<JPoint> visited;
  std::unordered_map<JPoint, int> pointToComponentId;
  int currentComponentId = 1;
  
  // Iterative DFS using a stack
  std::stack<JPoint> stack;
  for (const auto& point : adjList) {
    if (visited.find(point.first) == visited.end()) {
      stack.push(point.first);
      while (!stack.empty()) {
	JPoint currentPoint = stack.top();
	stack.pop();
	if (visited.find(currentPoint) == visited.end()) {
	  visited.insert(currentPoint);
	  pointToComponentId[currentPoint] = currentComponentId;
	  for (const auto& neighbor : adjList[currentPoint]) {
	    stack.push(neighbor);
	  }
	}
      }
      currentComponentId++;
    }
  }
  
  // setup columns to store the delaunay components
  FloatColPtr d_label = std::make_shared<FloatCol>();
  FloatColPtr d_size  = std::make_shared<FloatCol>();
  
  // count the number of nodes for each component
  std::unordered_map<int, size_t> dcount;
  for (const auto& c : pointToComponentId) {
    dcount[c.second]++;
  }
  
  // fill the data into the columns
  for (size_t i = 0; i < m_x_ptr->size(); i++) {
    JPoint p = {m_x_ptr->getData().at(i),m_y_ptr->getData().at(i)};    
    
    // find the point in the component labels
    auto idd = pointToComponentId.find(p);

    // point not found, so delaunay edges already deleted
    if (idd == pointToComponentId.end()) {
      d_label->PushElem(0);
      d_size->PushElem(1);

    // point found
    } else {
      // set the label
      d_label->PushElem(idd->second);
      
      // set the label count
      assert(dcount.find(idd->second) != dcount.end());
      d_size->PushElem(dcount[idd->second]);
    }
  }
  
  // form the data tag
  Tag dtag_label(Tag::CA_TAG, "delaunay_component", "");
  AddColumn(dtag_label, d_label);
  Tag dtag_size(Tag::CA_TAG, "delaunay_count", "");
  AddColumn(dtag_size, d_size);
  
  // draw the Delaunay PDF
  if (!pdf_delaunay.empty()) {

    if (m_verbose)
      std::cerr << "...setting up for outputing Delaunay triangulation to PDF: " << pdf_delaunay << std::endl;

#ifdef HAVE_CAIRO
    
    cairo_surface_t *surface = cairo_pdf_surface_create(pdf_delaunay.c_str(), width, height);
    cairo_t *cr = cairo_create(surface);
    
    // Set the color of the lines to black.
    cairo_set_source_rgba(cr, 0, 0, 0, 1);
    cairo_set_line_width(cr, 0.3);

    // draw the actual lines
    for (const auto& line : lines) {
      //cairo_move_to(cr, line.first.x(), line.first.y());
      //      cairo_line_to(cr, line.second.x(), line.second.y());
      cairo_move_to(cr, line.first.x, line.first.y);
      cairo_line_to(cr, line.second.x, line.second.y);
      
    }

    cairo_stroke(cr); // Stroke all the lines
    
    cairo_new_path(cr); // Start a new path

    // setup a color map
    std::unordered_map<int, Color> color_map;
    for (auto c : pointToComponentId) {
      if (color_map.find(c.second) == color_map.end())
	color_map[c.second] = {rand() % 256,
			       rand() % 256,
			       rand() % 256};
    }

    // draw the points (cells) colored by component
    for (const auto& pair : pointToComponentId) {
      const JPoint& point = pair.first;
      int componentId = pair.second;
      Color c = color_map[componentId];

      cairo_set_source_rgb(cr, c.red/255.0, c.green/255.0, c.blue/255.0);
      cairo_arc(cr, point.x, point.y, 1.0, 0.0, 2.0 * M_PI);
      cairo_fill(cr);
    }
    
    if (m_verbose) 
      std::cerr << "...finalizing PDF of Delaunay triangulation on " << AddCommas(draw_count) <<
	" unique lines, skipping " << AddCommas(skip_count) << " lines for having length < " << limit << std::endl;
    
    // Clean up
    cairo_destroy(cr);
    cairo_surface_destroy(surface);
#else
    std::cerr << "Warning: Cairo PDF library needs to be linked during cysift build, no PDF will be output." << std::endl;
#endif
  }

  ///////////////
  // VORONOI
  ///////////////
  if (!pdf_voronoi.empty()) {

    if (m_verbose)
      std::cerr << "...setting up for outputing Voronoi diagram to PDF: " << pdf_voronoi << std::endl;

#if defined(HAVE_CAIRO) && defined(HAVE_CGAL)
    
    std::vector<Point> points(ncells);
    for (size_t i = 0; i < ncells; ++i) {
      points[i] = Point(m_x_ptr->getData().at(i), m_y_ptr->getData().at(i));
    }
    
    DelaunayData dt;
    dt.insert(points.begin(), points.end());

    // Create the Cairo surface and context
    cairo_surface_t *surface = cairo_pdf_surface_create(pdf_voronoi.c_str(), width, height); 
    cairo_t *cr = cairo_create(surface);

    // draw the original points
    for (const auto& p : points) {
      cairo_set_source_rgba(cr, 0.1, 0.0, 0.0, 1.0); // Set color to red
      cairo_arc(cr, p.x(), p.y(), 2.0, 0.0, 2.0 * M_PI);
      cairo_fill(cr);
    }

    // Set the color of the lines to black.
    cairo_set_source_rgba(cr, 0, 0, 0, 1);
    cairo_set_line_width(cr, 0.1);
    
    // draw the Voronoi edges    
    for (DelaunayData::Finite_edges_iterator eit = dt.finite_edges_begin(); eit != dt.finite_edges_end(); ++eit) {
      
      // The line segment between the circumcenters of the two triangles adjacent to each edge is an edge in the Voronoi diagram.
      DelaunayData::Face_handle f1 = eit->first;
      int adjacent_index = dt.mirror_index(f1, eit->second);
      DelaunayData::Face_handle f2 = f1->neighbor(adjacent_index);
      Point voronoi_edge_start = dt.circumcenter(f1);
      Point voronoi_edge_end = dt.circumcenter(f2);
      
      // Draw the Voronoi edge
      //cairo_move_to(cr, voronoi_edge_start.x(), voronoi_edge_start.y());
      //cairo_line_to(cr, voronoi_edge_end.x(), voronoi_edge_end.y());
    }

    //cairo_stroke(cr);
    
    // Finish the PDF
    if (m_verbose) 
      std::cerr << "...finalizing PDF of Voronoi diagram" << std::endl;
    
    cairo_show_page(cr);
    cairo_destroy(cr);
    cairo_surface_destroy(surface);
#else
    std::cerr << "Warning: Cairo PDF library and CGAL geometry libraries needs to be included/linked during cysift build for Voronoi output. Otherwise no PDF will be output." << std::endl;
#endif    
  }
  
  return;
  
}
