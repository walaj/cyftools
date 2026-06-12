#include "cell_table.h"
#include "color_map.h"
#include "cell_selector.h"
#include <stack>


#ifdef HAVE_CAIRO
#include "cairo/cairo.h"
#include "cairo/cairo-pdf.h"
#endif

// block for attempting bufferered writing to improve paralellization
// not really doing much, CELL_BUFFER_LIMIT 1 basically turns this off,
// but still here in case useful later
#define CELL_BUFFER_LIMIT 1

//#ifdef __clang__
//std::vector<Cell> cell_buffer;
//#pragma omp threadprivate(cell_buffer)
//#endif

int CellTable::dfs(int cellIndex,
		   std::vector<int>& component_label,
		   int currentLabel,
		   const std::vector<std::vector<size_t>>& neighbors) const {

  std::stack<int> stack;
  stack.push(cellIndex);
  
  int size = 0;
  
  while (!stack.empty()) {
    int curr = stack.top();
    stack.pop();

    if (component_label[curr] == 0) {
      component_label[curr] = currentLabel;
      size++;
      
      for (int neighbor : neighbors[curr]) {
	if (component_label[neighbor] == 0) {
	  stack.push(neighbor);
	}
      }
    }
  }
  
  return size;
}

void CellTable::IslandFill(size_t n, int flag_from, bool invert_from,
			   int flag_to, bool invert_to) {
  
  validate();

  if (flag_from == 0) {
    std::cerr << "Error - cyftools island - flag from is zero, needs to be non-zero" << std::endl;
    return;
  } 
  
  if (flag_to == 0) {
    std::cerr << "Error - cyftools island - flag to is zero, needs to be non-zero" << std::endl;
    return;
  } 

  const size_t NUM_NEIGHBORS = 10;
  const float DIST_LIMIT = 100;
  
#ifdef HAVE_KNNCOLLE  

  size_t num_cells = CellCount();

  // build the vp tree for the entire network
  knncolle::VpTree<knncolle::distances::Euclidean, int, float> searcher_all
    = build_vp_tree();

  // for every cell, determine if it has flag-set neighbor
  if (m_verbose)
    std::cerr << "...getting whole slide nearest neighbors." << std::endl;

  std::vector<bool> marked_neighbor(num_cells, false);
  for (size_t i = 0; i < num_cells; i++) {

    // find the k nearest neighbors
    JNeighbors neigh = searcher_all.find_nearest_neighbors(i, 10);

    for (const auto& nn : neigh) {

      // is the neighbor a marked cell or not
      bool marked = IS_FLAG_SET(m_cflag_ptr->at(nn.first), flag_from);
      marked = invert_from ? !marked : marked;
      
      // then if within dist limited, neighbor is marked
      if (nn.second < DIST_LIMIT && marked) {
	marked_neighbor[i] = true;
	break;
      }
    }
  }

  // map to move from indicies in stroma to global indicies
  std::vector<int> stroma_to_all_index_map;
  std::vector<int> all_to_stroma_index_map(num_cells, -1);
  std::vector<bool> ix(num_cells, false);
  size_t stroma_index = 0;
  for (size_t i = 0; i < num_cells; i++) {

    // find if cell is marked
    bool marked = IS_FLAG_SET(m_cflag_ptr->at(i), flag_from);
    marked = invert_from ? !marked : marked;
    
    // if its NON-marked cell (often, stromal) cell,
    // 1) add the index from where would be if counted all cells
    // 2) add the stroma index to the all map
    if (!marked) { //IS_FLAG_SET(m_cflag_ptr->at(i), MARK_FLAG)) {
      ix[i] = true;
      stroma_to_all_index_map.push_back(i);
      all_to_stroma_index_map[i] = stroma_index;      
      stroma_index++;      
    }
  }
  
  const int num_stroma = stroma_index;
  assert(stroma_to_all_index_map.size() == num_stroma);
  
  // store the neighbors in stromal coordinates
  std::vector<std::vector<size_t>> neighbors;
  neighbors.resize(num_stroma);

  // build the vp tree for the stroma
  knncolle::VpTree<knncolle::distances::Euclidean, int, float> searcher = build_vp_tree(ix);
  
  // from the VP tree, build the nearest neighbors indicies vector
  if (m_verbose)
    std::cerr << "...getting unmarked-only nearest neighbors" << std::endl;
  
  stroma_index = 0;
  for (size_t i = 0; i < num_stroma; i++) {

    JNeighbors neigh = searcher.find_nearest_neighbors(i, NUM_NEIGHBORS);
    neighbors[i].reserve(NUM_NEIGHBORS);

    // make the stromal (unmarked) nearest neighbors graph
    for (const auto& nn : neigh) {
      assert(nn.first < num_stroma);
      if (nn.second < DIST_LIMIT)
	neighbors[i].push_back(nn.first);
    }
  }
  
  // Initialize connected component labels
  std::vector<int> component_label(num_stroma, 0);
  int currentLabel = 1;
  
  if (m_verbose)
    std::cerr << "...looping cells to fill unmarked islands with fewer than " <<
      AddCommas(n) << " cells" << std::endl;
  
  // loop the cells and start DFS connected component algorithm on stromal cells
  std::unordered_set<int> island_count;
  size_t filled_cell_count = 0;
  stroma_index = 0;
  for (int i = 0; i < num_cells; i++) {
    
    // find if cell is marked
    bool marked = IS_FLAG_SET(m_cflag_ptr->at(i), flag_from);
    marked = invert_from ? !marked : marked;
    
    // skip marked cells
    if (marked) // || already_visited) // IS_FLAG_SET(m_cflag_ptr->at(i), MARK_FLAG))
      continue;

    // find index of cell in table (all) coords
    int all_index = stroma_to_all_index_map.at(stroma_index);
    assert(all_index == i);

    // if not yet part of a component, do the DFS
    if (component_label[stroma_index] == 0) {
      int size = dfs(stroma_index, component_label, currentLabel, neighbors);
      
      //      std::cerr << "...DFS " << currentLabel << " with " << AddCommas(size) << " components " <<
      //	" and marked neighbor: " << marked_neighbor.at(all_index) << std::endl;

      // store which connected components should get flipped
      if (size < n && marked_neighbor.at(i)) {
	island_count.insert(currentLabel);
	filled_cell_count += size;
      }
      currentLabel++;
    }

    stroma_index++;
  }

  // which connected components have a tumor connection?
  std::unordered_set<size_t> marked_connect_component;
  for (size_t i = 0; i < num_stroma; i++) {
    int all_index = stroma_to_all_index_map.at(i);
    if (marked_neighbor[all_index])
      marked_connect_component.insert(component_label[i]);
  }

  // go through and actually change the flags
  // remember we are looping here in all coords but components are in stromal coords
  for (size_t i = 0; i < num_cells; i++) {

    int stromal_ix = all_to_stroma_index_map.at(i);
    if (island_count.count(component_label[stromal_ix]) &&
	marked_connect_component.count(component_label[stromal_ix])) {
      if (invert_to)
	CLEAR_FLAG((*m_cflag_ptr)[i], flag_to);
      else
	SET_FLAG((*m_cflag_ptr)[i], flag_to);	
    }
  }
  
  if (m_verbose) {
    std::cerr << "...filled " << AddCommas(island_count.size()) << " unmarked \"islands\" involving " <<
      AddCommas(filled_cell_count) << " cells" << std::endl;
  }
  
#else
  std::cerr << "Warning: tumor call function requires including header library knncolle (https://github.com/LTLA/knncolle)" <<
    " and preprocessor directive -DHAVE_KNNCOLLE" << std::endl;
#endif
  
}

#ifdef HAVE_KNNCOLLE
knncolle::VpTree<knncolle::distances::Euclidean, int, float>
CellTable::build_vp_tree(const std::vector<bool>& ix) const {

  // get the number of cells to build table on
  validate();
  const int ndim = 2;
  int nobs = CellCount(); // max number of cells we'll use in graph
  assert(nobs == ix.size() || ix.size() == 0);
  
  // column major the coordinate data
  std::vector<float> concatenated_data;
  concatenated_data.reserve(nobs * ndim);

  // add x and y
  size_t n = 0;  
  for (size_t i = 0; i < nobs; i++) {
    if (ix.size() == 0 || ix.at(i)) {
      concatenated_data.push_back(m_x_ptr->at(i));
      concatenated_data.push_back(m_y_ptr->at(i));      
      n++;
    }
  }

  // update the number of obs (cells) based on how many we actually added
  assert(n == concatenated_data.size() / 2);
  nobs = n; 
  
  if (m_verbose) 
    std::cerr << "...building KNN VP tree with " << AddCommas(concatenated_data.size()/2) << " points" << std::endl;
  
  JNeighbors output(nobs);
   
  // initialize the tree. Can choose from different algorithms, per knncolle library
  knncolle::VpTree<knncolle::distances::Euclidean, int, float> searcher(ndim, nobs, concatenated_data.data());
  //knncolle::AnnoyEuclidean<int, float> searcher(ndim, nobs, concatenated_data.data());
  //knncolle::Kmknn<knncolle::distances::Euclidean, int, float> searcher(ndim, nobs, concatenated_data.data());  

  if (m_verbose)
    std::cerr << "..KNN (spatial) vp tree built" << std::endl;
  
   return searcher;
   
}
#endif

void CellTable::Distances(const std::string& id) {

#ifdef HAVE_KNNCOLLE
  
  validate();

  int n = CellCount();
  
  // find which cells are "marked"
  std::vector<bool> ix(n, false);
  for (size_t i = 0; i <  n; i++) {
    if (IS_FLAG_SET(m_cflag_ptr->at(i), MARK_FLAG))
      ix[i] = true;
  }

  // clear the mark
  ClearCFlag(MARK_FLAG);

  // ensure there are some that are set
  int marked_cells = std::count(ix.begin(), ix.end(), true);
  if (marked_cells == 0) {
    std::cerr << "Warning - cyftools distance on id " << id << " - No marked cells, did you run cyftools filter -M first?" << std::endl;
    return;
  }
  if (m_verbose)
    std::cerr << "...cyftools dist - Finding distance to " << marked_cells << " marked cells and storing in label " << id << std::endl;
   
  // build the tree
  knncolle::VpTree<knncolle::distances::Euclidean, int, float> searcher = build_vp_tree(ix);

  // store the distances
  auto dist_ptr = std::make_shared<FloatCol>();
  dist_ptr->resize(n);
  
  // find the nearest neighbor
  for (size_t i = 0; i < n; i++) {

    if (i % 100000 == 0 && m_verbose) {
      std::cerr << "...cyftools dist - getting dist for cell " << i << std::endl;
    }

    // if its already marked, dist 0 by definition
    if (IS_FLAG_SET(m_cflag_ptr->at(i), MARK_FLAG)) {
      dist_ptr->SetNumericElem(0, i);
      continue;
    }
    
    float x = m_x_ptr->at(i);
    float y = m_y_ptr->at(i);
    float point[2] = {x, y};
    JNeighbors neigh = searcher.find_nearest_neighbors(point, 1);
    assert(neigh.size() == 1);

    // set the distance
    dist_ptr->SetNumericElem(neigh.at(0).second, i);
    
  }
  
  // add the columm
  Tag dtag(Tag::CA_TAG, "dist_" + id, "");
  AddColumn(dtag, dist_ptr);
  
#else
  std::cerr << "Warning - not able to calculate without a KNN tree implementation linked" << std::endl;
#endif  
}

void CellTable::CallTLS(cy_uint bcell_marker, cy_uint immune_marker,
			int min_cluster_size, int dist_max) {

  if (bcell_marker == 0) {
    std::cerr << "Error - cyftools tls - bcell_marker is zero, needs to be non-zero" << std::endl;
    return;
  } 
  if (immune_marker == 0) {
    std::cerr << "Error - cyftools tls - immune_marker is zero, needs to be non-zero" << std::endl;
    return;
  } 

  
  if (m_verbose)
    std::cerr << "...cyftools tls - bcell marker " << bcell_marker << " immune marker(s) " <<
      immune_marker << " min cluster size " << min_cluster_size << " dist max " << dist_max << std::endl;
  
  validate();

  int n = CellCount();
  
  // Step 1) Identify B-cells and mark them
  if (m_verbose) { std::cerr << "...cyftools tls - marking B-cells (flag = " << bcell_marker << ")" << std::endl; }
  for (size_t i = 0; i < n; i++) {

    cy_uint& cf  = (*m_cflag_ptr)[i];
    cy_uint pf  = m_pflag_ptr->at(i);    
    
    // clear the mark (shouldn't be set anyway)
    CLEAR_FLAG(cf, MARK_FLAG);
    
    // if the b-cell flag is set, set the mark
    if (IS_FLAG_SET(pf, bcell_marker)) {
      SET_FLAG(cf, MARK_FLAG);
    }
  }
  if (m_verbose) { std::cerr << "...cyftools tls - " << CountCFlag(MARK_FLAG) << " B-cells using " << bcell_marker << std::endl; }
  
  // Step 2) Find cells that are surrounded by B-cells, mark as TLS
  int num_neighbors = 25;
  float frac_pos = 0.5;
  if (m_verbose) { std::cerr << "...cyftools tls - KNN to find cells near B-cells" << std::endl; }  
  ClearCFlag(TLS_FLAG); // start with all TLS_FLAG off
  AnnotateCall(num_neighbors, frac_pos, dist_max, TLS_FLAG);
  ClearCFlag(MARK_FLAG); // reset the mark
  if (m_verbose) { std::cerr << "...cyftools tls - found " << CountCFlag(TLS_FLAG) << " cells near B-cells" << std::endl; }

  // check that we found something
  if (CountCFlag(TLS_FLAG) == 0) {
    std::cerr << "cyftools tls -- Warning - no b-cell clusters found, using B-cell flag " << bcell_marker << std::endl;
    return;
  }
  
  // Step 3) Clear the TLS flag if its not a B-cell
  //        , we are for now just finding B-cell clusters
  if (m_verbose) { std::cerr << "...cyftools tls - removing non-bcells from clusters" << std::endl; }
  for (size_t i = 0; i < n; i++) {
    cy_uint& cf  = (*m_cflag_ptr)[i];
    cy_uint pf  = m_pflag_ptr->at(i);
    
    // clear the TLS flag if not a B-cell
    if (!IS_FLAG_SET(pf, bcell_marker))
      CLEAR_FLAG(cf, TLS_FLAG);
  }
  if (m_verbose) { std::cerr << "...cyftools tls - found " << CountCFlag(TLS_FLAG) << " B-cells near B-cells" << std::endl; }

  // check that we found something
  if (CountCFlag(TLS_FLAG) == 0) {
    std::cerr << "cyftools tls -- Warning - clusters found after removing non-B-cells " << bcell_marker << std::endl;
    return;
  }

  
  // Step 4) Find cells that are near B-cell clusters (TLS_FLAG)
  //         This expands the area around the TLS
  num_neighbors = 100;
  frac_pos = 0.2;
  ClearCFlag(MARK_FLAG); // ensure mark is cleared
  CopyCFlag(TLS_FLAG, MARK_FLAG); // set mark as TLS for Annotate
  /// this will look for cells that are near MARK cells, and add TLS_FLAG
  // so after this, number of TLS_FLAG cells is >= number from before the call
  if (m_verbose) { std::cerr << "...cyftools tls - KNN to find cells near" << std::endl; }
  int dist_max_immune = 200;
  AnnotateCall(num_neighbors, frac_pos, dist_max_immune, TLS_FLAG);
  
  // Step 5) Remove TLS_FLAG from non-immune cells
  if (m_verbose) { std::cerr << "...cyftools tls - removing non-immune cells" << std::endl; }
  for (size_t i = 0; i < n; i++) {
    cy_uint& cf  = (*m_cflag_ptr)[i];
    cy_uint pf  = m_pflag_ptr->at(i);
    
    // clear the TLS flag if not an immune cell (any of flags in immune-Marker oK)
    if (!IS_FLAG_SET_OR(pf, immune_marker))
      CLEAR_FLAG(cf, TLS_FLAG);
  }
  if (m_verbose) { std::cerr << "...cyftools tls - found " << CountCFlag(TLS_FLAG) << " cells (nominal) in TLSes" << std::endl; }
  
  // Step 6) DBscan clustering to identify clusters of TLS
  ClearCFlag(MARK_FLAG); // ensure mark is cleared
  CopyCFlag(TLS_FLAG, MARK_FLAG); // dbscan will cluster now on TLS_FLAG
  if (m_verbose) { std::cerr << "...cyftools tls - dbscan clustering on TLS" << std::endl;  }
  // dbscan params
  float epsilon = 100.0f;
  int min_size = 25;
  clusterDBSCAN(epsilon, min_size, 200); //min_cluster_size);

  // Step 7) Clear the TLS flag and re-label TLS if part of a non-zero dbscan cluster
  ClearCFlag(MARK_FLAG + TLS_FLAG); // clear both flags
  
  // get the dbscan_cluster column
  std::unordered_map<std::string, FloatColPtr>::iterator it = m_table.find("dbscan_cluster");
  assert(it != m_table.end());
  FloatColPtr db_cl_ptr = it->second;
  
  // setup a new column for tls id
  FloatColPtr tls_id = std::make_shared<FloatCol>();
  tls_id->reserve(n);
  
  // loop cells and transfer cluster number to TLS
  if (m_verbose) { std::cerr << "...cyftools tls - transfering dbscan clusters to tls" << std::endl; }
  for (size_t i = 0; i < n; i++) {
    float dbcluster_num  = (*db_cl_ptr)[i];
    tls_id->PushElem(dbcluster_num);
  }
  assert(tls_id->size() == n);

  // Step 9) Convex hull to get cells inside
  //

  // get unique TLS ids
  std::set<float> unique_tls(tls_id->begin(), tls_id->end()); 

  // loop the clusters and make the convex hull
  if (m_verbose) { std::cerr << "...cyftools tls - looping clusters and making convex hull" << std::endl; } 
  for (const auto& cl : unique_tls) {
    
    // cluster 0 is holder for not a cluster
    if (cl == 0)
      continue;


    // fill polygon with the points in the TLS definition
    //    that will need to have hull around them
    std::vector<JPoint> polygon;
    polygon.reserve(n);
    for (size_t i = 0; i < n; i++) {
      if (tls_id->at(i) == cl)
	polygon.push_back(JPoint(m_x_ptr->at(i), m_y_ptr->at(i)));
    }

    // get the convex hull
    std::vector<JPoint> convex_hull = convexHull(polygon);

    if (polygon.empty())
      continue;
    
    // Make a polygon (roi) from the points
    Polygon tls_boundary_hull(convex_hull);

    // loop all cells and see if they are members, if so, add
    for (size_t i = 0; i < n; i++) {
      if (tls_boundary_hull.PointIn(m_x_ptr->at(i), m_y_ptr->at(i))) {
	if (IS_FLAG_SET_OR(m_pflag_ptr->at(i), immune_marker)) { // tls contains only immune cells
	  tls_id->SetValueAt(i, cl); // set the tls id
	}
      }
    }
    
  } // end cluster loop

  // Step 9
  // HIstogram the tls counts
  if (m_verbose) { std::cerr << "...cyftools tls - histograming the tls counts " << std::endl; } 
  std::unordered_map<float, size_t> histo;
  for (const auto& t : *tls_id) 
    histo[t]++;

  // set tls-id elems to 0 if below threshold
  for (auto& t : *tls_id)
    if (histo[t] < min_cluster_size)
      t = 0;
  
  // add the tls_id to the table
  Tag tls_tag(Tag::CA_TAG, "tls_id", "");
  AddColumn(tls_tag, tls_id);

  // Step 8) Turn the TLS flag back on to mark individual cells
  for (size_t i = 0; i < n; i++) {
    cy_uint& cf  = (*m_cflag_ptr)[i];
    if (tls_id->at(i) > 0)
      SET_FLAG(cf, TLS_FLAG);
  }

  
  return;
  
}  
  
void CellTable::AnnotateCall(int num_neighbors, float frac,
			     cy_uint dist, cy_uint flag_to_set) {

#ifdef HAVE_KNNCOLLE
  
  int nobs = CellCount();
  // vector to say if build mark flag is on or off for that cell
  // that is, only cells that are marked as BUILD_GRAPH_FLAG *ON* are added
  std::vector<bool> ix;
  ix.reserve(nobs);
  
  // not really exposed, so should never run for right now
  /*  if (build_tree_with_marked_only) {

    // not really exposed, so should never run for right now
    assert(false);
    
    // add true for cells with build flag marked
    size_t n = 0;
    for (size_t i = 0; i < m_cflag_ptr->size(); ++i) {
      ix.push_back(IS_FLAG_SET(m_cflag_ptr->at(i), BUILD_GRAPH_FLAG));
      if (ix.back())
	n++;
    }

    // warn and exit if no cells
    if (n == 0) {
      std::cerr << "Warning: AnnotateCall - no build flags set (from cyftools filter <> <> -B)" << std::endl;
      std::cerr << "         May need to either" << std::endl;
      std::cerr << "            1) run cyftools filter with -B set or change filter criteria" << std::endl;
      std::cerr << "            2) run cyftools annotate *without* the -B flag" << std::endl;
       
    }

    build_vp_tree(ix);
  }

  // build flag not set, so build all of them
  else {
  */    
  ix = std::vector<bool>(nobs, true);
  
  //}
  
  knncolle::VpTree<knncolle::distances::Euclidean, int, float> searcher = build_vp_tree(ix);

#ifdef HAVE_OMP  
#pragma omp parallel for num_threads(m_threads)
#else
  std::cerr << "OMP not included, no support for multithreading. Compile with to support" << std::endl;
#endif  
  for (size_t i = 0; i < nobs; ++i) {
    
    // verbose printing
    if (i % 500000 == 0 && m_verbose)
      std::cerr << "...annotating for cell " << setw(12) << AddCommas(i) << 
	" K " << setw(3) << num_neighbors << " Dist: " << setw(9) << dist <<
	" frac " << setw(4) << frac << " flag " << setw(4) << flag_to_set << 
	std::endl;
    
    JNeighbors neigh = searcher.find_nearest_neighbors(i, num_neighbors);

    // distance trim
    neigh.erase(
		std::remove_if(neigh.begin(), neigh.end(), 
			       [dist](const std::pair<int, float>& p) { return p.second > dist; }),
		neigh.end()
		);
    
    // finally do tumor stuff
    float tumor_cell_count = 0;
    for (size_t j = 0; j < neigh.size(); j++) {
      
      size_t indexr = neigh.at(j).first;
      
      //debug
      /*
      if (m_x_ptr->at(i) == 453 & m_y_ptr->at(i) == 5421)      
	std::cerr << " x " << m_x_ptr->at(indexr) << " y " <<
	  m_y_ptr->at(indexr) << " panck " << IS_FLAG_SET(node.m_flags.at(j), MARK_FLAG) <<
	  " cflag " << m_cflag_ptr->at(indexr) <<
	  " cflag2 " << cflag_vec[j] << 
	  " pflag " <<
	  m_pflag_ptr->at(indexr) <<
	  " id " << m_id_ptr->at(indexr) << std::endl;
      */
      
      if (IS_FLAG_SET(m_cflag_ptr->at(indexr), MARK_FLAG)) {
      	tumor_cell_count++;
      }
    }
    
    // debug
    /*if (m_x_ptr->at(i) == 453 & m_y_ptr->at(i) == 5421)
      std::cerr << " node size " << node.size() <<
	" frac " << frac << " tumor_cell_count " <<
	tumor_cell_count << std::endl;
    */
    
    if (tumor_cell_count / static_cast<float>(neigh.size()) >= frac)
      SET_FLAG((*m_cflag_ptr)[i], flag_to_set); 

  }// end for

  // clear the marks now that all done
  for (size_t i = 0; i < nobs; ++i) {
    CLEAR_FLAG((*m_cflag_ptr)[i], MARK_FLAG);
  }
  
#else
  std::cerr << "Warning: tumor call function requires including header library knncolle (https://github.com/LTLA/knncolle)" <<
    " and preprocessor directive -DHAVE_KNNCOLLE" << std::endl;
#endif
  
}

void CellTable::TumorMargin(float dist,
			    cy_uint tumor_flag,
			    cy_uint margin_flag) {

  validate();

  if (tumor_flag == 0) {
    std::cerr << "Error - cyftools annotate - tumor_flag is zero, needs to be non-zero" << std::endl;
    return;
  } 

  if (margin_flag == 0) {
    std::cerr << "Error - cyftools annotate - margin_flag is zero, needs to be non-zero" << std::endl;
    return;
  } 
  
  // build the tree
  if (m_verbose)
    std::cerr << "...building the KDTree" << std::endl;
  BuildKDTree();

  // loop the cells
  size_t countr = 0; // for progress reporting
  if (m_verbose)
    std::cerr << "...starting cell loop -- units = cells / 1000^2 pixels" << std::endl;

  const size_t num_cells = CellCount();
  
  // do a quick loop to make sure tumor is labeled
  size_t tumor_count = 0;
  for (size_t i = 0; i < num_cells; i++) {
    if (IS_FLAG_SET(m_cflag_ptr->at(i), tumor_flag))
      tumor_count++;
  }
  if (tumor_count == 0) {
    std::cerr << "*******************************************************************" << std::endl;
    std::cerr << "WARNING: No cells labeled as tumor. Likely need to run cysift tumor" << std::endl;
    std::cerr << "*******************************************************************" << std::endl;    
  }

  size_t margin_count = 0;
  
#ifdef HAVE_OMP  
#pragma omp parallel for num_threads(m_threads) schedule(dynamic, 100)
#else
  std::cerr << "OMP not included, no support for multithreading. Compile with to support" << std::endl;
#endif  
  for (size_t i = 0; i < num_cells; i++) {
    std::vector<float> cell_count;
    std::vector<size_t> total_cell_count;
    
    float x1 = m_x_ptr->at(i);
    float y1 = m_y_ptr->at(i);

    const float query[2] = { x1, y1 };

    std::vector<std::vector<size_t>> neighbors;
    std::vector<std::vector<double>> distances;
    cyf::KRange r(0.0, dist);

    // inclusive of this point
    m_kdtree->Search(query, r, neighbors, distances);

    // only one query point, so just reference that
    const std::vector<size_t>& ind = neighbors.at(0);

    const bool cell_is_tumor = IS_FLAG_SET(m_cflag_ptr->at(i), tumor_flag);

    // clear current margin flag
    CLEAR_FLAG((*m_cflag_ptr)[i], margin_flag); 
    
    // loop the nodes connected to each cell
    for (size_t n = 0; n < ind.size(); n++) {
      const uint32_t ncflag = m_cflag_ptr->at(ind.at(n)); // neighbor c flag
      if ( ( cell_is_tumor && ! IS_FLAG_SET(ncflag, tumor_flag)) || // cell is tumor, neighbor is stroma
	   (!cell_is_tumor &&   IS_FLAG_SET(ncflag, tumor_flag))) { // cell is stroma, neighbor is tumor
	SET_FLAG((*m_cflag_ptr)[i], margin_flag); // mark as margin
	assert(IS_FLAG_SET(m_cflag_ptr->at(i), margin_flag));
	margin_count++;
	break;
      }
    }
  }
  
  if (m_verbose)
    std::cerr << "...identified " << AddCommas(margin_count) << " cells at the margin within a distance of " << dist << std::endl;
  

}

void CellTable::BuildKDTree() {

  validate();

  if (m_verbose)
    std::cerr << "...building KD-tree (nanoflann)" << std::endl;

  // coordinate columns must outlive the tree (they persist on the CellTable)
  m_kdtree = std::make_unique<cyf::KDTree2D>(m_x_ptr->getData().data(),
                                             m_y_ptr->getData().data(),
                                             m_x_ptr->size());
}


int CellTable::RadialDensityKD(std::vector<cy_uint> inner, std::vector<cy_uint> outer,
			       std::vector<cy_uint> logor, std::vector<cy_uint> logand,
			       std::vector<std::string> label, std::vector<int> normalize_local,
			       std::vector<int> normalize_global) {

  // check the radial geometry parameters
  assert(inner.size());
  assert(inner.size() == outer.size());
  assert(inner.size() == logor.size());
  assert(inner.size() == logand.size());
  assert(inner.size() == label.size());
  assert(inner.size() == normalize_local.size());
  assert(inner.size() == normalize_global.size());  

  if (m_verbose)  {
    for (size_t i = 0; i < inner.size(); i++) {
      std::cerr << "...Radius: [" << std::setw(4) << inner.at(i) << "," << std::setw(4) << outer.at(i) <<
	"] OR: " << std::setw(7) << logor.at(i) << " AND: " << std::setw(7) << logand.at(i) << " - " << label.at(i) << std::endl;
    }
  }

  validate();
  
  // store the densities plus total cell count
  std::vector<std::shared_ptr<FloatCol>> dc(inner.size() + 1);
  for(auto& ptr : dc) {
    ptr = std::make_shared<FloatCol>();
    ptr->resize(m_pflag_ptr->size());
  }

  // build the tree
  if (m_verbose)
    std::cerr << "...building the KDTree" << std::endl;
  BuildKDTree();
  
  // get max radius to compute on
  cy_uint max_radius = 0;
  for (const auto& r : outer)
    if (max_radius < r)
      max_radius = r;

  // pre-compute the bools
  if (m_verbose)
    std::cerr << "...pre-computing flag results" << std::endl;
  std::vector<std::vector<int>> flag_result(inner.size(), std::vector<int>(m_pflag_ptr->size(), 0));
  for (size_t i = 0; i < m_pflag_ptr->size(); i++) {
    CellFlag mflag(m_pflag_ptr->getData().at(i));
    for (size_t j = 0; j < inner.size(); j++) {
      if (logor[j]  == std::numeric_limits<cy_uint>::max() ||
	  logand[j] == std::numeric_limits<cy_uint>::max()) { // special for ALL cells
	flag_result[j][i] = 1;
      } 
      else if ( (!logor[j] && !logand[j]) || mflag.testAndOr(logor[j], logand[j])) {
	flag_result[j][i] = 1; 
      }
    }
  }

  // get the global cell counts for each flag
  std::vector<size_t> flag_count(inner.size());
  assert(flag_result.size() == inner.size());
  for (size_t i = 0; i < flag_result.size(); i++) {
    flag_count[i] = std::accumulate(flag_result[i].begin(), flag_result[i].end(), 0);
  }

  // total cell count
  float total_image_cell_count = CellCount();
  
  // loop the cells
  size_t countr = 0; // for progress reporting
  if (m_verbose)
    std::cerr << "...starting cell loop -- units = cells / 1000^2 pixels" << std::endl;

  // inner and outer can be used in square mode only from now on
  //for (size_t i = 0; i < inner.size(); i++) {
  //  inner[i] = inner[i] * inner[i];
  // outer[i] = outer[i] * outer[i];
  //}
      
  // calculate the area
  std::vector<float> area(inner.size());
  for (size_t j = 0; j < inner.size(); ++j) {
    float outerArea = static_cast<float>(outer[j] * outer[j]) * 3.1415926535f;
    float innerArea = static_cast<float>(inner[j] * inner[j]) * 3.1415926535f;
    area[j] = outerArea - innerArea;
  }
  
  // display the header for the columns being display
  if (m_verbose) {
    std::cerr << std::string(47, ' ');
      for (size_t j = 5; j < dc.size(); j++) {
	std::cerr << std::setw(12) << label.at(j) << " ";
	if (j > 10)
	  break;
      }
      std::cerr << " ... " << std::endl;
    std::cerr << std::string(47, ' ');
      for (size_t j = 5; j < dc.size(); j++) {
	std::cerr << std::setw(12) << logor.at(j) << " ";
	if (j > 10)
	  break;
      }
      std::cerr << " ... " << std::endl;
      
  }

  //static std::vector<float> cell_count;
  //static std::vector<size_t> total_cell_count;  
  //    #pragma omp threadprivate(cell_count)
  //    #pragma omp threadprivate(total_cell_count)  

  // this will be inclusive of this point

  assert(m_kdtree);
#ifdef HAVE_OMP
#pragma omp parallel for num_threads(m_threads) schedule(dynamic, 100)
#else
  std::cerr << "OMP not included, no support for multithreading. Compile with to support" << std::endl;
#endif
  for (size_t i = 0; i < m_pflag_ptr->size(); i++) {
    std::vector<float> cell_count;
    std::vector<size_t> total_cell_count;  
    
    // Initialize the counts for each radial condition
    cell_count.resize(inner.size());
    //std::fill(cell_count.begin(), cell_count.end(), 0.0f);
    total_cell_count.resize(inner.size());
    //std::fill(total_cell_count.begin(), total_cell_count.end(), 0);

    float x1 = m_x_ptr->getData().at(i);
    float y1 = m_y_ptr->getData().at(i);

    const float query[2] = { x1, y1 };

    std::vector<std::vector<size_t>> neighbors;
    std::vector<std::vector<double>> distances;
    cyf::KRange r(0.0, max_radius);

    // inclusive of this point
    m_kdtree->Search(query, r, neighbors, distances);

    // only one query point, so just reference that
    const std::vector<size_t>& ind = neighbors.at(0);
    const std::vector<double>& dist = distances.at(0);

    // loop the nodes connected to each cell
    for (size_t n = 0; n < ind.size(); n++) {
      // test if the connected cell meets the flag criteria
      // n.first is cell_id of connected cell to this cell

      /*      if (i == 5000)
	std::cerr << " n " << n << " ind[n] " << ind[n] <<
	  " flag_result[9] " << flag_result[9][ind[n]] << " label " <<
	  label[9] << " flag[ind[n]] " << m_pflag_ptr->getData().at(ind[n]) << std::endl;
      */
      
      for (size_t j = 0; j < inner.size(); j++) {
	float d = dist[n]; 
	if (d >= inner[j] && d <= outer[j])
	  total_cell_count[j]++;
	if (flag_result[j][ind[n]])
	  cell_count[j] += ((d >= inner[j]) && (d <= outer[j]));
      }
    }
    
    // do the density calculation for each condition
    // remember, i is iterator over cells, j is over conditions
    for (size_t j = 0; j < area.size(); ++j) {
      float value = cell_count[j] * 1000000 / area[j]; // density per 1000 square um
      if (normalize_local[j] != 0) // normalize to cell count in the radius
	value = total_cell_count[j] == 0 ? 0 : value / total_cell_count[j];
      //if (normalize_global[j] != 0 && flag_count[j] != 0) // normalize to cell count in the image
      //	value = total_cell_count[j] == 0 ? 0 : value / flag_count[j] * total_image_cell_count;
      dc[j]->SetNumericElem(value, i);
    }

    const int print_interval = 100000;
    if (m_verbose && i % print_interval == 0) {
      countr +=print_interval;
      std::cerr << std::fixed << "...cell " << std::setw(11) << AddCommas(i) 
		<< " thr " << std::setw(2) << omp_get_thread_num() 
		<< " %done: " << std::setw(3) << static_cast<int>(static_cast<float>(countr) / m_pflag_ptr->size() * 100) 
		<< " count: ";
      
      for (size_t j = 5; j < dc.size(); j++) {
	std::cerr << std::setw(12) << static_cast<int>(std::round(dc.at(j)->getData().at(i))) << " ";
	if (j > 10)
	  break;
      }
      std::cerr << " ... " << std::endl;
    }
    
  } // end the main cell loop
  
  if (m_verbose)
    std::cerr << "...adding the density column" << std::endl;
  
  for (size_t i = 0; i < label.size(); ++i) {
    
    // form the data tag
    Tag dtag(Tag::CA_TAG, label[i], "");
    
    AddColumn(dtag, dc[i]);
    
  }

  return 0;
}

void CellTable::Select(CellSelector select, 
		       SelectOpMap criteria,
		       bool or_toggle,
		       float radius) {

  // add a dummy to ensure that OutputTable recognizes that m_cells_to_write is non-empty
  m_cells_to_write.insert(-1);
  validate();
  
  if (m_verbose)
    std::cerr << "...starting loop to tag cells as include/exclude" << std::endl;
  
  size_t count = CellCount();
  for (size_t i = 0; i < count; i++) {

    bool write_cell = false;
    
    ///////
    // FLAGS
    ///////
    // NB: even if flags are all empty, default should be to trigger a "write_cell = true"
    cy_uint pflag = m_pflag_ptr->getData().at(i);
    cy_uint cflag = m_cflag_ptr->getData().at(i);
    if (select.TestFlags(pflag, cflag)) {
      write_cell = true;
    }
    
    ///////
    // FIELD
    ///////
    bool flag_write_cell = write_cell;
    write_cell = or_toggle ? false : write_cell; // if or toggle is on, then ignore flag criteria and set true

    for (const auto& c : criteria) {

      // to-do -- will segfault if column not in table
      FloatColPtr ptr = m_table.at(c.first);
      
      float value = ptr->getData().at(i);
    
      for (const auto& cc : c.second) { // loop the vector of {optype, float}
	
	switch (cc.first) { // switching on the optype
	case optype::GREATER_THAN: write_cell          = write_cell = or_toggle ? (write_cell || (value >  cc.second)) : (write_cell && (value >  cc.second));  break;
	case optype::LESS_THAN: write_cell             = write_cell = or_toggle ? (write_cell || (value <  cc.second)) : (write_cell && (value <  cc.second));  break;
	case optype::GREATER_THAN_OR_EQUAL: write_cell = write_cell = or_toggle ? (write_cell || (value >= cc.second)) : (write_cell && (value >= cc.second));  break;
	case optype::LESS_THAN_OR_EQUAL: write_cell    = write_cell = or_toggle ? (write_cell || (value <= cc.second)) : (write_cell && (value <= cc.second));  break;
	case optype::EQUAL_TO: write_cell              = write_cell = or_toggle ? (write_cell || (value == cc.second)) : (write_cell && (value == cc.second));  break;
	default: assert(false);
	}
      }
    } // end criteria loop

    if (flag_write_cell && write_cell)
      m_cells_to_write.insert(m_id_ptr->getData().at(i));
    
  } // end cell loop

  // second loop to actually write
  
  // build the tree
  if (m_verbose)
    std::cerr << "...building the KDTree" << std::endl;
  BuildKDTree();
  
  if (m_verbose)
    std::cerr << "...starting second loop to include cells within radius of include" << std::endl;

#ifdef HAVE_OMP  
#pragma omp parallel for num_threads(m_threads)
#else
  std::cerr << "OMP not included, no support for multithreading. Compile with to support" << std::endl;
#endif  
  for (size_t i = 0; i < count; i++) {

    if (i % 100000 == 0 && m_verbose)
      std::cerr << "...working on cell " << AddCommas(i) << std::endl;
    
    // if this cell is write, no KD lookup needed
    int this_cid = m_id_ptr->getData().at(i);
    if (m_cells_to_write.count(this_cid)) {
      continue;
    }
    
    // find the point
    float x1 = m_x_ptr->getData().at(i);
    float y1 = m_y_ptr->getData().at(i);
    //point_t pt = { x1, y1 };

    // this will be inclusive of this point
    assert(m_kdtree);
    const float query[2] = { x1, y1 };

    std::vector<std::vector<size_t>> neighbors;
    std::vector<std::vector<double>> distances;
    cyf::KRange r(0.0, radius);

    m_kdtree->Search(query, r, neighbors, distances);
    
    // loop the inds, and check if neighbors a good cell
    for (const auto& j : neighbors.at(0)) {
      int neigh_cid = m_id_ptr->getData().at(j);
      if (m_cells_to_write.count(neigh_cid)) {
	m_cells_to_write.insert(this_cid);
	continue;
      }
    }
  }
  return;
}
