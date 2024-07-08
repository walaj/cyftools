#include "cell_table.h"
#include "color_map.h"
#include "cell_selector.h"
#include <stack>

#ifdef HAVE_UMAPPP
#include "umappp/Umap.hpp"
#include "umappp/NeighborList.hpp"
#include "umappp/combine_neighbor_sets.hpp"
#include "umappp/find_ab.hpp"
#include "umappp/neighbor_similarities.hpp"
#include "umappp/optimize_layout.hpp"
#include "umappp/spectral_init.hpp"
#endif

#ifdef HAVE_KMEANS
#include "kmeans/Kmeans.hpp"
#endif

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

void CellTable::UMAP(int num_neighbors) { 

#ifdef HAVE_UMAPPP
  
  // get the number of markers
  int ndim = 0;
  for (const auto& t : m_header.GetDataTags())
    if (t.type == Tag::MA_TAG)
      ndim++;

  // number of cells
  int nobs = CellCount();
  
  if (m_verbose)
    std::cerr << "...finding K nearest-neighbors (marker) on " << AddCommas(nobs) << " cells" << std::endl;  
  
  // column major the marker data
  std::vector<float> concatenated_data;

  // concatenate the data
  for (const auto& t : m_header.GetDataTags()) {

    // skip if not a marker tag
    if (t.type != Tag::MA_TAG)
      continue;
    
    auto it = m_table.find(t.id);
    if (it != m_table.end()) {
      FloatColPtr numeric_column_ptr = it->second;
      assert(numeric_column_ptr);
      
      // Check if the column is a NumericColumn
      const auto& data = it->second->getData();
      concatenated_data.insert(concatenated_data.end(), data.begin(), data.end());
    }
  }
  
  // convert to row major?
  column_to_row_major(concatenated_data, nobs, ndim);

  if (m_verbose)
    std::cerr << "...setting up KNN graph on " << AddCommas(nobs) << " points (" <<
      ndim << "-dimensional)" << std::endl;
  
  if (m_verbose)
    std::cerr << "...building KNN (marker) graph" << std::endl;
  
  // initialize the tree. Can choose from different algorithms, per knncolle library
  knncolle::VpTree<knncolle::distances::Euclidean, int, float> searcher(ndim, nobs, concatenated_data.data());
  //knncolle::Annoy<knncolle::distances::Euclidean, int, float> searcher(ndim, nobs, concatenated_data.data());
  //knncolle::Kmknn<knncolle::distances::Euclidean, int, float> searcher(ndim, nobs, concatenated_data.data());

  umappp::NeighborList<float> nlist;
  nlist.resize(nobs);

#ifdef HAVE_OMP  
#pragma omp parallel for num_threads(m_threads)
#else
  std::cerr << "OMP not included, no support for multithreading. Compile with to support" << std::endl;
#endif  
  for (size_t i = 0; i < nobs; ++i) {
    if (i % 50000 == 0 && m_verbose)
      std::cerr << "...processing cell "
		<< std::setw(12) << std::left << AddCommas(i) // Set width to 12 and left justify
		<< " on thread "
		<< std::setw(2) << std::left << omp_get_thread_num() // Set width to 2 and left justify
		<< " K "
		<< std::setw(3) << std::left << num_neighbors // Set width to 3 and left justify
		<< std::endl;      
    
    JNeighbors neigh = searcher.find_nearest_neighbors(i, num_neighbors);

#ifdef HAVE_OMP  
#pragma omp critical
#endif  
    {
      nlist[i] = neigh;
    }
  }

  //debug
  // size_t i = 0;
  // for (const auto& nn : nlist) {
  //   std::cerr << " i " << i << std::endl;
  //   i++;
  //   for (const auto& mm: nn) {
  //     std::cerr << " neighbor " << mm.first << " - Distance " << mm.second << std::endl;
  //   }
  //   std::cerr << std::endl;
  // }

  // run umap
  umappp::Umap<float> umapr;
  umapr.set_initialize(umappp::InitMethod::RANDOM);
  if (m_threads > 1)
    umapr.set_parallel_optimization(true);
  umapr.set_num_threads(m_threads);
  umapr.set_seed(42);
  std::vector<float> embedding(nobs*2);
  
  if (m_verbose)
    std::cerr << "...initializing UMAP " << std::endl;
  
  auto status = umapr.initialize(std::move(nlist), 2, embedding.data());
  
  if (m_verbose)
    std::cerr << "...running UMAP" << std::endl;
  status.run(); 

  if (m_verbose)
    std::cerr << "...done with umap" << std::endl;

  // rescale to -1 to 1
  // find abs value
  float max_abs = 0.0f;
  for (int i = 0; i < nobs*2; i++) {
    max_abs = std::max(max_abs, std::abs(embedding[i]));
  }
  assert(max_abs > 0);
    
  // transfer to Column
  FloatColPtr comp1 = make_shared<FloatCol>(); 
  FloatColPtr comp2 = make_shared<FloatCol>();

  comp1->resize(nobs); 
  comp2->resize(nobs); 

  for (size_t i = 0; i < nobs; i++) {
    comp1->SetNumericElem(embedding[i*2    ] / max_abs, i);
    comp2->SetNumericElem(embedding[i*2 + 1] / max_abs, i);
  }

  //
  Tag gtag1(Tag::CA_TAG, "umap1", "NN:" + std::to_string(num_neighbors));  
  AddColumn(gtag1, comp1);

  Tag gtag2(Tag::CA_TAG, "umap2", "NN:" + std::to_string(num_neighbors));  
  AddColumn(gtag2, comp2);

#else
  std::cerr << "Warning: Need to include umappp header library (https://github.com/LTLA/umappp) and adding -DHAVE_UMAPPP preprocessor directive" << std::endl;
#endif
  
}

void CellTable::UMAPPlot(const std::string& file, int width, int height) const {

  if (m_table.find("umap1") == m_table.end() ||
      m_table.find("umap1") == m_table.end()) {
    std::cerr << "Warning: Need to build UMAP first before plotting" << std::endl;
    return;
  }

  ColorMap colormap;
  
#ifdef HAVE_CAIRO
  FloatColPtr u1_ptr = m_table.at("umap1");
  FloatColPtr u2_ptr = m_table.at("umap2");
  assert(u1_ptr->size() == u2_ptr->size());

  // clusters from kmeans
  std::vector<int> clusters(u1_ptr->size());

 
#ifdef HAVE_KMEANS

  // store as concatenated data
  std::vector<double> concatenated_umap_data;
  const auto& u1d = u1_ptr->getData();
  concatenated_umap_data.insert(concatenated_umap_data.end(), u1d.begin(), u1d.end());
  const auto& u2d = u2_ptr->getData();
  concatenated_umap_data.insert(concatenated_umap_data.end(), u2d.begin(), u2d.end());

  const int ndim = 2;
  const int nobs = u1_ptr->size();
  const int nclusters = 6;

  if (nclusters > 12) {
    std::cerr << "Requesting more clusters than available colors. Max 12. Will recycle colors" << std::endl;
  }

  kmeans::InitializeRandom rd;
  kmeans::HartiganWong hw;
  hw.set_max_iterations(100);
  hw.set_num_threads(m_threads);
  kmeans::Kmeans km;
  km.set_seed(42);
  km.set_num_threads(m_threads);
  
  // auto res = km.run(ndim, nobs, concatenated_umap_data.data(), nclusters, &rd, &hw);
  
  // kmeans::Kmeans km;
  // km.set_seed(42);
  // km.set_num_threads(m_threads);

  if (m_verbose)
    std::cerr << "...running k-means clustering with " << nclusters << " clusters" << std::endl;
  auto res = kmeans::Kmeans().run(ndim, nobs, concatenated_umap_data.data(), nclusters);
  
  clusters = res.clusters;
  
  // set the best colormap based on number of clusters
  colormap = getColorMap(nclusters);
  
  if (m_verbose) {
    std::cerr << "...done with kmeans" << std::endl;
    std::cerr << "   Iterations: " << res.details.iterations << std::endl;
    std::cerr << "   Status: " << res.details.status << std::endl;
    for (size_t i = 0; i < nclusters; i++)
      std::cerr << "   Cluster " << i << " -- size " << res.details.sizes.at(i) << std::endl;
  }
  
#elif defined(HAVE_MLPACK)
  
  
#else
  std::cerr << "Warning: Unable to run umap *clustering* without including CppKmeans (https://github.com/LTLA/CppKmeans)" <<
    " and -DHAVE_KMEANS preprocessor directive" << std::endl;
#endif
  
  // open the PDF for drawing
  const int cwidth = 500;
  cairo_surface_t *surface = cairo_pdf_surface_create(file.c_str(), cwidth, cwidth);
  cairo_t *cr = cairo_create(surface);

  for (size_t i = 0; i < u1_ptr->size(); i++) {
    Color c = colormap.at(clusters.at(i) % 12);
    cairo_set_source_rgba(cr, c.redf(), c.greenf(), c.bluef(), 0.3);
    cairo_arc(cr, u1_ptr->at(i)*cwidth/2+cwidth/2, u2_ptr->at(i)*cwidth/2+cwidth/2, 0.5, 0.0, 2.0 * M_PI);
    cairo_fill(cr);
  }

  // Clean up and close
  cairo_destroy(cr);
  cairo_surface_destroy(surface);
    
  
#else
  std::cerr << "Error: Need to build / link with Cairo library to build PDF." << std::endl;
#endif
  
}

/*void CellTable::KNN_spatial(int num_neighbors, int dist) {

#ifdef HAVE_KNNCOLLE
  
  // number of cells
  int nobs = CellCount();
  
  if (m_verbose)
    std::cerr << "...finding K nearest-neighbors (spatial) on " << AddCommas(nobs) << " cells" << std::endl;  
  
  // column major the coordinate data
  std::vector<float> concatenated_data;
  
  const int ndim = 2;

  assert(m_table.find("x") != m_table.end());
  assert(m_table.find("y") != m_table.end());
  const FloatColPtr x_ptr = dynamic_pointer_cast<FloatCol>(m_table.at("x"));
  const FloatColPtr y_ptr = dynamic_pointer_cast<FloatCol>(m_table.at("y"));
  
  const auto& x_data = x_ptr->getData(); // just a ptr, not a copy
  if (m_verbose) 
    std::cerr << "...adding " << AddCommas(x_ptr->size()) << " points on x" << std::endl;
  concatenated_data.insert(concatenated_data.end(), x_data.begin(), x_data.end());
  
  const auto& y_data = y_ptr->getData();
  if (m_verbose) 
    std::cerr << "...adding " << AddCommas(y_ptr->size()) << " points on y" << std::endl;
   concatenated_data.insert(concatenated_data.end(), y_data.begin(), y_data.end());
   
   // convert to row major?
   column_to_row_major(concatenated_data, nobs, ndim);
   
   // get the cell id colums
   assert(m_table.find("id") != m_table.end());
   const IDColPtr id_ptr = dynamic_pointer_cast<IDCol>(m_table.at("id"));
   
   if (m_verbose)
     std::cerr << "...setting up KNN graph (spatial) on " << AddCommas(nobs) << " points" << std::endl;
   
  umappp::NeighborList<float> output(nobs);

  if (m_verbose)
    std::cerr << "...building KNN (spatial) graph" << std::endl;

  // store the final graph
  auto graph = std::make_shared<GraphColumn>();
  graph->resize(nobs);
  
  
  // track number of cases where all N nearest neighbors are within the
  // distance cutoff, implying that there are likely additional cells within
  // the desired radius that are cutoff
  size_t lost_cell = 0;
  
  // initialize the tree. Can choose from different algorithms, per knncolle library
  knncolle::VpTree<knncolle::distances::Euclidean, int, float> searcher(ndim, nobs, concatenated_data.data());
  //knncolle::AnnoyEuclidean<int, float> searcher(ndim, nobs, concatenated_data.data());
  //knncolle::Kmknn<knncolle::distances::Euclidean, int, float> searcher(ndim, nobs, concatenated_data.data());  
    
  // archive the header
  assert(m_archive);
  (*m_archive)(m_header);

  // create the cells and print
  size_t numRows = CellCount();

  // setup for converting to Cell
  assert(m_table.find("cflag") != m_table.end());
  assert(m_table.find("pflag") != m_table.end());
  const IntColPtr cflag_ptr = dynamic_pointer_cast<IntCol>(m_table.at("cflag"));
  const IntColPtr pflag_ptr = dynamic_pointer_cast<IntCol>(m_table.at("pflag"));  
  std::vector<FloatColPtr> col_ptr;
  for (const auto& t : m_header.GetDataTags()) {
    assert(m_table.find(t.id) != m_table.end());
    col_ptr.push_back(dynamic_pointer_cast<FloatCol>(m_table.at(t.id)));
  }
  
#pragma omp parallel for num_threads(m_threads)
  for (size_t i = 0; i < nobs; ++i) {

    // verbose printing

    if (i % 50000 == 0 && m_verbose)
      std::cerr << "...working on cell " <<
	AddCommas(i) << " with thread " <<
	omp_get_thread_num() << " K " << num_neighbors << " Dist: " << dist <<
	std::endl;
    
    JNeighbors neigh = searcher.find_nearest_neighbors(i, num_neighbors);


    CellNode node;

    // remove less than distance
    if (dist > 0) {
      
      size_t osize = neigh.size();
      JNeighbors neigh_trim;
      for (const auto& nnn : neigh) {
	if (nnn.second < dist)
	  neigh_trim.push_back(nnn);
      }
      
      // print a warning if we trimmed off too many neighbors
      if (osize == neigh_trim.size()) {  
#pragma omp critical
	{
	  lost_cell++;
	  if (lost_cell % 500 == 0)
	    std::cerr << "osize " << osize << " Lost cell " << AddCommas(lost_cell) << " of " << AddCommas(nobs) << std::endl;
	}
      }
      
      neigh = neigh_trim;

    }

    // add the pheno flags
    std::vector<cy_uint> pflag_vec(neigh.size());
    for (size_t j = 0; j < neigh.size(); j++) {
      pflag_vec[j] = static_cast<IntCol*>(pflag_ptr.get())->getData().at(neigh.at(j).first);
    }

    assert(pflag_vec.size() == neigh.size());
    
    // Now, set the node pointer to CellID rather than 0-based
    for (auto& cc : neigh) {
      cc.first = id_ptr->getData().at(cc.first);
    }
    node = CellNode(neigh, pflag_vec);

    // build the Cell object
    ////////////////////////
    Cell cell;
    cell.id          = id_ptr->getData().at(i);
    cell.pflag  = pflag_ptr->getData().at(i);
    cell.cflag   = cflag_ptr->getData().at(i);
    cell.x           = x_ptr->getData().at(i);
    cell.y           = y_ptr->getData().at(i);    

    // fill the Cell data columns
    for (const auto& c : col_ptr) {
      cell.cols.push_back(c->getData().at(i));
    }

    // fill the Cell graph data
    //node.FillSparseFormat(cell.m_spatial_ids, cell.m_spatial_dist);
    //assert(pflag_vec.size() == cell.m_spatial_ids.size());
    //cell.m_spatial_flags = pflag_vec;
    
    
    // Add the cell to the thread-local buffer.
#ifdef __clang__
    cell_buffer.push_back(cell);
#endif

    //if (cell_buffer.size() >= CELL_BUFFER_LIMIT) {
#pragma omp critical
    {
#ifdef __clang__
      for (const auto& buffered_cell : cell_buffer) {
      	(*m_archive)(buffered_cell);
      }
#else
      (*m_archive)(cell);      
#endif
    }
    
#ifdef __clang__    
    cell_buffer.clear();
#endif
    
  }// end for
  
  // After the loop, dump any remaining cells in the buffer.


#ifdef __clang__
#pragma omp critical
  {
    for (const auto& buffered_cell : cell_buffer) {
      (*m_archive)(buffered_cell);
    }
  }
  cell_buffer.clear();
#endif
  
  if (m_verbose)
    std::cerr << "...done with graph construction" << std::endl;
  
  return;
  
  Tag gtag(Tag::GA_TAG, "spat", "NN:" + std::to_string(num_neighbors));
  //gtag.addValue("NN", std::to_string(num_neighbors));
  AddColumn(gtag, graph);

#else
  std::cerr << "Warning: KNN spatial function requires including header library knncolle (https://github.com/LTLA/knncolle)" <<
    " and preprocessor directive -DHAVE_KNNCOLLE" << std::endl;
#endif

}*/

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
    std::cerr << "...cyftools dist - Finding distance to " << marked_cells << " marked cells" << std::endl;
   
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
  clusterDBSCAN(epsilon, min_size, 50); //min_cluster_size);

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
    if (i % 50000 == 0 && m_verbose)
      std::cerr << "...annoting for cell " << setw(12) << AddCommas(i) << 
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

#ifdef HAVE_KNNCOLLE

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

    arma::mat query(2, 1);
    query(0, 0) = x1;
    query(1, 0) = y1;
    
    std::vector<std::vector<size_t>> neighbors;
    std::vector<std::vector<double>> distances;
    mlpack::Range r(0.0, dist);

    // this will be inclusive of this point
    ml_kdtree->Search(query, r, neighbors, distances);

    // only one query point, so just reference that
    const std::vector<size_t>& ind = neighbors.at(0);
    const std::vector<double>& dist = distances.at(0);

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
  
#else
  std::cerr << "Warning: tumor call function requires including header library knncolle (https://github.com/LTLA/knncolle)" <<
    " and preprocessor directive -DHAVE_KNNCOLLE" << std::endl;
#endif
  
  
}

void CellTable::BuildKDTree() {

  validate();
  
#ifdef HAVE_MLPACK

  if (m_verbose)
    std::cerr << "...building ML pack KDTree" << std::endl;
  
  arma::mat data(2, m_x_ptr->size());
  for (size_t i = 0; i < m_x_ptr->size(); ++i)
    {
      data(0, i) = m_x_ptr->at(i);
      data(1, i) = m_y_ptr->at(i);
    }
  
  //assert(ml_kdtree == nullptr);
  ml_kdtree = new mlpack::RangeSearch<mlpack::EuclideanDistance>(data);
  assert(ml_kdtree);
  
#else
    std::cerr << "Warning: Unable to build KD-tree, need to include MLPack library " << 
    " and add preprocessor directive -DHAVE_MLPACK" << std::endl;
#endif
  
}

void CellTable::MoranI(const std::vector<cy_uint>& flags) {

  /*
  validate();
  
  // build the tree
  if (m_verbose)
    std::cerr << "...building the KDTree" << std::endl;
  BuildKDTree();


  // pre-compute the bools
  if (m_verbose)
    std::cerr << "...pre-computing flag results" << std::endl;
  std::vector<std::vector<int>> flag_result(inner.size(), std::vector<int>(m_pflag_ptr->size(), 0));
  for (size_t i = 0; i < m_pflag_ptr->size(); i++) {
    CellFlag mflag(m_pflag_ptr->getData().at(i));
    for (size_t j = 0; j < inner.size(); j++) {
      if (mflag.testAndOr(flags[j], 0)) {
	flag_result[j][i] = 1; 
      }
    }
  }


  std::vector<float> flag_moran(flags.size(), 0.0f);
  for (size_t i = 0; i < flag_count.size(); i++) {
    for (size_t x = 0; x < m_x_ptr->size(); x++) {
      CellFlag mflag(m_pflag_ptr->getData().at(x));
      if (mflag.testAndOr(flags.at(i))) {
	for (size_t y = 0; y < m_y_ptr->size(); y++) {
	  if (x != y && )
	    }
      }
    }
  }

  size_t countr = 0;
#pragma omp parallel for num_threads(m_threads) schedule(dynamic, 100)
  for (size_t i = 0; i < m_pflag_ptr->size(); i++) {
    
    float x1 = m_x_ptr->getData().at(i);
    float y1 = m_y_ptr->getData().at(i);

    arma::mat query(2, 1);
    query(0, 0) = x1;
    query(1, 0) = y1;
    
    std::vector<std::vector<size_t>> neighbors;
    std::vector<std::vector<double>> distances;
    mlpack::Range r(0.0, max_radius);

    // this will be inclusive of this point
    ml_kdtree->Search(query, r, neighbors, distances);

    // only one query point, so just reference that
    const std::vector<size_t>& ind = neighbors.at(0);
    const std::vector<double>& dist = distances.at(0);

    // vec to store morans I
    std::vector<float> flag_variance(flags.size(), 0.0f);    
    std::vector<int>   flag_count(flags.size(), 0);
    std::vector<float> flag_mean(flags.size(), 0.0f);

    // calculate the total number of cells for each flag type
    for (size_t j = 0; j < flags.size(); j++) {
      for (const auto& n : ind) {
	if (flag_result[j][n])
	  flag_count[j]++;
      }
    }
    
    // find the mean number of cells
    flag_mean[j] = static_cast<float>(flag_count.at(j)) / static_cast<float>(ind.size());
    
    // find the variance
    flag_variance[j] = (1 - flag_mean[j]) * (1 - flag_mean[j]) * flag_count[j];

    if (m_verbose && i % 5000 == 0) {
      countr += 5000;
      std::cerr << std::fixed << "...cell " << std::setw(11) << AddCommas(i) 
		<< " thr " << std::setw(2) << omp_get_thread_num() 
		<< " %done: " << std::setw(3) << static_cast<int>(static_cast<float>(countr) / m_pflag_ptr->size() * 100) 
		<< " moran: ";
      
      for (size_t j = 5; j < dc.size(); j++) {
	std::cerr << std::setw(12) << static_cast<int>(std::round(dc.at(j)->getData().at(i))) << " ";
	if (j > 10)
	  break;
      }
      std::cerr << " ... " << std::endl;
    }
    
  } // end the main cell loop
  */  
  
}

int CellTable::RadialDensityKD(std::vector<cy_uint> inner, std::vector<cy_uint> outer,
			       std::vector<cy_uint> logor, std::vector<cy_uint> logand,
			       std::vector<std::string> label, std::vector<int> normalize_local,
			       std::vector<int> normalize_global) {
#ifdef HAVE_MLPACK
  
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
      std::cerr << "Radius: [" << std::setw(4) << inner.at(i) << "," << std::setw(4) << outer.at(i) <<
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

  assert(ml_kdtree);
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

    arma::mat query(2, 1);
    query(0, 0) = x1;
    query(1, 0) = y1;
    
    std::vector<std::vector<size_t>> neighbors;
    std::vector<std::vector<double>> distances;
    mlpack::Range r(0.0, max_radius);

    // this will be inclusive of this point
    ml_kdtree->Search(query, r, neighbors, distances);

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
  
    if (m_verbose && i % 5000 == 0) {
      countr += 5000;
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

#else
  std::cerr << "Warning: Unable to run without including kdtree header library (https://github.com/crvs/KDTree)" << std::endl;
#endif
  
  return 0;
}

void CellTable::Select(CellSelector select, 
		       SelectOpMap criteria,
		       bool or_toggle,
		       float radius) {

#ifdef HAVE_MLPACK

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
    assert(ml_kdtree);
    arma::mat query(2, 1);
    query(0, 0) = x1;
    query(1, 0) = y1;

    std::vector<std::vector<size_t>> neighbors;
    std::vector<std::vector<double>> distances;
    mlpack::Range r(0.0, radius);

    // this will be inclusive of this point
    ml_kdtree->Search(query, r, neighbors, distances);
    
    // loop the inds, and check if neighbors a good cell
    for (const auto& j : neighbors.at(0)) {
      int neigh_cid = m_id_ptr->getData().at(j);
      if (m_cells_to_write.count(neigh_cid)) {
	m_cells_to_write.insert(this_cid);
	continue;
      }
    }
  }
#else
  std::cerr << "Warning -- not able to do radial selection without compiling with mlpack" << std::endl;
#endif
  
  return;
}
