#include "cell_table.h"
#include "color_map.h"

#ifdef HAVE_UMAPPP
#include "umappp/Umap.hpp"
#include "umappp/NeighborList.hpp"
#include "umappp/combine_neighbor_sets.hpp"
#include "umappp/find_ab.hpp"
#include "umappp/neighbor_similarities.hpp"
#include "umappp/optimize_layout.hpp"
#include "umappp/spectral_init.hpp"
#include "knncolle/knncolle.hpp"
#endif

#ifdef HAVE_KMEANS
#include "kmeans/Kmeans.hpp"
#endif

#ifdef HAVE_CAIRO
#include "cairo/cairo.h"
#include "cairo/cairo-pdf.h"
#endif

#ifdef HAVE_KDTREE
#include "KDTree.hpp"
#endif

// block for attempting bufferered writing to improve paralellization
// not really doing much, CELL_BUFFER_LIMIT 1 basically turns this off,
// but still here in case useful later
#define CELL_BUFFER_LIMIT 1

#ifdef __clang__
std::vector<Cell> cell_buffer;
#pragma omp threadprivate(cell_buffer)
#endif

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
      auto column_ptr = it->second;
      auto numeric_column_ptr = std::dynamic_pointer_cast<FloatCol>(column_ptr);
      assert(numeric_column_ptr);
      
      // Check if the column is a NumericColumn
      const auto& data = numeric_column_ptr->getData();
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
  
#pragma omp parallel for num_threads(m_threads)
  for (size_t i = 0; i < nobs; ++i) {
    if (i % 50000 == 0 && m_verbose)
      std::cerr << "...working on cell " <<
	AddCommas(i) << " with thread " <<
	omp_get_thread_num() << " K " << num_neighbors << std::endl;
    
    JNeighbors neigh = searcher.find_nearest_neighbors(i, num_neighbors);
    
#pragma omp critical
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
  std::shared_ptr<FloatCol> comp1 = make_shared<FloatCol>(); 
  std::shared_ptr<FloatCol> comp2 = make_shared<FloatCol>();

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
  FloatColPtr u1_ptr = dynamic_pointer_cast<FloatCol>(m_table.at("umap1"));
  FloatColPtr u2_ptr = dynamic_pointer_cast<FloatCol>(m_table.at("umap2"));  
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
    cairo_arc(cr, u1_ptr->GetNumericElem(i)*cwidth/2+cwidth/2, u2_ptr->GetNumericElem(i)*cwidth/2+cwidth/2, 0.5, 0.0, 2.0 * M_PI);
    cairo_fill(cr);
  }

  // Clean up and close
  cairo_destroy(cr);
  cairo_surface_destroy(surface);
    
  
#else
  std::cerr << "Error: Need to build / link with Cairo library to build PDF." << std::endl;
#endif
  
}

void CellTable::KNN_spatial(int num_neighbors, int dist) {

#ifdef HAVE_KNNCOLLE
  
  // number of cells
  int nobs = CellCount();
  
  if (m_verbose)
    std::cerr << "...finding K nearest-neighbors (spatial) on " << AddCommas(nobs) << " cells" << std::endl;  
  
  // column major the coordinate data
  std::vector<float> concatenated_data;
  
  const int ndim = 2;
  
  auto x_ptr = m_table.at("x");
  auto y_ptr = m_table.at("y");
  
  const auto& x_data = std::dynamic_pointer_cast<FloatCol>(x_ptr)->getData(); // just a ptr, not a copy
  if (m_verbose) 
    std::cerr << "...adding " << AddCommas(x_data.size()) << " points on x" << std::endl;
  concatenated_data.insert(concatenated_data.end(), x_data.begin(), x_data.end());
  
  const auto& y_data = std::dynamic_pointer_cast<FloatCol>(y_ptr)->getData();
  if (m_verbose) 
    std::cerr << "...adding " << AddCommas(y_data.size()) << " points on y" << std::endl;
   concatenated_data.insert(concatenated_data.end(), y_data.begin(), y_data.end());
   
   // convert to row major?
   column_to_row_major(concatenated_data, nobs, ndim);
   
   // get the cell id colums
   auto id_ptr = m_table.at("id");
   
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
  auto cflag_ptr = m_table.at("cflag");
  auto pflag_ptr = m_table.at("pflag");  
  std::vector<ColPtr> col_ptr;
  for (const auto& t : m_header.GetDataTags()) {
    col_ptr.push_back(m_table.at(t.id));
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
      pflag_vec[j] = static_cast<IntCol*>(pflag_ptr.get())->GetNumericElem(neigh.at(j).first);
    }

    assert(pflag_vec.size() == neigh.size());
    
    // Now, set the node pointer to CellID rather than 0-based
    for (auto& cc : neigh) {
      cc.first = id_ptr->GetNumericElem(cc.first);
    }
    node = CellNode(neigh, pflag_vec);

    // build the Cell object
    ////////////////////////
    Cell cell;
    cell.m_id   = static_cast<IntCol*>(id_ptr.get())->GetNumericElem(i);
    cell.m_pheno_flag = static_cast<IntCol*>(pflag_ptr.get())->GetNumericElem(i);
    cell.m_cell_flag = static_cast<IntCol*>(cflag_ptr.get())->GetNumericElem(i);    
    cell.m_x    = static_cast<FloatCol*>(x_ptr.get())->GetNumericElem(i);
    cell.m_y    = static_cast<FloatCol*>(y_ptr.get())->GetNumericElem(i);

    // fill the Cell data columns
    for (const auto& c : col_ptr) {
      cell.m_cols.push_back(static_cast<FloatCol*>(c.get())->GetNumericElem(i));
    }

    // fill the Cell graph data
    node.FillSparseFormat(cell.m_spatial_ids, cell.m_spatial_dist);
    assert(pflag_vec.size() == cell.m_spatial_ids.size());
    cell.m_spatial_flags = pflag_vec;

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
  
}

void CellTable::TumorCall(int num_neighbors, float frac,
			  cy_uint orflag, cy_uint andflag, cy_uint dist) {

#ifdef HAVE_KNNCOLLE
  
  // number of cells
  int nobs = CellCount();
  
  if (m_verbose)
    std::cerr << "...finding K nearest-neighbors (spatial) on " << AddCommas(nobs) << " cells" << std::endl;
  
  // column major the coordinate data
  std::vector<float> concatenated_data;
  
  const int ndim = 2;
  
  auto x_ptr = m_table.at("x");
  auto y_ptr = m_table.at("y");
  auto pflag_ptr = m_table.at("pflag");
  auto cflag_ptr = m_table.at("cflag");
  
  const auto& x_data = std::dynamic_pointer_cast<FloatCol>(x_ptr)->getData(); // just a ptr, not a copy
  if (m_verbose) 
    std::cerr << "...adding " << AddCommas(x_data.size()) << " points on x" << std::endl;
  concatenated_data.insert(concatenated_data.end(), x_data.begin(), x_data.end());
  
  const auto& y_data = std::dynamic_pointer_cast<FloatCol>(y_ptr)->getData();
  if (m_verbose) 
    std::cerr << "...adding " << AddCommas(y_data.size()) << " points on y" << std::endl;
   concatenated_data.insert(concatenated_data.end(), y_data.begin(), y_data.end());
   
   // convert to row major?
   column_to_row_major(concatenated_data, nobs, ndim);
   
   if (m_verbose)
     std::cerr << "...setting up KNN graph (spatial) for tumor calling on " << AddCommas(nobs) << " points" << std::endl;

   JNeighbors output(nobs);

  if (m_verbose)
    std::cerr << "...building KNN (spatial) graph with " <<
      num_neighbors << " nearest neigbors and dist limit " << dist << std::endl;

  // initialize the tree. Can choose from different algorithms, per knncolle library
  knncolle::VpTree<knncolle::distances::Euclidean, int, float> searcher(ndim, nobs, concatenated_data.data());
  //knncolle::AnnoyEuclidean<int, float> searcher(ndim, nobs, concatenated_data.data());
  //knncolle::Kmknn<knncolle::distances::Euclidean, int, float> searcher(ndim, nobs, concatenated_data.data());  

  if (m_verbose)
    std::cerr << " threads " << m_threads <<
      " OR flag " << orflag << " AND flag " << andflag << std::endl;
  
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
    
    // add the pheno flags
    std::vector<cy_uint> pflag_vec(neigh.size());
    for (size_t j = 0; j < neigh.size(); j++) {
      pflag_vec[j] = static_cast<IntCol*>(pflag_ptr.get())->GetNumericElem(neigh.at(j).first);
    }

    assert(pflag_vec.size() == neigh.size());

    // make the node
    node = CellNode(neigh, pflag_vec);

    // finally do tumor stuff
    node.sort_ascending_distance();
    float tumor_cell_count = 0;
    for (size_t j = 0; j < node.size(); j++) {
      CellFlag cellflag(node.m_flags.at(j));
      if (cellflag.testAndOr(orflag, andflag))
	tumor_cell_count++;
    }
    if (tumor_cell_count / static_cast<float>(node.size()) >= frac)
      static_cast<IntCol*>(cflag_ptr.get())->SetNumericElem(1, i);
    
  }// end for

#else
  std::cerr << "Warning: tumor call function requires including header library knncolle (https://github.com/LTLA/knncolle)" <<
    " and preprocessor directive -DHAVE_KNNCOLLE" << std::endl;
#endif
  
}

void CellTable::BuildKDTree() {

#ifdef HAVE_KDTREE
  
  pointVec points;
  const auto x_ptr = m_table.find("x");
  const auto y_ptr = m_table.find("y");
  assert(x_ptr != m_table.end());
  assert(y_ptr != m_table.end());

  // build the points vector
  for (size_t i = 0; i < CellCount(); i++) {
    points.push_back({ x_ptr->second->GetNumericElem(i), y_ptr->second->GetNumericElem(i)});
  }
  
  m_kdtree = new KDTree(points);

  /*
  std::cerr << "...querying tree" << std::endl;
  size_t i = 0;
  for (const auto& p : points) {
    i++;
    if (i % 10000 == 0)
      std::cerr << AddCommas(i) << std::endl;
    std::vector<size_t> ind = tree.neighborhood_indices(p, 10);
    if (i == 234) {
      std::cerr << "Query P: " << p[0] << "," << p[1] << std::endl;
      for (const auto& l : ind)
	std::cerr << "ind: " << l << " Point " << points[l][0] << "," << points[l][1] << std::endl;
    }
  }
  */  

#else
  std::cerr << "Warning: Unable to build KD-tree, need to include KDTree header library (https://github.com/crvs/KDTree) " <<
    " and add preprocessor directive -DHAVE_KDTREE" << std::endl;
#endif
  
}


int CellTable::RadialDensityKD(std::vector<cy_uint> inner, std::vector<cy_uint> outer,
			     std::vector<cy_uint> logor, std::vector<cy_uint> logand,
			     std::vector<std::string> label) {
#ifdef HAVE_KDTREE
  
  // check the radial geometry parameters
  assert(inner.size());
  assert(inner.size() == outer.size());
  assert(inner.size() == logor.size());
  assert(inner.size() == logand.size());
  assert(inner.size() == label.size());    

  if (m_verbose)  {
    for (size_t i = 0; i < inner.size(); i++) {
      std::cerr << "In " << inner.at(i) << " Out " << outer.at(i) <<
	" log " << logor.at(i) << " logand " << logand.at(i) << " label " << label.at(i) << std::endl;
    }
  }

  // get the flag column
  shared_ptr<IntCol> fc = std::dynamic_pointer_cast<IntCol>(m_table["pflag"]);
  assert(fc);
  assert(fc->size());

  // check sizes are good
  assert(fc->size() == CellCount());

  // check if already exists in the table
  for (const auto& l : label) {
    if (m_table.find(l) != m_table.end()) {
      throw std::runtime_error("Adding density meta column " + l + " already in table");
    }
  }
  
  // store the densities
  std::vector<std::shared_ptr<FloatCol>> dc(inner.size());
  for(auto& ptr : dc) {
    ptr = std::make_shared<FloatCol>();
    ptr->resize(fc->size());
  }

  // get the cell id colums
  auto id_ptr = m_table.at("id"); 
  assert(id_ptr->size() == fc->size());

  // flip it so that the id points to the zero-based index
  /*  std::unordered_map<uint32_t, size_t> inverse_lookup;
  if (m_verbose)
    std::cerr << "...done loading, doing the flip" << std::endl;
  uint32_t max_cell_id = 0;
  for (size_t i = 0; i < id_ptr->size(); ++i) {
    uint32_t value = id_ptr->GetNumericElem(i);
    inverse_lookup[value] = i;
    if (value > max_cell_id)
      max_cell_id = value;
      }*/

  if (m_verbose)
    std::cerr << "...radial density: starting loop" << std::endl;

  // build the tree
  if (m_verbose)
    std::cerr << "...building the KDTree" << std::endl;
  BuildKDTree();

  //
  const auto x_ptr = m_table.find("x");
  const auto y_ptr = m_table.find("y");
  assert(x_ptr != m_table.end());
  assert(y_ptr != m_table.end());

  // get max radius to compute on
  cy_uint max_radius = 0;
  for (const auto& r : outer)
    if (max_radius < r)
      max_radius = r;

  // pre-compute the bools
  std::vector<std::vector<int>> flag_result(inner.size(), std::vector<int>(fc->size(), 0));
  for (size_t i = 0; i < fc->size(); i++) {
    CellFlag mflag(fc->GetNumericElem(i));
    for (size_t j = 0; j < inner.size(); j++) {
      if ( (!logor[j] && !logand[j]) || mflag.testAndOr(logor[j], logand[j])) {
	flag_result[j][i] = 1; 
      }
    }
  }
  
    // loop the cells
#pragma omp parallel for num_threads(m_threads)
  for (size_t i = 0; i < fc->size(); i++) {
    
    // initialize the counts for each radial condition
    std::vector<float> cell_count(inner.size());

    float x1 = x_ptr->second->GetNumericElem(i);
    float y1 = y_ptr->second->GetNumericElem(i);
    point_t pt = { x1, y1 }; 

    // this will be inclusive of this point
    assert(m_kdtree);
    std::vector<size_t> inds = m_kdtree->neighborhood_indices(pt, max_radius);

    // loop the nodes connected to each cell
    for (const auto& n : inds) {

      float x2 = x_ptr->second->GetNumericElem(n);
      float y2 = y_ptr->second->GetNumericElem(n);
      float dx = x2 - x1;
      float dy = y2 - y1;
      float dist = std::sqrt(dx*dx + dy*dy);
      //float dist = euclidean_distance(x1, y1, x2, y2);
      
      // test if the connected cell meets the flag criteria
      // n.first is cell_id of connected cell to this cell
      for (size_t j = 0; j < inner.size(); j++) {
	//std::cerr << "inner " << inner[j] << " outer " << outer[j] << std::endl;
	//std::cerr << flag_result[j][n] << " logor " << logor[j] << " logand " << logand[j] << " flag " << fc->GetNumericElem(n) << std::endl;
	if (flag_result[j][n])
	  cell_count[j] += ((dist >= inner[j]) && (dist <= outer[j]));
      }
    }
    
    // calculate the density
    std::vector<float> area(inner.size());
    for (size_t j = 0; j < inner.size(); ++j) {
      float outerArea = static_cast<float>(outer[j]) * static_cast<float>(outer[j]) * 3.1415926535f;
      float innerArea = static_cast<float>(inner[j]) * static_cast<float>(inner[j]) * 3.1415926535f;
      area[j] = outerArea - innerArea;
    }

    // do the density calculation for each condition
    // remember, i is iterator over cells, j is over conditions
    for (size_t j = 0; j < area.size(); ++j) {
      if (!inds.empty()) {
	float value = cell_count[j] * 1000000 / area[j]; // density per 1000 square pixels
	dc[j]->SetNumericElem(value, i);
      } else {
	dc[j]->SetNumericElem(0, i); 
      }
    }

    /*
    // build the Cell object
    ////////////////////////
    Cell cell;
    cell.m_id   = id_ptr->GetNumericElem(i);
    cell.m_flag = fc->GetNumericElem(i);
    cell.m_x    = x_ptr->second->GetNumericElem(i);
    cell.m_y    = y_ptr->second->GetNumericElem(i);

    // fill the Cell data columns
    std::vector<ColPtr> col_ptr;
    for (const auto& t : m_header.GetDataTags()) {
      col_ptr.push_back(m_table.at(t.id));
    }

    for (const auto& c : col_ptr) {
      cell.m_cols.push_back(c->GetNumericElem(i));
    }

    // fill the Cell graph data
    CellNode node;    
    node.FillSparseFormat(cell.m_spatial_ids, cell.m_spatial_dist, cell.m_spatial_flags);

    // write
    (*m_archive)(cell);
    */
    if (m_verbose && i % 5000 == 0) {
      std::cerr << "...radial density: computing on cell " << AddCommas(i) << " and thread " << omp_get_thread_num() <<
	" and found densities ";
      for (size_t j = 0; j < dc.size(); j++) {
	std::cerr << dc.at(j)->GetNumericElem(i) << " ";
	if (j > 5)
	  break;
      }
      std::cerr << std::endl;
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
