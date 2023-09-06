#include "cell_table.h"
#include "color_map.h"
#include "cell_selector.h"

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
    cairo_arc(cr, u1_ptr->getData().at(i)*cwidth/2+cwidth/2, u2_ptr->getData().at(i)*cwidth/2+cwidth/2, 0.5, 0.0, 2.0 * M_PI);
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
  validate();
  
  const std::vector<float>& x_data = m_x_ptr->getData(); // just a ptr, not a copy
  if (m_verbose) 
    std::cerr << "...adding " << AddCommas(x_data.size()) << " points on x" << std::endl;
  concatenated_data.insert(concatenated_data.end(), x_data.begin(), x_data.end());

  const std::vector<float>& y_data = m_y_ptr->getData(); // just a ptr, not a copy
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
    std::cerr << " OR flag " << orflag << " AND flag " << andflag << std::endl;
  
#pragma omp parallel for num_threads(m_threads)
  for (size_t i = 0; i < nobs; ++i) {
    
    // verbose printing
    if (i % 50000 == 0 && m_verbose)
      std::cerr << "...working on cell " << AddCommas(i) << 
	" K " << num_neighbors << " Dist: " << dist <<
	std::endl;
    
    JNeighbors neigh = searcher.find_nearest_neighbors(i, num_neighbors);
    
    CellNode node;
    
    // add the pheno flags
    std::vector<cy_uint> pflag_vec(neigh.size());
    for (size_t j = 0; j < neigh.size(); j++) {
      pflag_vec[j] = m_pflag_ptr->getData().at(neigh.at(j).first);
    }

    assert(pflag_vec.size() == neigh.size());

    // make the node
    node = CellNode(neigh, pflag_vec);

    // finally do tumor stuff
    node.sort_ascending_distance();
    float tumor_cell_count = 0;
    for (size_t j = 0; j < node.size(); j++) {
      CellFlag cellflag(node.m_flags.at(j));
      //      std::cerr << " cell flag " << node.m_flags.at(j) << " test OR " << orflag << " and " << andflag <<
      //	" test " << cellflag.testAndOr(orflag, andflag) << std::endl;
      if (cellflag.testAndOr(orflag, andflag))
	tumor_cell_count++;
    }
    if (tumor_cell_count / static_cast<float>(node.size()) >= frac)
      m_cflag_ptr->SetNumericElem(1, i);
    
  }// end for

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
      data(0, i) = m_x_ptr->getData().at(i);
      data(1, i) = m_y_ptr->getData().at(i);
    }
  
  //assert(ml_kdtree == nullptr);
  ml_kdtree = new mlpack::RangeSearch<mlpack::EuclideanDistance>(data);
  assert(ml_kdtree);
  
#else
    std::cerr << "Warning: Unable to build KD-tree, need to include MLPack library " << 
    " and add preprocessor directive -DHAVE_MLPACK" << std::endl;
#endif
  
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
  
  // store the densities
  std::vector<std::shared_ptr<FloatCol>> dc(inner.size());
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
      if ( (!logor[j] && !logand[j]) || mflag.testAndOr(logor[j], logand[j])) {
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
  
#pragma omp parallel for num_threads(m_threads) schedule(dynamic, 100)
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
      float value = cell_count[j] * 1000000 / area[j]; // density per 1000 square pixels
      if (normalize_local[j] != 0) // normalize to cell count in the radius
	value = total_cell_count[j] == 0 ? 0 : value / total_cell_count[j];
      if (normalize_global[j] != 0 && flag_count[j] != 0) // normalize to cell count in the image
	value = total_cell_count[j] == 0 ? 0 : value / flag_count[j] * total_image_cell_count;
      dc[j]->SetNumericElem(value, i);
    }
  
    if (m_verbose && i % 5000 == 0) {
      countr += 5000;
      std::cerr << std::fixed << "...cell " << std::setw(11) << AddCommas(i) 
		<< " thr " << std::setw(2) << omp_get_thread_num() 
		<< " %done: " << std::setw(3) << static_cast<int>(static_cast<float>(countr) / m_pflag_ptr->size() * 100) 
		<< " density: ";
      
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

#pragma omp parallel for num_threads(m_threads)
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

  return;
}
