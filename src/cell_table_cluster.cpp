#include "cell_table.h"

#include <unordered_map>
#include <unordered_set>
#include <vector>

// DBSCAN clustering. Formerly used mlpack::DBSCAN (which pulled in Armadillo +
// BLAS/LAPACK); now implemented directly on the header-only nanoflann KD-tree
// (cyf::KDTree2D), so this file needs no heavy dependencies. (Filename kept for
// now to avoid build churn; it no longer uses mlpack.)

void CellTable::clusterDBSCAN(float epsilon,
			      size_t min_size,
			      size_t min_cluster_size
			      ) {

  validate();

  // collect the coordinates of marked cells (DBSCAN runs only on MARK_FLAG cells)
  std::vector<float> mx, my;
  std::vector<size_t> idx;            // map local index -> original cell index
  mx.reserve(CellCount());
  my.reserve(CellCount());
  idx.reserve(CellCount());
  for (size_t i = 0; i < m_cflag_ptr->size(); i++) {
    if (IS_FLAG_SET(m_cflag_ptr->at(i), MARK_FLAG)) {
      mx.push_back(m_x_ptr->at(i));
      my.push_back(m_y_ptr->at(i));
      idx.push_back(i);
    }
  }

  const size_t count = idx.size();
  if (count == 0) {
    std::cerr << "cyftools dbscan -- no cells to cluster. Did you run cyftools filter -M first?" << std::endl;
    return;
  }

  if (m_verbose)
    std::cerr << "...running dbscan (nanoflann) on " << AddCommas(count) <<
      " cells. Epsilon: " << epsilon << " Min points: " <<
      min_size << " min cluster size: " << min_cluster_size << std::endl;

  // KD-tree over the marked points (coordinate arrays must outlive the tree)
  cyf::KDTree2D tree(mx.data(), my.data(), count);

  const size_t NOISE = SIZE_MAX;
  std::vector<size_t> assignments(count, NOISE);
  std::vector<char>   visited(count, 0);

  auto region_query = [&](size_t p, std::vector<size_t>& out) {
    const float q[2] = { mx[p], my[p] };
    std::vector<std::vector<size_t>> nb;
    std::vector<std::vector<double>> ds;
    tree.Search(q, cyf::KRange(0.0, epsilon), nb, ds);   // radius search (incl. self)
    out = std::move(nb[0]);
  };

  size_t next_cluster = 0;
  std::vector<size_t> nbrs, qnbrs;
  for (size_t p = 0; p < count; p++) {
    if (visited[p]) continue;
    visited[p] = 1;

    region_query(p, nbrs);
    if (nbrs.size() < min_size) {
      assignments[p] = NOISE;          // may be reclaimed as a border point below
      continue;
    }

    // start a new cluster and expand it (classic DBSCAN seed-set growth)
    const size_t cid = next_cluster++;
    assignments[p] = cid;
    std::vector<size_t> seeds = nbrs;
    for (size_t k = 0; k < seeds.size(); k++) {
      const size_t q = seeds[k];
      if (!visited[q]) {
        visited[q] = 1;
        region_query(q, qnbrs);
        if (qnbrs.size() >= min_size)
          seeds.insert(seeds.end(), qnbrs.begin(), qnbrs.end());
      }
      if (assignments[q] == NOISE)
        assignments[q] = cid;
    }
  }

  // counts per cluster, to drop clusters smaller than min_cluster_size
  std::unordered_map<size_t, size_t> cluster_hist;
  for (size_t a : assignments)
    cluster_hist[a]++;

  // store cluster ids in a new column (1-based; 0 = unclustered / too small)
  FloatColPtr fc = std::make_shared<FloatCol>();
  fc->resize(CellCount());            // zero-filled

  std::unordered_set<size_t> assignment_set;
  for (size_t i = 0; i < count; i++) {
    const size_t a = assignments[i];
    if (a == NOISE || cluster_hist[a] < min_cluster_size) {
      fc->SetValueAt(idx[i], 0);
    } else {
      fc->SetValueAt(idx[i], static_cast<float>(a + 1));
      assignment_set.insert(a);
    }
  }

  if (m_verbose)
    std::cerr << "...produced " << AddCommas(assignment_set.size()) << " clusters" << std::endl;

  Tag ttag1(Tag::CA_TAG, "dbscan_cluster", "");
  assert(fc->size() == CellCount());
  AddColumn(ttag1, fc);
}
