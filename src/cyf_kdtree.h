#pragma once
//
// cyf_kdtree.h — a header-only 2D KD-tree (nanoflann) that replaces mlpack's
// RangeSearch for cyftools' spatial queries. Removing mlpack also removes
// Armadillo, ensmallen, and BLAS/LAPACK from the build.
//
// KDTree2D::Search mirrors the shape of mlpack::RangeSearch::Search for a single
// 2D query point, so the call sites in cell_table_graph.cpp change only trivially
// (the query becomes a float[2], mlpack::Range becomes cyf::KRange). Distances are
// returned as ACTUAL Euclidean distances (nanoflann works in squared L2), because
// RadialDensityKD bins neighbors by true radius.
//
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "nanoflann.hpp"

namespace cyf {

// Mirrors mlpack::Range(lo, hi); only the upper bound is used as the radius.
struct KRange {
  double lo, hi;
  KRange(double l, double h) : lo(l), hi(h) {}
};

// nanoflann dataset adaptor over two parallel coordinate arrays.
struct PointCloud2D {
  const float* x = nullptr;
  const float* y = nullptr;
  std::size_t  n = 0;
  inline std::size_t kdtree_get_point_count() const { return n; }
  inline float kdtree_get_pt(std::size_t i, std::size_t dim) const { return dim == 0 ? x[i] : y[i]; }
  template <class BBOX> bool kdtree_get_bbox(BBOX&) const { return false; }
};

class KDTree2D {
public:
  // The coordinate arrays must outlive the tree (cyftools holds them in the
  // CellTable's x/y columns for the tree's lifetime).
  KDTree2D(const float* x, const float* y, std::size_t n)
      : m_cloud{x, y, n},
        m_index(2, m_cloud, nanoflann::KDTreeSingleIndexAdaptorParams(16)) {
    m_index.buildIndex();
  }

  // Drop-in for mlpack RangeSearch::Search(query, range, neighbors, distances)
  // for one query point: neighbors[0] = indices within r.hi (inclusive of the
  // query point itself if present), distances[0] = their Euclidean distances.
  void Search(const float query[2], const KRange& r,
              std::vector<std::vector<std::size_t>>& neighbors,
              std::vector<std::vector<double>>& distances) const {
    neighbors.assign(1, std::vector<std::size_t>{});
    distances.assign(1, std::vector<double>{});
    const float radius2 = static_cast<float>(r.hi * r.hi);  // nanoflann radius is squared L2
    std::vector<nanoflann::ResultItem<std::uint32_t, float>> matches;
    nanoflann::SearchParameters params;
    params.sorted = false;
    const std::size_t nfound = m_index.radiusSearch(query, radius2, matches, params);
    (void)nfound;
    neighbors[0].reserve(matches.size());
    distances[0].reserve(matches.size());
    for (const auto& m : matches) {
      neighbors[0].push_back(static_cast<std::size_t>(m.first));
      distances[0].push_back(std::sqrt(static_cast<double>(m.second)));
    }
  }

  std::size_t size() const { return m_cloud.n; }

private:
  PointCloud2D m_cloud;
  using Metric = nanoflann::L2_Simple_Adaptor<float, PointCloud2D>;
  nanoflann::KDTreeSingleIndexAdaptor<Metric, PointCloud2D, 2, std::uint32_t> m_index;
};

} // namespace cyf
