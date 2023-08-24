#include "cell_synth.h"

#include <random>
#include "cell_table.h"

struct Point {
    double x;
    double y;
};

std::vector<Point> generateCluster(double centerX, double centerY, int numPoints, double sigma) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> d(0, sigma);

    std::vector<Point> cluster;
    for (int i = 0; i < numPoints; ++i) {
        cluster.push_back({centerX + d(gen), centerY + d(gen)});
    }
    return cluster;
}

CellSynth::CellSynth(size_t w, size_t h) : m_width(w), m_height(h) {

}

void CellSynth::WriteTable(const std::string& outfile) {

  assert(m_table);
  
  m_table->SetupOutputWriter(outfile);

  m_table->OutputTable();

}

void CellSynth::Clusters(size_t num_clusters,
			 size_t points_per_cluster,
			 double sigma,
			 size_t num_markers) {

  std::vector<Point> allPoints;

  size_t num_points = num_clusters * points_per_cluster;

  // fill up a fake table with length num_points
  m_table = new CellTable(num_points);
  
  std::random_device rd;
  std::mt19937 gen(m_seed);
  std::uniform_real_distribution<> disWidth(0, m_width);
  std::uniform_real_distribution<> disHeight(0, m_height);
  
  // make fake marker data
  for (size_t i = 0; i < num_markers; i++) {
    std::shared_ptr<FloatCol> mm = std::make_shared<FloatCol>();
    mm->resize(num_points);
    // add to the table
    Tag ttag(Tag::MA_TAG, "marker" + std::to_string(i+1), "");
    m_table->AddColumn(ttag, mm);
  }

  std::shared_ptr<FloatCol> x_ptr = std::make_shared<FloatCol>();
  std::shared_ptr<FloatCol> y_ptr = std::make_shared<FloatCol>();
  x_ptr->resize(m_table->CellCount());
  y_ptr->resize(m_table->CellCount());

  size_t j = 0;
  for (int i = 0; i < num_clusters; ++i) {
    double clusterCenterX = disWidth(gen);
    double clusterCenterY = disHeight(gen);
    std::cerr << " cluster at " << clusterCenterX << "," << clusterCenterY << std::endl;
    std::cerr << " width " << m_width << " hegiht " << m_height << std::endl;
    auto cluster = generateCluster(clusterCenterX, clusterCenterY, points_per_cluster, sigma);

    // add the points from this cluster
    for (const auto& c : cluster) {
      x_ptr->SetValueAt(j, c.x);
      y_ptr->SetValueAt(j, c.y);
      j++;
    }
  }

  // add the columns back
  m_table->AddICPXYColumn(x_ptr, "x");
  m_table->AddICPXYColumn(y_ptr, "y");  
  m_table->validate();
}
