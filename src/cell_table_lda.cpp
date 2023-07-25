#include "cell_table.h"
#include "lda_model.h"
#include "color_map.h"

#ifdef HAVE_CAIRO
#include "cairo/cairo.h"
#include "cairo/cairo-pdf.h"
#endif


#ifdef HAVE_LDAPLUSPLUS

std::ostream& operator<<(std::ostream& os, const LDAModel& result) {

    os << "alpha size: " << result.alpha.size() << "\n";
    os << "beta size: " << result.beta.size() << "\n";
    if(!result.beta.empty()) {
        os << "beta[0] size: " << result.beta[0].size() << "\n";
    }
    return os;
}

// Constructor that takes Eigen objects and converts them to std::vector
LDAModel::LDAModel(const Eigen::VectorXd& alpha_eigen,
		   const Eigen::MatrixXd& beta_eigen,
		   const std::vector<std::string>& markers) {

    markers_to_run = markers;
    
    // Convert Eigen::VectorXd to std::vector<double>
    alpha = std::vector<double>(alpha_eigen.data(), alpha_eigen.data() + alpha_eigen.size());

    // Convert Eigen::MatrixXd to std::vector<std::vector<double>>
    beta.resize(beta_eigen.rows(), std::vector<double>(beta_eigen.cols()));
    for(int i = 0; i < beta_eigen.rows(); ++i) 
      for(int j = 0; j < beta_eigen.cols(); ++j)
        beta[i][j] = beta_eigen(i, j);
  }

void CellTable::LDA_write_model(const std::string& model_out) const {

  // store the model
  if (!model_out.empty()) {
    
    if (m_verbose) 
      std::cerr << "...storing the model to " << model_out << std::endl;

    assert(m_ldamodel);
    m_ldamodel->JSONArchive(model_out);
    
    //std::ofstream ofs(model_out);
    //cereal::PortableBinaryOutputArchive oarchive(ofs);
    //oarchive(m_ldamodel);
  }

  return;
}

void CellTable::LDA_load_model(const std::string& model_file) {

  if (!check_readable(model_file))  {
    std::cerr << "Error: Model file " << model_file << " not readable/exist" << std::endl;
    return;
  }

  m_ldamodel = new LDAModel();
  
  std::ifstream is(model_file);
  cereal::JSONInputArchive iarchive(is);
  iarchive(*m_ldamodel);
  assert(m_ldamodel);

  if (m_verbose)
    std::cerr << *m_ldamodel << std::endl;
  
  //std::ifstream is(model_file, std::ios::binary);
  //cereal::PortableBinaryInputArchive iarchive(is);
  //iarchive(m_ldamodel);
  
  return;
}

void CellTable::LDA_score_cells(const std::string& pdffile,
				int topic_highlight,
				float cont_cutoff) { 

  if (m_verbose)
    std::cerr << "...scoring cells with LDA model" << std::endl;

  if (!m_ldamodel || !m_ldamodel->initialized()) {
    std::cerr << "Error: Need to create or load LDA model first" << std::endl;
    return;
  }
  
  // Build the model
  std::shared_ptr<ldaplusplus::parameters::ModelParameters<double>> model = std::make_shared<ldaplusplus::parameters::ModelParameters<double>>(m_ldamodel->toModelParameters());
  ldaplusplus::LDA<double> lda_r = ldaplusplus::LDABuilder<double>()
    .initialize_topics_from_model(model);

  // Create an Eigen matrix to hold the document data
  size_t n_words = m_ldamodel->markers_to_run.size();
  size_t n_docs = CellCount();
  assert(n_words > 0);
  assert(n_docs > 0);
  Eigen::MatrixXi X(n_words, n_docs);
  
  // Fill the matrix with the marker data
  int i = 0;
  assert(m_ldamodel);
  for (const auto& s : m_ldamodel->markers_to_run) {
    assert(m_table.find(s) != m_table.end());
    shared_ptr<FloatCol> fc = std::dynamic_pointer_cast<FloatCol>(m_table.at(s));
    assert(fc->size());
    
    for (int j = 0; j < CellCount(); j++) {
      X(i, j) = static_cast<int>(fc->GetNumericElem(j));
    }
    i++;
  }
  
  // predict
  if (m_verbose)
    std::cerr << "...getting topic assignment using existing LDA model" << std::endl;
  
  Eigen::MatrixXd Zr = lda_r.transform(X); // cols is docs, rows is topics
  for (int i = 0; i < Zr.cols(); ++i) {
    double sum = Zr.col(i).sum();
    Zr.col(i) = Zr.col(i) /= sum;
  }

  // add to the output to the data
  for (size_t i = 0; i < Zr.rows(); ++i) {
    std::shared_ptr<FloatCol> fc = make_shared<FloatCol>();
    fc->resize(Zr.cols());
    for (size_t j = 0; j < Zr.cols(); j++) {
      fc->SetNumericElem(Zr(i, j), j);
    }
    //
    Tag ttag1(Tag::CA_TAG, "topic" + std::to_string(i+1),"");
    AddColumn(ttag1, fc);
  }
  
  //////
  // PLOT
  //////
  if (!pdffile.empty()) {
#ifdef HAVE_CAIRO

    const float scale_factor = 1.0f;
    const float radius_size = 6.0f;
    const float alpha_val = 0.5f;
    
    // get the x y coordinates of the cells
    const auto x_ptr = m_table.find("x");
    const auto y_ptr = m_table.find("y");
    assert(x_ptr != m_table.end());
    assert(y_ptr != m_table.end());
    
    // open the PDF for drawing
    const int width  = x_ptr->second->Max();
    const int height = y_ptr->second->Max();    
    //cairo_surface_t *surface = cairo_pdf_surface_create(pdffile.c_str(), width, height);
    //cairo_t *cr = cairo_create(surface);

    // open PNG for drawing
    cairo_surface_t *surfacep = cairo_image_surface_create (CAIRO_FORMAT_ARGB32, width*scale_factor, height*scale_factor);
    cairo_t *crp = cairo_create (surfacep);
    cairo_set_source_rgb(crp, 0, 0, 0); // background color

    ColorMap cm = getColorMap(m_ldamodel->getNumTopics());

    topic_highlight -= 1;
    if (topic_highlight >= 0) {
      if (topic_highlight > m_ldamodel->getNumTopics())
	throw std::runtime_error("Error: selected topic greater than number of topics");
      for (size_t i = 0; i < cm.size(); i++) {
	if (i != topic_highlight)
	  cm[i] = {255,255,255};
	else
	  cm[i] = {255,0,0};
      }
    }
    
    // setup for drawing loop
    assert(Zr.cols() == CellCount());
    assert(Zr.rows() == m_ldamodel->getNumTopics());
    constexpr float TWO_PI = 2.0 * M_PI;
    
    // loop and draw
    size_t count = 0;
    for (size_t j = 0; j < Zr.cols(); j++) { // loop the cells

      // Draw the arc segment
      const float x = x_ptr->second->GetNumericElem(j);
      const float y = y_ptr->second->GetNumericElem(j);
      
      float start_angle = 0.0;
      
      // Loop over rows (topics)
      for(int i=0; i < Zr.rows(); i++) {
	
        // Calculate how much this topic contributes to the document
        float contribution = Zr(i, j);
	
	if (contribution < cont_cutoff)
	  continue;
	
	count++;
	
        // Calculate the arc length corresponding to this contribution
        float arc_length = contribution * TWO_PI;

        // Set color for this topic
        Color c = cm[i];
        //cairo_set_source_rgba(cr, c.redf(), c.greenf(), c.bluef(), 0.3);
        //cairo_arc(cr, x, y, 5, start_angle, start_angle + arc_length);

	//png
        cairo_set_source_rgba(crp, c.redf(), c.greenf(), c.bluef(), alpha_val);
        cairo_arc(crp, x*scale_factor, y*scale_factor, radius_size, start_angle, start_angle + arc_length);
	
	//if (j % 10000 == 0) {
	//  std::cerr << "X " << x << " Y " << y << " contribtion " << contribution <<
	//    " arc length " << arc_length << std::endl;
	//}
	
        // Close the path
        //cairo_line_to(cr, x, y);
        //cairo_fill(cr);

        cairo_line_to(crp, x*scale_factor, y*scale_factor);
        cairo_fill(crp);

	
        // Update the start angle for the next iteration
        start_angle += arc_length;
      }
    }
    
    //std::cerr << " drawing " << AddCommas(count) << " fragments " << std::endl;
    
    // Clean up and close
    //cairo_destroy(cr);
  //cairo_surface_destroy(surface);

    std::vector<std::string> labels;
    for (size_t i = 0; i < Zr.rows(); i++) {
      labels.push_back("Topic " + std::to_string(i+1));
    }

    /// LEGEND
    int legend_width = 1000*scale_factor; // Width of each color box
    int legend_height = 300*scale_factor; // Height of each color box
    int font = 216*scale_factor;
    // Starting position for the legend in the upper right corner
    int legend_padding = 20;
    int legend_x = width*scale_factor - legend_width - legend_padding;
    int legend_y = legend_padding;
    add_legend_cairo(crp, font, legend_width, legend_height,
		     legend_x, legend_y, cm, labels);
    /*    
    for (int i = 0; i < 10; ++i) {
      
      // Set color for this entry
      Color c = cm[i];
      cairo_set_source_rgb(crp, c.redf(), c.greenf(), c.bluef());
      
      // Draw the color box
      cairo_rectangle(crp, legend_x, legend_y + i*(legend_height+legend_padding), legend_width, legend_height);
      cairo_fill(crp);
      
      // Draw the label
      cairo_set_source_rgb(crp, 0, 0, 0); // Set color to black for the text
      cairo_move_to(crp, legend_x, legend_y + i*(legend_height+legend_padding) + legend_height);
      cairo_show_text(crp, labels[i].c_str());
    }
    */
    cairo_destroy (crp);
    cairo_surface_write_to_png (surfacep, pdffile.c_str());
    cairo_surface_destroy (surfacep);
    
#else
    std::cerr << "Unable to make PDF -- need to include / link Cairo library" << std::endl;
#endif
  }
  
  /*  int j = 0; //document to print
  std::cerr << "cols " << Zr.cols() << " rows " << Zr.rows() << std::endl;
  for (int i = 0; i < Zr.rows(); ++i) {
    std::cerr << "Score of topic (re)" << i
	      << " for document " << j 
	      << ": " << Zr(i, j) << std::endl;
  }
  std::cerr << "..end" << std::endl;
  */
  
  // // Assuming that each row of beta represents a topic
  // for(int topicIdx = 0; topicIdx < model->beta.rows(); topicIdx++) {
    
  //   // Get the word distribution for the current topic
  //   Eigen::Matrix<double, 1, Eigen::Dynamic> wordDistribution = model->beta.row(topicIdx);

  //   std::cerr << "...Topic " << (topicIdx) << " rows " << wordDistribution.rows() << std::endl;
  //   // Loop over each word's probability in the topic
  //   for(int wordIdx = 0; wordIdx < wordDistribution.cols(); wordIdx++) {
  //     double probability = wordDistribution[wordIdx];
  //     std::cerr << "Word " << wordIdx << " probability: " << probability << "\n";
  //   }
  // }
  
  
}

void CellTable::LDA_create_model(const std::vector<std::string>& marker_cols,
				 size_t n_topics,
				 size_t n_iterations) {

  size_t n_words = marker_cols.size();
  size_t n_docs = CellCount();

  if (m_verbose)
    std::cerr << "...setting up LDA with " << AddCommas(n_docs) <<
      " cells; " <<
      n_words << " columns; " <<
      n_topics << " topics" << std::endl;
  
  // Create an Eigen matrix to hold the document data
  Eigen::MatrixXi X(n_words, n_docs);
  
  // Fill the matrix with the marker data
  int i = 0;
  for (const auto& s : marker_cols) {
    assert(m_table.find(s) != m_table.end());
    shared_ptr<FloatCol> fc = std::dynamic_pointer_cast<FloatCol>(m_table.at(s));
    assert(fc->size());
    
    for (int j = 0; j < CellCount(); j++) {
      X(i, j) = static_cast<int>(fc->GetNumericElem(j));
    }
    i++;
  }

  // all the parameters below are the default and can be omitted
  ldaplusplus::LDA<double> lda = ldaplusplus::LDABuilder<double>()
    .initialize_topics_random(
        X.rows(),   // X.rows() is the number of words in the vocab
        n_topics  // how many topics we want to infer
    )
    .set_iterations(n_iterations)
    .set_workers(m_threads);
  
  /*
  ldaplusplus::LDA<double> lda = ldaplusplus::LDABuilder<double>()
    .set_fast_supervised_e_step(
				10,   // expectation iterations
				1e-2, // expectation tolerance
				1,    // C parameter of fsLDA (see the paper)
				0.01, // percentage of documents to compute likelihood for
				42    // the randomness seed
				)
    .set_fast_supervised_online_m_step(
				       n_classes,   // number of classes in the dataset
				       0.01, // the regularization penalty
				       128,  // the minibatch size
				       0.9,  // momentum for SGD training
				       0.01, // learning rate for SGD
				       0.9   // weight for the LDA natural gradient
				       )
    .initialize_topics_seeded(
			      X,  // the documents to seed from
			      n_topics, // the number of topics
			      42  // the randomness seed
			      )
    .initialize_eta_zeros(n_topics) // initialize the supervised parameters
    .set_iterations(15);
  */

  
  // add a listener to calculate and print the likelihood for every iteration
  // and a progress for every 128 documents (for every minibatch)
  double likelihood = 0;
  int count_likelihood = 0;
  int count = 0;
  lda.get_event_dispatcher()->add_listener(
					   [this, &likelihood, &count, &count_likelihood](std::shared_ptr<ldaplusplus::events::Event> ev) {
					     // an expectation has finished for a document
					     if (ev->id() == "ExpectationProgressEvent") {
					       count++; // seen another document
					       
					       // aggregate the likelihood if computed for this document
					       auto expev =
						 std::static_pointer_cast<ldaplusplus::events::ExpectationProgressEvent<double> >(ev);
					       if (expev->likelihood() < 0) {
						 likelihood += expev->likelihood();
						 count_likelihood ++;
					       }
					     }
					     
					     // A whole pass from the corpus has finished print the approximate per
					     // document likelihood and reset the counters
					     else if (ev->id() == "EpochProgressEvent") {
					       if (this->m_verbose)
						 std::cerr << "Per document likelihood ~= "
							   << likelihood / count_likelihood << std::endl;
					       likelihood = 0;
					       count_likelihood = 0;
					       count = 0;
					     }
					   }
					   );

  // run the training for 15 iterations (we could also manually run each
  // iteration using partial_fit())
  if (m_verbose)
    std::cerr << "...running LDA" << std::endl;
  lda.fit(X);

  // Extract the top words of the unsupervised model
  const std::shared_ptr<ldaplusplus::parameters::ModelParameters<> > model =
    lda.model_parameters<ldaplusplus::parameters::ModelParameters<> >();

  // add the model to this object
  m_ldamodel = new LDAModel(model->alpha, model->beta, marker_cols);
  m_ldamodel->cmd = m_cmd;

}

//end HAVE_LDAPLUSPLU
#endif 

