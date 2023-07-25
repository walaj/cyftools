#include "cell_lda.h"

#include <iostream>

#include <ldaplusplus/LDA.hpp>
#include <ldaplusplus/LDABuilder.hpp>
#include <ldaplusplus/NumpyFormat.hpp>
#include <ldaplusplus/events/ProgressEvents.hpp>

CellLDA::CellLDA(size_t nt, size_t nc) {
  n_topics = nt;
  n_classes = nc;
}

void CellLDA::run(CellTable& tab) {

  n_topics = 10;
  n_classes = 10;
  n_docs = 1000;
  //size_t n_words = 100;



  std::cerr << "...created documents" << std::endl;
  
  // all the parameters below are the default and can be omitted
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

  std::cerr << "...lda plus plus ran" << std::endl;
  
  // add a listener to calculate and print the likelihood for every iteration
  // and a progress for every 128 documents (for every minibatch)
  double likelihood = 0;
  int count_likelihood = 0;
  int count = 0;
  lda.get_event_dispatcher()->add_listener(
					   [&likelihood, &count, &count_likelihood](std::shared_ptr<ldaplusplus::events::Event> ev) {
					     // an expectation has finished for a document
					     if (ev->id() == "ExpectationProgressEvent") {
					       count++; // seen another document
					       /*if (count % 128 == 0) {
						 std::cout << count << std::endl;
						 }*/
					       
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
					       std::cout << "Per document likelihood ~= "
							 << likelihood / count_likelihood << std::endl;
					       likelihood = 0;
					       count_likelihood = 0;
					       count = 0;
					     }
					   }
					   );


  // run the training for 15 iterations (we could also manually run each
  // iteration using partial_fit())
  lda.fit(X, y);
  
}
