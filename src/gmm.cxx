
#include <common.h>
#include <gmm.h>
#include <utility.h>

// gmm parameter estimation. (M step).
int gmm_mstep(const vnl_matrix <double> & data, 
	      const vnl_matrix<double> & gmm_labels,
	      const vnl_vector<unsigned> & alpha,
	      unsigned whatground,
	      GMMType & gmm,
	      bool cov_shared,
	      unsigned cov_type) 

{
     unsigned label = 0;
     unsigned n_samples = data.rows();
     unsigned n_channels = data.cols();
     vnl_vector<double> sample_centered;

     // clear to zero before computation. 
     gmm.n_pts = 0;
     for (unsigned comp_id = 0; comp_id < gmm.n_comp; comp_id ++) {
	  gmm.comp[comp_id].label = 0;
	  gmm.comp[comp_id].numPts = 0;
	  gmm.comp[comp_id].mu.fill( 0 );
	  gmm.comp[comp_id].cov.fill( 0 );
     }

     // estimate mu. First compute sum. 
     for (unsigned sample_id = 0; sample_id < n_samples; sample_id ++) {
	  // assert (gmm_labels.get_row(sample_id).sum() == 1);
	  if (alpha[sample_id] == whatground) {
	       gmm.n_pts ++;
	       for (unsigned k = 0; k < gmm.n_comp; k ++) {
		    gmm.comp[k].numPts += gmm_labels[sample_id][k];
		    gmm.comp[k].mu +=  gmm_labels[sample_id][k] * data.get_row(sample_id);
	       }
	  }
     } 

     // comptue mean from the sum.
     for (unsigned k = 0; k < gmm.n_comp; k ++) {
     	  gmm.comp[k].mu /= gmm.comp[k].numPts;
     }

     // compute covariance matrix.
     if (cov_shared) {
	  vnl_matrix<double> big_cov(n_channels, n_channels, 0);
	  for (unsigned sample_id = 0; sample_id < n_samples; sample_id ++) {
	       if (alpha[sample_id] == whatground) {
		    for (unsigned k = 0; k < gmm.n_comp; k ++) {
			 sample_centered = data.get_row(sample_id) - gmm.comp[k].mu;
			 big_cov += gmm_labels[sample_id][k] * outer_product(sample_centered, sample_centered);
		    } // k
	       } // alpha
	  } // sample_id
	  big_cov /= gmm.n_pts;

	  for (unsigned k = 0; k < gmm.n_comp; k ++) {
	       gmm.comp[k].cov = big_cov;
	  }
     }

     else {	  // each comp has own cov mat.
	  for (unsigned sample_id = 0; sample_id < n_samples; sample_id ++) {
	       if (alpha[sample_id] == whatground) {
		    for (unsigned k = 0; k < gmm.n_comp; k ++) {
			 sample_centered = data.get_row(sample_id) - gmm.comp[k].mu;
			 gmm.comp[k].cov += gmm_labels[sample_id][k] * outer_product(sample_centered, sample_centered);
		    }
	       }
	  }

	  // normalize
	  for (unsigned k = 0; k < gmm.n_comp; k ++) {
	       gmm.comp[k].cov /= gmm.comp[k].numPts;
	  }
     }

     // regularize the covariance matrix. 
     if (cov_type == 0 ) {
	  // full covariance matrix. No regularization.
     }
     if (cov_type >= 1) {
	  // diagonal matrix. Set non-diagonal to zero.
	  for (unsigned k = 0; k < gmm.n_comp; k ++) {
	       for (unsigned i = 0; i < n_channels; i ++) {
		    for (unsigned j = 0; j < n_channels; j++) {
			 if (i != j) {