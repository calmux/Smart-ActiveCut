
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
			      gmm.comp[k].cov(i,j) = 0;
			 }
		    }
	       }
	  } // k
     }
     if (cov_type >= 3) {
	  // identity matrix.
	  for (unsigned k = 0; k < gmm.n_comp; k ++) {
	       gmm.comp[k].cov.set_identity();
	  }
     }

     // Compute inverse of cov.
     for (unsigned k = 0; k < gmm.n_comp; k ++) {
	  // check singularity.
	  vnl_diag_matrix<double> diag_filler(n_channels, EPS);
	  // gmm.comp[k].cov.print(std::cout);
	  while (vnl_rank(gmm.comp[k].cov) < n_channels) {
	       gmm.comp[k].cov += diag_filler;
	       printf("gmm_mstep(): inv_cov[%i](FG): add on diag: %f\n", k, diag_filler[0]);
	       diag_filler *= 2;
	  }
	  gmm.comp[k].inv_cov = vnl_matrix_inverse<double>(gmm.comp[k].cov);	  
	  gmm.comp[k].det_cov = vnl_determinant(gmm.comp[k].cov);
     }

     // update pi. (do not need this code since we have a function for updating
     // the random variable pi given the hidden variables z.

     // for (unsigned k = 0; k < gmm.n_comp; k ++) {
     // 	  gmm.pi[k] = gmm.comp[k].numPts;
     // }
     // gmm.pi /= gmm.pi.sum();

     return 0;
}

// estimate the label posterior, given the gmm parameters. (E step).
double gmm_estep(const vnl_matrix <double> & data, 
		 vnl_sparse_matrix <double> & con_map,
		 vnl_matrix<double> & gmm_labels,
		 const vnl_vector<unsigned> & alpha,
		 unsigned whatground,
		 const GMMType & gmm)
{
     unsigned n_samples = gmm_labels.rows();
     unsigned n_channels = data.cols();
     unsigned n_comp = gmm.n_comp;
     double m = 0;
     // vnl_vector <double> zero_vec(gmm_labels.cols(), 0);

     vnl_vector<double> sample_c; // centered.
     vnl_vector<double> exp_term (gmm.n_comp, 0);
     vnl_matrix<double> gmm_oldlabels(gmm_labels);
     unsigned nbr_id = 0;
     vnl_sparse_matrix<double>::row con_map_row;

#pragma omp parallel for schedule(dynamic, 1000)  default(none) private(con_map_row, sample_c, exp_term, nbr_id, m) shared(alpha, whatground, gmm, con_map, n_channels, gmm_labels, data, n_samples, n_comp)
     for (unsigned sample_id = 0; sample_id < n_samples; sample_id ++) {     
	  exp_term.set_size(gmm.n_comp);
	  exp_term.fill(0);
     	  if (alpha[sample_id] == whatground) {
	       con_map_row = con_map.get_row(sample_id);		    
	       for (unsigned k = 0; k < gmm.n_comp; k ++) {
		    sample_c = data.get_row(sample_id) - gmm.comp[k].mu;
		    exp_term[k] = (- 0.5 * n_channels) * log(2 * PI)
			 - 0.5 * log(gmm.comp[k].det_cov) 
			 - 0.5 * dot_product(gmm.comp[k].inv_cov * sample_c, sample_c)
			 + log(gmm.pi(sample_id, k));

		    vnl_vector<double> z(gmm_labels.cols(), 0);
		    z[k] = 1;
		    for (vnl_sparse_matrix<double>::row::const_iterator col_iter = con_map_row.begin(); col_iter != con_map_row.end();  ++ col_iter) {
			 // go the neighbor sample_id
			 nbr_id = (*col_iter).first;
			 if (alpha[nbr_id] == whatground) { // in same ground.
			      exp_term[k] += gmm.beta * dot_product(z, gmm_labels.get_row(nbr_id));
			      // the above line does not take into account the
			      // patches. Int assume each patch is just a single voxel. To
			      // account for that, need to follow the build_nliks func.
			 }
		    } // con_iter
	       } // k
	       m = exp_term.max_value();
	       exp_term = exp_term - m;

	       for (unsigned k = 0; k < gmm.n_comp; k ++) {
		    gmm_labels(sample_id, k) = exp(exp_term[k]);
	       }	       
	       // normalize so it sum to 1. (unfilled values are just zero so it
	       // doesn't matter) Since each row may be for FG or BG. the col
	       // num may be larger than necessary. make sure unused elements
	       // does not contribute to the sum.
	       gmm_labels.scale_row(sample_id, 1/ gmm_labels.get_row(sample_id).extract(n_comp, 0).sum() );
	  } // alpha
     } // sample_id

     // sum of absolute values / 2 == the change of labels. 
     return (gmm_labels - gmm_oldlabels).array_one_norm() * 0.5;
}

int update_pi(const vnl_vector<unsigned> & alpha,
	      const vnl_matrix<double> & atlas,
	      GMMType & gmm,
	      double lambda)
{
     for (unsigned n = 0; n < atlas.rows(); n ++) {
	  // Important: For points not in this gmm, we still need a pi value
	  // (will be usedin the graphcuts and gmm_eval_ll func). So, no check
	  // for whatground.
	  for (unsigned k = 0; k < gmm.n_comp; k ++) {
	       gmm.pi(n, k) = pow(atlas(n, k) + EPS, lambda) * double(gmm.comp[k].numPts)/double(gmm.n_pts);
	  }

	  // normalize so pi sums to one. 
	  gmm.pi.scale_row(n, 1/gmm.pi.get_row(n).sum() );
     } // n
     return 0;
}

// query the log-likelihood of one data point. 
double gmm_eval_ll(const vnl_vector<double> & data,
		   const GMMType & gmm,
		   unsigned sample_id)
{
     vnl_vector<double> sample_c; // centered.
     vnl_vector<double> exp_term (gmm.n_comp, 0);
     unsigned n_channels = gmm.comp[0].mu.size();
     double sum_term = 0, m = 0;

     for (unsigned k = 0; k < gmm.n_comp; k ++) {
	  sample_c = data - gmm.comp[k].mu;
	  exp_term[k] = (- 0.5 * n_channels) * log(2 * PI)
	       - 0.5 * log(gmm.comp[k].det_cov) 
	       - 0.5 * dot_product(gmm.comp[k].inv_cov * sample_c, sample_c);
     }
     m = exp_term.max_value();
     exp_term = exp_term - m; // to avoid underflow.
     sum_term = 0;
     for (unsigned k = 0; k < gmm.n_comp; k ++) {
	  sum_term += gmm.pi(sample_id, k) * exp(exp_term[k]);
     }	       
     
     return log(sum_term) + m;
}

int kmeans(const vnl_matrix <double> & data,
	   const vnl_vector<unsigned> & alpha,
	   vnl_matrix<double> & gmm_labels, // initial gmm labels.
	   ParType & par,
	   unsigned whatground)
{
     unsigned n_comp = 0;
     if (whatground == ALPHA_FG) n_comp = par.gmm_fg.n_comp;
     else if (whatground == ALPHA_BG) n_comp = par.gmm_bg.n_comp;
     else {
	  printf("whatground must be either FG or BG!.\n");
	  exit (1);
     }

     // cluster centers. 
     std::vector<vnl_matrix<double> > cc(par.kmeansruns);
     vnl_vector<double> mean_ssr(par.kmeansruns, 0);
     double mean_ssr_old = 0;

     // init by kmeans++
#pragma omp parallel for default(none) private(mean_ssr_old) shared(mean_ssr, data, alpha, par, whatground, gmm_labels, n_comp, cc)

     for (unsigned r = 0; r < par.kmeansruns; r ++) {
	  printf("kmeans run: %i begin\n", r);
	  mean_ssr_old = 0;
	  printf("n_comp = %i, par.n_channels = %i\n", n_comp, par.n_channels);
	  (cc[r]).set_size(n_comp, par.n_channels);
	  (cc[r]).fill(0);
	  kmeans_init(data, alpha, cc[r], par, whatground, par.seed + r);

	  do {
	       // update labels. 
	       mean_ssr_old = mean_ssr[r];
	       mean_ssr[r] = kmeans_updatelabels(data, cc[r], alpha, gmm_labels, par, whatground);

	       // update cluster centers.
	       kmeans_updatecc(data, alpha, cc[r], gmm_labels, par, whatground);
	       // mean_ssr = compute_mean_ssr(data, alpha, cc, gmm_labels, par, whatground);

	       if (par.verbose >= 2) 