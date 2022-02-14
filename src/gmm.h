
int gmm_mstep(const vnl_matrix <double> & data, 
	      const vnl_matrix<double> & gmm_labels, // initial gmm labels.
	      const vnl_vector<unsigned> & alpha,
	      unsigned whatground,
	      GMMType & gmm,
     	      bool cov_shared,
     	      unsigned cov_type);

// estimate the label posterior, given the gmm parameters. (E step).
double gmm_estep(const vnl_matrix <double> & data, 
		 vnl_sparse_matrix <double> & con_map,
		 vnl_matrix<double> & gmm_labels,
		 const vnl_vector<unsigned> & alpha,
		 unsigned whatground,
		 const GMMType & gmm);

// query the log-likelihood of one data point. 
double gmm_eval_ll(const vnl_vector<double> & data,
		   const GMMType & gmm,
		   unsigned sample_id);

int kmeans_init(const vnl_matrix <double> & data,
		const vnl_vector<unsigned> & alpha,
		vnl_matrix<double> & cc,
		ParType & par,
		unsigned whatground,
		unsigned seed);

int kmeans(const vnl_matrix <double> & data,
	   const vnl_vector<unsigned> & alpha,
	   vnl_matrix<double> & gmm_labels, // initial gmm labels.
	   ParType & par,
	   unsigned whatground);

int compute_dist(vnl_vector<double> & pdf,
		 const vnl_matrix<double> data, 
		 const vnl_vector<unsigned> & alpha,
		 const vnl_matrix<double> & cc,
		 unsigned cur_comp_id,
		 unsigned whatground);

double kmeans_updatelabels(const vnl_matrix <double> & data,
			   const vnl_matrix<double> & cc,
			   const vnl_vector<unsigned> & alpha,
			   vnl_matrix<double> & gmm_labels, // initial gmm labels.
			   ParType & par,
			   unsigned whatground);

int kmeans_updatecc(const vnl_matrix <double> & data,
		    const vnl_vector<unsigned> & alpha,
		    vnl_matrix<double> & cc,
		    const vnl_matrix<double> & gmm_labels, // initial gmm labels.
		    ParType & par,
		    unsigned whatground);

double compute_mean_ssr(const vnl_matrix <double> & data,