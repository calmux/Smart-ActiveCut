double compute_beta(const vnl_matrix <double> & data, 
		    vnl_sparse_matrix <double> & nlinks_map);
int print_par(ParType par);
int print_gmm(const GMMType & gmm);
int print_vnl_matrix(const vnl_matrix<double> & mat, unsigned r, unsigned c);
int print_vnl_vec(const vnl_vector<double> & vec, unsigned N);

int save_gmm_labelmap(vnl_matrix<double> & gmm_labels, // initial gmm labels.
		      ParType & par,
		      const vnl_vector<unsigned> & alpha,
		      unsigned whatground,
		      ImageType3DU::Pointer lindexPtr,
		