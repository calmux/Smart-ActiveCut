int load_data(ImageType4DF::Pointer dataPtr,
	      ImageType3DI::Pointer maskPtr,
	      vnl_matrix<double> & D, 
	      ImageType3DU::Pointer lindexPtr);

int load_constraints(ImageType3DU::Pointer parcelPtr,
		     ImageType3DC::Pointer initPtr,
		     vnl_vector<unsigned> & hard_constraints,
		     const ParType & par);

int load_priors(ImageType3DU::Pointer parcelPtr,
		std::string prior