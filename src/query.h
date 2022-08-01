
// the struct for attributes of connected components
typedef struct ConnectedCompInfor{
     int labelValue;
     int volumeSize;
     float avgPredictProb;
     float sumPredictProb;
} myConnectedCompInfor;

typedef struct CountCCP{
     int countNumVoxels;
     float curSumPredictProb;
} myCountCCP;

// the pair data structure for returning indices in sorting
typedef std::pair<int,int> mypair;
typedef std::pair<float,int> mypairFI;

// define unsigned char type for morphological processing
typedef unsigned char PixelType3DUC;
typedef itk::Image<PixelType3DUC , 3> ImageType3DUC;
typedef itk::ImageFileWriter< ImageType3DUC >  WriterType3DUC;


int logistic(vnl_vector<double> & score_map,
	     const vnl_matrix<double> & data,
	     vnl_sparse_matrix <double> & con_map,
	     const vnl_vector<unsigned> & alpha,
	     const ParType & par);


int logistic_init(vnl_vector<double> & score_map,
		  const vnl_matrix<double> & data,
		  vnl_sparse_matrix <double> & con_map, // can not set const due to a vnl bug.
		  const vnl_vector<unsigned> & alpha,
		  const ParType & par,
		  double smalleta,
		  unsigned n_iter);

void preprocessingBinaryPredicProb(ImageType3DI::Pointer binaryPredictProb);

void morphologicalProcessingTrueCandidatesVolume(ImageType3DUC::Pointer trueCandidatesVolume,
						 ImageType3DUC::Pointer trueCandidatesVolumeErosion, 