#include <common.h>
#include <gmm.h>
#include <graphcuts.h>
#include <utility.h>
#include <loadfiles.h>
#include <query.h>

namespace po = boost::program_options;

int main(int argc, char* argv[])
{
     std::string input_image_file, init_image_file, alpha_file, tlink_image_file, nlink_image_file, mask_file, gmm_fg_file, gmm_bg_file, priorfg_file, priorbg_file, fg_label_file, bg_label_file, score_file, cand_file;
     
     // all parameters including GMM.
     ParType par;
     unsigned maxem = 0;;
     unsigned cov_type = 0;
     bool cov_shared = false;
     bool noactlearn = false;
     par.baseline = false;
     

     // program options.
     po::options_description mydesc("Options can only used at commandline");
     mydesc.add_options()
	  ("help,h", "Grab cut segmentation of 3D multi-channel traumatic brain injury (TBI) images.")

	  ("ncompbg,b", po::value<unsigned>(&par.gmm_bg.n_comp)->default_value(4),
	   "number of Gaussian components in background.")

	  ("ncompfg,f", po::value<unsigned>(&par.gmm_fg.n_comp)->default_value(2),
	   "number of Gaussian components in foreground.")

	  ("maxem,m", po::value<unsigned>(&maxem)->default_value(30),
	   "Max number of EM iterations.")

	  ("kmeansruns,r", po::value<unsigned>(&par.kmeansruns)->default_value(5),
	   "number of kmeans runs.")

	  ("betaf", po::value<double>(&par.gmm_fg.beta)->default_value(1),
	   "smoothness (MRF) constraint for foreground gmm segmentation.")

	  ("betab", po::value<double>(&par.gmm_bg.beta)->default_value(1),
	   "Smoothness (MRF) constraint for background gmm segmentation.")

	  ("gamma,g", po::value<double>(&par.gamma)->default_value(1),
	   "alpha smoothness constraint.")

	  ("beta0", po::value<double>(&par.beta0)->default_value(1),
	   "An additional parameter in front of beta. ")

	  ("eta", po::value<double>(&par.eta)->default_value(1),
	   "Smoothness for the query score")

	  ("neighbors,n", po::value<unsigned>(&par.n_nbrs)->default_value(6),
	   "number of neighbors of a voxel for graphcuts N-Links. Legal values are 6, 18 and 26. ")

	  ("covtype", po::value<unsigned>(&cov_type)->default_value(1),
	   "Type of covariance matrix. 0 for full matrix, 1 for diagonal, 3 for identity.")

	  ("covshared", po::bool_switch(&cov_shared), 
	   "Whether all components of GMM has same covariance matrix. Default is no.")

	  ("data,d", po::value<std::string>(&input_image_file),
	   "Input all-channel image file. A 4D gipl or nii or nii.gz file.")

	  ("init,i", po::value<std::string>(&init_image_file),
	   "A 3D volume image giving the user's initial input. Inside the box will be init'd unkown, and outside will be init'd as background.")

	  ("priorfg", po::value<std::string>(&priorfg_file)->default_value("priorfg.nii.gz"),
	   "4D file for the foreground (bleeding/edema) prior probability. Give a zero image if no prior available.")

	  ("priorbg", po::value<std::string>(&priorbg_file)->default_value("priorbg.nii.gz"),
	   "4D file for the background (GM/WM/CSF) prior probability. Give a zero image if no prior available.")

	  ("lambdafg", po::value<double>(&par.lambda_fg)->default_value(1),
	   "the confidence on the foreground atlas. choose a lambda > 1 means more confidence on the atlas. ")

	  ("lambdabg", po::value<double>(&par.lambda_bg)->default_value(1),
	   "the confidence on the background atlas. choose a lambda > 1 means more confidence on the atlas. ")

	  ("mask,p", po::value<std::string>(&mask_file)->default_value("mask.nii.gz"),
	   "Mask file. Outside the brain is zero, inside is one.")

	  ("alphafile,a", po::value<std::string>(&alpha_file)->default_value("alpha.nii.gz"),
	   "The output segmentation file. Will be binary volume with same size as input data.")

	  ("gmmfg", po::value<std::string>(&gmm_fg_file)->default_value("gmmfg.nii.gz"),
	   "Foreground GMM posterior. ")

	  ("fglabel", po::value<std::string>(&fg_label_file)->default_value("gmmfg_label.nii.gz"),
	   "Foreground GMM segmentation label image.")

	  ("gmmbg", po::value<std::string>(&gmm_bg_file)->default_value("gmmbg.nii.gz"),
	   "Background GMM posterior.")

	  ("bglabel", po::value<std::string>(&bg_label_file)->default_value("gmmbg_label.nii.gz"),
	   "Background GMM segmentation label image.")

	  ("queryscores,q", po::value<std::string>(&score_file)->default_value("queryscores.nii.gz"),
	   "query scores file.")

	  ("cand,c", po::value<std::string>(&cand_file)->default_value("cand.nii.gz"),
	   "output candidate component for user to check.")

	  ("predth", po::value<float>(&par.pred_th)->default_value(0.75),
	   "Threshold used to get connected components from predictive probability.")

	  ("qscoreth", po::value<float>(&par.qscore_th)->default_value(1.3),
	   "Threshold for query score. Below that score, the active learning will stop.")

	  ("baseline", po::bool_switch(&par.baseline), 
	   "set baseline for baseline testing that allow background --> foreground and no active learning.")

	  ("noactlearn", po::bool_switch(&noactlearn), 
	   "If doing self-training and active learning after EM and graphcuts converges. default is yes.")

	  ("logmax", po::value<unsigned>(&par.logimax)->default_value(5),
	   "Max number of iterations in logistic() func.")

	  ("seed,s", po::value<unsigned>(&par.seed)->default_value(0),
	   "Seed for random number generator.")

	  ("verbose,v", po::value<unsigned>(&par.verbose)->default_value(0),
	   "verbose level in [0, 3]. ")
	  ;

     po::variables_map vm;        
     po::store(po::parse_command_line(argc, argv, mydesc), vm);
     po::notify(vm);    

     try {
	  if ( (vm.count("help")) | (argc == 1) ) {
	       std::cout << "Usage: grabcut [options]\n";
	       std::cout << mydesc << "\n";
	       return 0;
	  }
     }
     catch(std::exception& e) {
	  std::cou