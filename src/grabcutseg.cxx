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

	  ("eta", po::value<double>(&par.et