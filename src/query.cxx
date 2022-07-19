
#include <common.h>
#include <gmm.h>
#include <query.h>
#include <loadfiles.h>
// include libraries and define some struts for CCP
#include <vector>       // std::vector

// the comparator for sorting
bool myCompFunction (int i,int j) { return (i>j); }

// the pair data structure for returning indices in sorting
bool myComparatorforCcp ( const mypair& l, const mypair& r)
   { return l.first > r.first; }
bool myComparatorforCcpFI ( const mypairFI& l, const mypairFI& r)
   { return l.first > r.first; }
 

// for mathematical morphology operations - erosion
#include "itkBinaryErodeImageFilter.h"
#include "itkBinaryDilateImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkConnectedComponentImageFilter.h"
// for duplicate image
#include "itkImageDuplicator.h"

using namespace std;

int logistic(vnl_vector<double> & score_map,
	     const vnl_matrix<double> & data,
	     vnl_sparse_matrix <double> & con_map, // can not set const due to a vnl bug.
	     const vnl_vector<unsigned> & alpha,
	     const ParType & par)
{
     unsigned nbr_id = 0;
     vnl_sparse_matrix<double>::row con_map_row;
     vnl_vector<double> exp_term (2, 0);
     vnl_vector<double> score_map_old(score_map);
     double changed_score = 1e7;
     unsigned iter = 0;
     while(changed_score > 0.001 && iter < par.logimax) {
	  iter ++;
	  score_map_old = score_map; // save previous scores.
	  for (unsigned n = 0; n < par.n_samples; n ++) {
	       // compute difference of log-likelihood.
	       exp_term[0] = gmm_eval_ll(data.get_row(n), par.gmm_fg, n);
	       exp_term[1] = gmm_eval_ll(data.get_row(n), par.gmm_bg, n);

	       // compute diff of prior.  we don't bother looping the 2
	       // iterations FG and BG. just repleat it twice.
	       con_map_row = con_map.get_row(n);		    
	       for (vnl_sparse_matrix<double>::row::const_iterator col_iter = con_map_row.begin(); col_iter != con_map_row.end();  ++ col_iter) {
		    nbr_id = (*col_iter).first;
		    // if score = 1 (FG)
		    exp_term[0] += par.eta * (1 * score_map[nbr_id] + 0 *(1-score_map[nbr_id]));
		    // the above line does not take into account the
		    // patches. Int assume each patch is just a single voxel. To
		    // account for that, need to follow the build_nliks func.
	       }

	       for (vnl_sparse_matrix<double>::row::const_iterator col_iter = con_map_row.begin(); col_iter != con_map_row.end();  ++ col_iter) {
		    nbr_id = (*col_iter).first;
		    // if score == 0 (BG)
		    exp_term[1] += par.eta * (0 * score_map[nbr_id] + 1 * (1-score_map[nbr_id]));
	       }

	       score_map[n] = exp_term[0] - exp_term[1];
	       score_map[n] = 1 / (1 + exp(- score_map[n]));
	  } // n i.e. sample id.
	  changed_score = (score_map_old - score_map).two_norm() / score_map.two_norm();
	  if (par.verbose >= 0) {
	       printf("logistic(): iteration %i, changed_score: %f\n", iter, changed_score);
	  }
     }
     return 0;
}


int logistic_init(vnl_vector<double> & score_map,
		  const vnl_matrix<double> & data,