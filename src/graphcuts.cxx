#include <common.h>
#include <graph.h>
#include <gmm.h>
#include <utility.h>

int build_adjmat(ImageType3DU::Pointer lindexPtr,
		 const vnl_matrix <double> & data, 
		 vnl_sparse_matrix <double> & con_map,
		 ParType & par)
{
     ImageType3DU::IndexType lindexIdx, nbr_idx;

     // unsigned n_samples = data.rows();
     con_map.set_size(par.n_samples, par.n_samples);
     con_map.clear();

     // define neighborhood iterator
     typedef itk::ConstantBoundaryCondition<ImageType3DU>  BoundaryConditionType;
     typedef itk::NeighborhoodIterator< ImageType3DU, BoundaryConditionType > NeighborhoodIteratorType;
     NeighborhoodIteratorType::RadiusType radius;
     radius.Fill(1);
     NeighborhoodIteratorType lindexIt(radius, lindexPtr, lindexPtr->GetLargestPossibleRegion());
     unsigned int nei_set_array[] = {4, 10, 12, 14, 16, 22, // 6 neighborhood
				     1, 3, 5, 7, 9, 11, 15, 17, 19, 21, 23, 25, // 18 neighborhood
				     0, 2, 6, 8, 18, 20, 24, 26}; // 26 neighborhood

     if (par.n_nbrs != 6 && par.n_nbrs != 18 && par.n_nbrs != 26) {
	  printf("graphcuts(): number of neighbors must be 6, 18, or 26. Other values may give inacruate results!\n");
	  exit(1);
     }

     int cur_vox_id = 0, nbr_vox_id = 0;
     // if current voxel is boundary between patch i and j, add one on
     // con_map(i,j). 
     unsigned offset = 0;
     for (lindexIt.GoToBegin(); !lindexIt.IsAtEnd(); ++ lindexIt) {
	  if ( lindexIt.GetCenterPixel() > 0) {
	       cur_vox_id = lindexIt.GetCenterPixel() - 1;
	       
	       for (unsigned neiIdx = 0; neiIdx < par.n_nbrs; neiIdx ++) {
		    offset = nei_set_array[neiIdx];
		    nbr_vox_id =  lindexIt.GetPixel(offset) - 1;

		    // make sure pair (i,j) only count once, since they are
		    // ordered. Only upper triangular part of the matrix is
		    // filled.
		    if ( (nbr_vox_id >= 0) && (nbr_vox_id > cur_vox_id) ) {
			 con_map.put(cur_vox_id, nbr_vox_id, con_map(cur_vox_id, nbr_vox_id) + 1);
		    }
	       } // neiIdx
	  } // in mask
     }
     return 0;
}

int build_nlinks(const vnl_matrix <double> & data, 
		 const vnl_sparse_matrix <double> & con_map,
		 vnl_sparse_matrix <double> & nlinks_map,
		 ParTyp