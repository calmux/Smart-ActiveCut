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
     NeighborhoodIteratorT