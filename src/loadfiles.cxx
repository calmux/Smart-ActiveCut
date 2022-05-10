
#include <common.h>
#include <gmm.h>
#include <graphcuts.h>
#include <utility.h>

int load_data(ImageType4DF::Pointer dataPtr,
	      ImageType3DI::Pointer maskPtr,
	      vnl_matrix<double> & D,
	      ImageType3DU::Pointer lindexPtr)
{
     unsigned spixel_id = 0;
     ImageType4DF::SizeType dataSize = dataPtr->GetLargestPossibleRegion().GetSize();
     IteratorType4DF dataIt(dataPtr, dataPtr->GetLargestPossibleRegion());
     ImageType4DF::IndexType dataIdx;

     ImageType3DI::SizeType maskSize = maskPtr->GetLargestPossibleRegion().GetSize();
     IteratorType3DI maskIt(maskPtr, maskPtr->GetLargestPossibleRegion());
     ImageType3DI::IndexType maskIdx;

     IteratorType3DU lindexIt(lindexPtr, lindexPtr->GetLargestPossibleRegion());

     unsigned n_channels = dataSize[3];
     unsigned n_samples = 0;

     // calculate number of voxels in brain mask. this code is stupid but it works. 
     for (maskIt.GoToBegin(); !maskIt.IsAtEnd(); ++ maskIt) {
	  if (maskIt.Get() > 0) {
	       n_samples ++;
	  }
     }

     printf("load_data(): Total number of voxels in mask: %i.\n", n_samples);
     
     D.set_size(n_samples, n_channels);
     D.fill( 0 );
     
     unsigned n = 0; // voxel linear index. 
     for (maskIt.GoToBegin(), lindexIt.GoToBegin(); !maskIt.IsAtEnd(); ++ maskIt, ++ lindexIt) {
	  if (maskIt.Get() > 0) {
	       maskIdx = maskIt.GetIndex();
	       dataIdx[0] = maskIdx[0];
	       dataIdx[1] = maskIdx[1];
	       dataIdx[2] = maskIdx[2];
	       for (unsigned channel_id = 0; channel_id < n_channels; channel_id ++) {
		    dataIdx[3] = channel_id;
		    D[n][channel_id] += dataPtr->GetPixel(dataIdx);
	       }
	       n++;
	       lindexIt.Set(n);
	  } // in mask
     }
     return 0;
}

int load_constraints(ImageType3DU::Pointer lindexPtr,
		     ImageType3DC::Pointer initPtr,
		     vnl_vector<unsigned> & hard_constraints,
		     const ParType & par)
{
     unsigned vox_id = 0;
     ImageType3DU::SizeType lindexSize = lindexPtr->GetLargestPossibleRegion().GetSize();
     IteratorType3DU lindexIt(lindexPtr, lindexPtr->GetLargestPossibleRegion());
     ImageType3DU::IndexType lindexIdx;

     IteratorType3DC initIt(initPtr, initPtr->GetLargestPossibleRegion());
     ImageType3DC::IndexType initIdx;

     // the code is simplified by removing supervoxel support.
     for (lindexIt.GoToBegin(), initIt.GoToBegin(); !initIt.IsAtEnd(); ++ lindexIt, ++ initIt) {
	  if (lindexIt.Get() > 0) { // in mask
	       vox_id = lindexIt.Get() - 1;
	       if (initIt.Get() == HC_FG) { // fg
		    hard_constraints[vox_id] = HC_FG;
	       }
	       // either init'd background, or user given background.
	       else if (initIt.Get() == HC_BG) { 
		    hard_constraints[vox_id] = HC_BG;
	       }
	       else if  (initIt.Get() == HC_BG_NEW) {
		    hard_constraints[vox_id] = HC_BG_NEW;
	       }
	       else { // unknown region will be init'd as foreground!
		    hard_constraints[vox_id] = HC_UNKNOWN;
	       }
	  } // in mask
	  
     }
     return 0;
}

int load_priors(ImageType3DU::Pointer lindexPtr,
		std::string prior_file,
		vnl_matrix<double> & priors,
		const ParType & par)
{
     unsigned vox_id = 0;
     ImageType3DU::SizeType lindexSize = lindexPtr->GetLargestPossibleRegion().GetSize();
     IteratorType3DU lindexIt(lindexPtr, lindexPtr->GetLargestPossibleRegion());
     ImageType3DU::IndexType lindexIdx;

     unsigned n_comp = 0;

     // printf("fg: %i, bg: %i\n", par.gmm_fg.n_comp, par.gmm_bg.n_comp);
     if (par.gmm_fg.n_comp > par.gmm_bg.n_comp)
	  n_comp = par.gmm_fg.n_comp;
     else
	  n_comp = par.gmm_bg.n_comp;
     
     priors.set_size(par.n_samples, n_comp);
     priors.fill(0);

     // read 4D prior image. (last dim is num of Gaussian components)
     ReaderType4DF::Pointer priorReader = ReaderType4DF::New();
     priorReader->SetFileName(prior_file);
     priorReader->Update();
     ImageType4DF::Pointer priorPtr = priorReader->GetOutput();
     ImageType4DF::SizeType priorSize = priorPtr->GetLargestPossibleRegion().GetSize();
     IteratorType4DF priorIt(priorPtr, priorPtr->GetLargestPossibleRegion());