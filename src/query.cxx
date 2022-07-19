
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
		  vnl_sparse_matrix <double> & con_map, // can not set const due to a vnl bug.
		  const vnl_vector<unsigned> & alpha,
		  const ParType & par,
		  double smalleta,
		  unsigned n_iter)
{
     unsigned nbr_id = 0;
     vnl_sparse_matrix<double>::row con_map_row;
     vnl_vector<double> score_map_old(score_map);
     vnl_vector<double> exp_term (2, 0);
     unsigned iter = 0;
     double changed_score = 1e7;

     for (iter = 0; iter < n_iter; iter ++) {
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
		    exp_term[0] += smalleta * (1 * score_map[nbr_id] + 0 *(1-score_map[nbr_id]));
		    // the above line does not take into account the
		    // patches. Int assume each patch is just a single voxel. To
		    // account for that, need to follow the build_nliks func.
	       }

	       for (vnl_sparse_matrix<double>::row::const_iterator col_iter = con_map_row.begin(); col_iter != con_map_row.end();  ++ col_iter) {
		    nbr_id = (*col_iter).first;
		    // if score == 0 (BG)
		    exp_term[1] += smalleta * (0 * score_map[nbr_id] + 1 * (1-score_map[nbr_id]));
	       }

	       score_map[n] = exp_term[0] - exp_term[1];
	       score_map[n] = 1 / (1 + exp(- score_map[n]));
	  } // n i.e. sample id.

	  changed_score = (score_map_old - score_map).two_norm() / score_map.two_norm();

	  if (par.verbose >= 1) {
	       printf("logistic_init(): small eta = %f, teration %i, changed score = %f\n", smalleta, iter, changed_score);
	  }
     }
     return 0;
}

void preprocessingBinaryPredicProb(ImageType3DI::Pointer binaryPredictProb)
{
     std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
     std::cout << "Do preprocessing for binary image." << std::endl;

     const ImageType3DI::SizeType sizeOfImage = binaryPredictProb->GetLargestPossibleRegion().GetSize();
     
     int width = sizeOfImage[0];
     int height = sizeOfImage[1];
     int slice = sizeOfImage[2];
     
     ImageType3DI::RegionType outRegion = binaryPredictProb->GetLargestPossibleRegion();

     int i,j,k;

     ImageType3DUC::Pointer binaryPredictProbUC = ImageType3DUC::New();
     binaryPredictProbUC->SetRegions(outRegion);
     binaryPredictProbUC->Allocate();
     binaryPredictProbUC->FillBuffer( 0 );

     ImageType3DUC::IndexType imageIndexUC;
     
     ImageType3DI::PixelType imageValueDI;
     for(k=0; k<slice; k++)
     {
	  for(j=0; j<height; j++)
	  {
	       for(i=0; i<width; i++)
	       {
		    imageValueDI = 0;
		    imageIndexUC[0] = i;
		    imageIndexUC[1] = j;
		    imageIndexUC[2] = k;

		    imageValueDI = binaryPredictProb->GetPixel( imageIndexUC );
		    if(imageValueDI > 0)
		    {
			binaryPredictProbUC->SetPixel(imageIndexUC, 1 );
		    }
	       }
	  }
     }

     // do large erosion

     const unsigned int Dimension = 3;

     typedef itk::BinaryBallStructuringElement<PixelType3DUC, Dimension> StructuringElementType;
     typedef itk::BinaryErodeImageFilter<ImageType3DUC, ImageType3DUC, StructuringElementType> ErodeFilterType;
	  
     ErodeFilterType::Pointer binaryErode = ErodeFilterType::New();

     StructuringElementType structuringElement;

     structuringElement.SetRadius( 3 ); // 3x3 structuring element
     structuringElement.CreateStructuringElement();
     binaryErode->SetKernel( structuringElement );
     binaryErode->SetInput( binaryPredictProbUC );
     binaryErode->SetErodeValue( 1 );

     ImageType3DUC::Pointer tempBinaryPredictErosion;
     tempBinaryPredictErosion = binaryErode->GetOutput();

     // std::string test_save_5 = "test_original_binary_Predict.nii.gz";
     // WriterType3DUC::Pointer writer5 = WriterType3DUC::New();
     // writer5->SetFileName( test_save_5 );
     // writer5->SetInput( binaryPredictProbUC ) ;
     // writer5->Update();

     ///////////////////////////////////////////
     // do dilation (same radius of structure as erosion)
     typedef itk::BinaryDilateImageFilter<ImageType3DUC, ImageType3DUC, StructuringElementType> DilateFilterType;

     DilateFilterType::Pointer binaryDilate = DilateFilterType::New();

     StructuringElementType structuringElementDilation;

     structuringElementDilation.SetRadius( 3 ); // 3x3 structuring element
     structuringElementDilation.CreateStructuringElement();

     binaryDilate->SetKernel( structuringElementDilation );
     binaryDilate->SetInput( tempBinaryPredictErosion );
     binaryDilate->SetDilateValue( 1 );
     ImageType3DUC::Pointer tempImagePointer = binaryDilate->GetOutput();

     // use extract image filter instead of saving the dilated volume - it works! but slow
     /*
     typedef itk::ExtractImageFilter< ImageType3DUC, ImageType3DUC > FilterType;
     FilterType::Pointer filter = FilterType::New();
     filter->SetExtractionRegion(outRegion);
     filter->SetInput(tempImagePointer);
     filter->SetDirectionCollapseToIdentity(); // This is required.
     filter->Update();
     ImageType3DUC::Pointer tempBinaryPredictErosionThenDilation = filter->GetOutput();
     */

     // old way: saving the dilated volume -- it works by slow.
     
     ImageType3DUC::Pointer tempBinaryPredictErosionThenDilation;
     tempBinaryPredictErosionThenDilation = binaryDilate->GetOutput();
     std::string test_save_5 = "test_erosion_and_dilation_binary_Predict.nii.gz";
     WriterType3DUC::Pointer writer5 = WriterType3DUC::New();
     writer5->SetFileName( test_save_5 );
     writer5->SetInput( tempBinaryPredictErosionThenDilation ) ;
     writer5->Update();
     

     // modify the content of binaryPredictProb
     
     ImageType3DUC::PixelType imageValueUC;

     for(k=0; k<slice; k++)
     {
	  for(j=0; j<height; j++)
	  {
	       for(i=0; i<width; i++)
	       {
		    imageValueDI = 0;
		    imageIndexUC[0] = i;
		    imageIndexUC[1] = j;
		    imageIndexUC[2] = k;

		    imageValueUC = tempBinaryPredictErosionThenDilation->GetPixel( imageIndexUC );

		    binaryPredictProb->SetPixel(imageIndexUC, imageValueUC );

	       }
	  }
     }

}

void morphologicalProcessingTrueCandidatesVolume(ImageType3DUC::Pointer trueCandidatesVolume,
						 ImageType3DUC::Pointer trueCandidatesVolumeErosion, 
						 ImageType3DUC::Pointer trueCandidatesVolumeDilationMinusErosion)
{
     std::cout << "Do morphological processing for true candidates volume to get masks." << std::endl;

     const ImageType3DI::SizeType sizeOfImage = trueCandidatesVolume->GetLargestPossibleRegion().GetSize();
     
     int width = sizeOfImage[0];
     int height = sizeOfImage[1];
     int slice = sizeOfImage[2];

     ImageType3DI::RegionType outRegion = trueCandidatesVolume->GetLargestPossibleRegion();

     int i,j,k;

     ImageType3DUC::IndexType imageIndexUC;
     ImageType3DUC::PixelType imageValueUCerosion, imageValueUCdelision;

     const unsigned int Dimension = 3;

     // do erosion = trueCandidatesVolumeErosion
     /*
     typedef itk::BinaryBallStructuringElement<PixelType3DUC, Dimension> StructuringElementType;
     typedef itk::BinaryErodeImageFilter<ImageType3DUC, ImageType3DUC, StructuringElementType> ErodeFilterType;
	  
     ErodeFilterType::Pointer binaryErode = ErodeFilterType::New();

     StructuringElementType structuringElement;

     structuringElement.SetRadius( 1 ); // 3x3 structuring element
     structuringElement.CreateStructuringElement();
     binaryErode->SetKernel( structuringElement );
     binaryErode->SetInput( trueCandidatesVolume );
     binaryErode->SetErodeValue( 1 );

     ImageType3DUC::Pointer tempTrueCandidatesVolumeErosion;
     tempTrueCandidatesVolumeErosion = binaryErode->GetOutput();

     std::string test_save_4 = "test_true_candidates_erosion.nii.gz";
     WriterType3DUC::Pointer writer4 = WriterType3DUC::New();
     writer4->SetFileName( test_save_4 );
     writer4->SetInput( tempTrueCandidatesVolumeErosion ) ;
     writer4->Update();
     */

     // no erosion
     ImageType3DUC::Pointer tempTrueCandidatesVolumeErosion;
     tempTrueCandidatesVolumeErosion = trueCandidatesVolume;

     for(k=0; k<slice; k++)
     {
	  for(j=0; j<height; j++)
	  {
	       for(i=0; i<width; i++)
	       {
		    imageValueUCerosion = 0;
		    imageIndexUC[0] = i;
		    imageIndexUC[1] = j;
		    imageIndexUC[2] = k;
		    imageValueUCerosion = tempTrueCandidatesVolumeErosion->GetPixel( imageIndexUC );

		    trueCandidatesVolumeErosion->SetPixel(imageIndexUC, imageValueUCerosion );
	       }
	  }
     }
     

     // do dilation

     typedef itk::BinaryBallStructuringElement<PixelType3DUC, Dimension> StructuringElementTypeDilation;
     typedef itk::BinaryDilateImageFilter<ImageType3DUC, ImageType3DUC, StructuringElementTypeDilation> DilateFilterType;
     DilateFilterType::Pointer binaryDilate = DilateFilterType::New();

     StructuringElementTypeDilation structuringElementDilation;
     structuringElementDilation.SetRadius( 1 ); // 3x3 structuring element
     structuringElementDilation.CreateStructuringElement();

     binaryDilate->SetKernel( structuringElementDilation );
     binaryDilate->SetInput( trueCandidatesVolume );
     binaryDilate->SetDilateValue( 1 );

     // use extract image filter instead of saving the dilated volume - don't use this!
     /*
     ImageType3DUC::Pointer tempImagePointer = binaryDilate->GetOutput();
     typedef itk::ExtractImageFilter< ImageType3DUC, ImageType3DUC > FilterType;
     FilterType::Pointer filter = FilterType::New();
     filter->SetExtractionRegion(outRegion);
     filter->SetInput(tempImagePointer);
     filter->SetDirectionCollapseToIdentity(); // This is required.
     filter->Update();
     ImageType3DUC::Pointer trueCandidatesVolumeDilation = filter->GetOutput();
     */

     // use this to make sure no segmentation fault and no wrong result.
     ImageType3DUC::Pointer trueCandidatesVolumeDilation = ImageType3DUC::New();
     trueCandidatesVolumeDilation->SetRegions(outRegion);
     trueCandidatesVolumeDilation->Allocate();
     trueCandidatesVolumeDilation->FillBuffer( 0 );
     trueCandidatesVolumeDilation = binaryDilate->GetOutput();
     std::string test_save_5 = "test_true_candidates_dilation.nii.gz";
     WriterType3DUC::Pointer writer5 = WriterType3DUC::New();
     writer5->SetFileName( test_save_5 );
     writer5->SetInput( trueCandidatesVolumeDilation ) ;
     writer5->Update();
     

     // dilation - erosion = trueCandidatesVolumeDilationMinusErosion

     for(k=0; k<slice; k++)
     {
	  for(j=0; j<height; j++)
	  {
	       for(i=0; i<width; i++)
	       {
		    imageValueUCerosion = 0;
		    imageValueUCdelision = 0;
		    imageIndexUC[0] = i;
		    imageIndexUC[1] = j;
		    imageIndexUC[2] = k;
		    imageValueUCerosion = trueCandidatesVolumeErosion->GetPixel( imageIndexUC );
		    imageValueUCdelision = trueCandidatesVolumeDilation->GetPixel( imageIndexUC );

		    if(imageValueUCerosion < 1 && imageValueUCdelision > 0)
		    {
			 trueCandidatesVolumeDilationMinusErosion->SetPixel(imageIndexUC, 1 );
		    }
	       }
	  }
     }

     // std::string test_save_6 = "test_true_candidates_dilationMinusErosion.nii.gz";
     // WriterType3DUC::Pointer writer6 = WriterType3DUC::New();
     // writer6->SetFileName( test_save_6 );
     // writer6->SetInput( trueCandidatesVolumeDilationMinusErosion ) ;
     // writer6->Update();

}


bool computeShapeScoreforCcpCandidates (ImageType3DI::Pointer ccpPredictProb,
					ImageType3DC::Pointer trimap,
					const int& numTopCcptoEvaluate,
					std::vector< std::pair<int, int> >& ccpVectorPairs,
					ConnectedCompInfor allConnectedComponents[],
					const float& thresholdScoreforTopCandidates)


{
     std::cout << "Compute shape information for Top ranked connected componets:" << std::endl;

     const ImageType3DI::SizeType sizeOfImage = ccpPredictProb->GetLargestPossibleRegion().GetSize();
     
     int width = sizeOfImage[0];
     int height = sizeOfImage[1];
     int slice = sizeOfImage[2];

     ImageType3DI::RegionType outRegion = ccpPredictProb->GetLargestPossibleRegion();

     // print the sorted connected components based on volume
     int i,j,k,n;
     int curLabelValue = 0;
     ImageType3DI::IndexType ccpIndex;
     ImageType3DI::PixelType ccpLabelValue;

     const unsigned int Dimension = 3;

     float curScoreCriteria = 0;
     float curShapeScore = 0;
     float shapeScoreScaling = 1.0;
     float curVolumeMultiplyProb = 0;
     
     //mypairFI * const ccpTopPair = (mypairFI*)_alloca(numTopCcptoEvaluate * sizeof(mypairFI)); // got error in linux
     mypairFI * const ccpTopPair = new mypairFI[numTopCcptoEvaluate];

     //for(n=0; n<numTopCcptoEvaluate; n++)
     for(n=0; n<numTopCcptoEvaluate; n++)
     {
	  //cout << "Volume: " << ccpVectorPairs[n].first << " Label: "<<  ccpVectorPairs[n].second << std::endl;
	  cout << " Label: "<<  ccpVectorPairs[n].second << " |  Volume: " << ccpVectorPairs[n].first ;
	  curLabelValue = ccpVectorPairs[n].second;
	  cout << " | Avg Prob.: "<< allConnectedComponents[curLabelValue].avgPredictProb << std::endl ;

	  // for each connected components
	  // step 1: extract the label
	  

	  ImageType3DUC::Pointer curCCPLabel = ImageType3DUC::New();
	  curCCPLabel->SetRegions(outRegion);
	  curCCPLabel->Allocate();
	  curCCPLabel->FillBuffer( 0 );

	  for(k=0; k<slice; k++)
	  {
	       for(j=0; j<height; j++)
	       {
		    for(i=0; i<width; i++)
		    {

			 ccpLabelValue = 0;
			 ccpIndex[0] = i;
			 ccpIndex[1] = j;
			 ccpIndex[2] = k;
			 ccpLabelValue = ccpPredictProb->GetPixel( ccpIndex );
			 
			 if(ccpLabelValue == curLabelValue)
			 {
			      curCCPLabel->SetPixel(ccpIndex, 1 );
			 }
		    }
	       }
	  }