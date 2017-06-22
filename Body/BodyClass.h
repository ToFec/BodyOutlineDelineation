/*
 * TracheaClass.h
 *
 * @Author:  Jose Dolz
 * @  Date:  Jan 16, 2014
 * @E-Mail:  jose.dolz@aquilab.com
 * 
 */

#ifndef __BODYCLASS_H
#define __BODYCLASS_H

#include <itkImage.h>
#include <itkBinaryThresholdImageFilter.h>
#include <itkThresholdImageFilter.h>
#include <itkCastImageFilter.h>
#include <itkRescaleIntensityImageFilter.h>
#include <itkExtractImageFilter.h>
#include <itkBinaryImageToShapeLabelMapFilter.h>
#include <itkLabelMapToBinaryImageFilter.h>
#include <itkConnectedThresholdImageFilter.h>
#include <itkInvertIntensityImageFilter.h>
#include <itkBinaryBallStructuringElement.h>
#include <itkBinaryErodeImageFilter.h>
#include <itkBinaryDilateImageFilter.h>
#include <itkConstantPadImageFilter.h>
#include <itkImageFileWriter.h>
#include <itkCropImageFilter.h>
#include <ctime>
#include <fstream>
#include <iostream>

template<typename ImageType, typename SegmentationImageType>
class BodyClass
{
public:
	BodyClass(){};
	virtual ~BodyClass(){};

	typename SegmentationImageType::Pointer		run_BodyExtraction(ImageType* inputImage);
	typename SegmentationImageType::SizeType		getBodyBoundingBox();

private:
	typename SegmentationImageType::Pointer	extractBody(ImageType* inputImage);
	
	int computeBestThresholdValue(ImageType* inputImage );
	typename SegmentationImageType::Pointer getBody(typename ImageType::Pointer inputImage, int threshold);

	typename SegmentationImageType::Pointer fillHolesByInvertFilter(SegmentationImageType* inputMask);

	typename SegmentationImageType::Pointer cleanBodySegmentedVolume(SegmentationImageType* inputMask);

	typename SegmentationImageType::SizeType	bodySize;

	typename SegmentationImageType::Pointer removeEmptySlices(typename SegmentationImageType::Pointer);

	//SegmentationImageType::Pointer	bodySegmentation;
};

template<typename ImageType, typename SegmentationImageType>
typename SegmentationImageType::Pointer BodyClass<ImageType, SegmentationImageType>::run_BodyExtraction(ImageType* inputImage)
{
		return extractBody(inputImage);
}

template<typename ImageType, typename SegmentationImageType>
typename SegmentationImageType::Pointer BodyClass<ImageType, SegmentationImageType>::extractBody(ImageType* inputImage)
{

	typedef itk::ThresholdImageFilter <ImageType>
	ThresholdImageFilterType;

	time_t tstart, tend;
	tstart = time(0);

	typename ThresholdImageFilterType::Pointer thresholdFilterLow
	= ThresholdImageFilterType::New();
	thresholdFilterLow->SetInput(inputImage);
	if (std::numeric_limits<typename ImageType::PixelType>::lowest() > -1024)
	{
		thresholdFilterLow->ThresholdBelow(std::numeric_limits<typename ImageType::PixelType>::lowest());
		thresholdFilterLow->SetOutsideValue(std::numeric_limits<typename ImageType::PixelType>::lowest());
	}
	else
	{
		thresholdFilterLow->ThresholdBelow(-1024);
		thresholdFilterLow->SetOutsideValue(-1024);
	}



	typename ThresholdImageFilterType::Pointer thresholdFilterHigh
	= ThresholdImageFilterType::New();
	thresholdFilterHigh->SetInput(thresholdFilterLow->GetOutput());

	if (std::numeric_limits<typename ImageType::PixelType>::max() < 1024)
	{
		thresholdFilterHigh->ThresholdAbove(std::numeric_limits<typename ImageType::PixelType>::max());
		thresholdFilterHigh->SetOutsideValue(std::numeric_limits<typename ImageType::PixelType>::max());
	}
	else
	{
		thresholdFilterHigh->ThresholdAbove(1024);
		thresholdFilterHigh->SetOutsideValue(1024);
	}
	thresholdFilterHigh->UpdateLargestPossibleRegion();

	int threshold = computeBestThresholdValue(thresholdFilterHigh->GetOutput());

	typename ImageType::SizeType extendRegion;
	extendRegion[0] = 10;
	extendRegion[1] = 10;
	extendRegion[2] = 0;
	typename ImageType::PixelType constantPixel = -1024;

	typedef itk::ConstantPadImageFilter <ImageType, ImageType> ConstantPadImageFilterType;
	typename ConstantPadImageFilterType::Pointer padFilter = ConstantPadImageFilterType::New();
	padFilter->SetInput(thresholdFilterHigh->GetOutput());
	padFilter->SetPadLowerBound(extendRegion);
	padFilter->SetPadUpperBound(extendRegion);
	padFilter->SetConstant(constantPixel);
	padFilter->Update();

	typename SegmentationImageType::Pointer segmentedBody = getBody(padFilter->GetOutput(), threshold);
	typedef itk::CropImageFilter <SegmentationImageType, SegmentationImageType> CropImageFilterType;
	typename SegmentationImageType::SizeType cropSize;
	cropSize[0] = 10;
	cropSize[1] = 10;
	cropSize[2] = 0;
	typename CropImageFilterType::Pointer cropFilter = CropImageFilterType::New();
	cropFilter->SetInput(segmentedBody);
	cropFilter->SetBoundaryCropSize(cropSize);
	cropFilter->Update();

	tend = time(0);
	std::cout << "It took "<< difftime(tend, tstart) <<" second(s)."<< std::endl;

	return cropFilter->GetOutput();
}

template<typename ImageType, typename SegmentationImageType>
typename SegmentationImageType::SizeType BodyClass<ImageType, SegmentationImageType>::getBodyBoundingBox()
{
	return bodySize;
}

template<typename ImageType, typename SegmentationImageType>
int BodyClass<ImageType, SegmentationImageType>::computeBestThresholdValue(ImageType* inputImage )
{
	typename ImageType::SizeType  sizeMask   = inputImage->GetLargestPossibleRegion().GetSize();
	typename ImageType::IndexType startMask  = inputImage->GetLargestPossibleRegion().GetIndex();

	startMask[2] = int((sizeMask[2] + 0.0)/ 2);
	sizeMask[2] = 0;

	// Create the region object to extract
	typename ImageType::RegionType desiredRegionG;
	desiredRegionG.SetSize(  sizeMask );
	desiredRegionG.SetIndex( startMask );

	typedef itk::Image<short, 2>                                             SliceType;

	typedef itk::ExtractImageFilter< ImageType,
									 SliceType >			    			 SliceExtractorTypeSS;

	typename SliceExtractorTypeSS::Pointer slicegraySS		      = SliceExtractorTypeSS::New();

	slicegraySS->SetExtractionRegion( desiredRegionG );
	slicegraySS->SetInput( inputImage );

#if ITK_VERSION_MAJOR >= 4
	slicegraySS->SetDirectionCollapseToIdentity();   // Mandatory in ITK Versions >= 4
#endif
	slicegraySS->Update();

	// Run over the slice to check the suitable value for the threshold

	SliceType::IndexType index;
	index[0]	=	index[1]	=	0;

	int min = slicegraySS->GetOutput()->GetPixel(index);
	int max = slicegraySS->GetOutput()->GetPixel(index);
	int diff = 0;
	int prevValue = slicegraySS->GetOutput()->GetPixel(index);

	int sizeDiag	=	slicegraySS->GetOutput()->GetLargestPossibleRegion().GetSize()[0];

	sizeDiag	=	(sizeDiag + 0.0) / 2;

	for(int idDiag = 1; idDiag < sizeDiag ; idDiag++)
	{
		SliceType::IndexType index;
		index.Fill(idDiag);
		int currentValue = slicegraySS->GetOutput()->GetPixel(index);

		if( abs(currentValue - prevValue) > diff)
		{
			min	=	prevValue;
			max	=	currentValue;
			diff	=	abs(currentValue - prevValue);
		}
		prevValue	=	currentValue;
	}

	int temp;
	if(min > max)
	{
		temp	=	max;
		max		=	min;
		min		=	temp;

	}
	int threshold = max - int(diff * 0.5);

	return threshold;
}

template<typename ImageType, typename SegmentationImageType>
typename SegmentationImageType::Pointer BodyClass<ImageType, SegmentationImageType>::getBody(typename ImageType::Pointer inputImage, int threshold)
{
	// 1 - Apply threshold
	typedef itk::BinaryThresholdImageFilter <ImageType, SegmentationImageType>
	BinaryThresholdImageFilterType;

	typename BinaryThresholdImageFilterType::Pointer thresholdFilter = BinaryThresholdImageFilterType::New();
	thresholdFilter->SetInput(inputImage);
	thresholdFilter->SetInsideValue(255);
	thresholdFilter->SetOutsideValue(0);
	thresholdFilter->SetLowerThreshold(threshold);
	thresholdFilter->Update();

	// 2 - Keep only the biggest label
	typedef itk::BinaryImageToShapeLabelMapFilter<SegmentationImageType> BinaryImageToShapeLabelMapFilterType;
	typename BinaryImageToShapeLabelMapFilterType::Pointer binaryImageToShapeLabelMapFilter = BinaryImageToShapeLabelMapFilterType::New();
	binaryImageToShapeLabelMapFilter->SetInput(thresholdFilter->GetOutput());
	binaryImageToShapeLabelMapFilter->Update();

	int numLabels			=	binaryImageToShapeLabelMapFilter->GetOutput()->GetNumberOfLabelObjects();
	int largestLabelID		=	-1;
	int largestLabelSize	=	-1;
	std::vector<int>	labelsArray;
	// Find the largest Label
	for(int i = 0; i < numLabels; i++)
	{
		typename BinaryImageToShapeLabelMapFilterType::OutputImageType::LabelObjectType* labelObject = binaryImageToShapeLabelMapFilter->GetOutput()->GetNthLabelObject(i);
		 int numPixels =  labelObject->GetNumberOfPixels();
		 labelsArray.push_back(labelObject->GetLabel());
		 if(numPixels > largestLabelSize)
		 {
			 largestLabelSize	=	numPixels;
			 largestLabelID		=	labelObject->GetLabel();
		 }
	}

	// Remove all the other Labels
	std::vector<int>::iterator iterLabel  =  labelsArray.begin();

	for(; iterLabel != labelsArray.end() ; iterLabel++)
	{
		 if((*iterLabel) != largestLabelID)
		 {
			binaryImageToShapeLabelMapFilter->GetOutput()->RemoveLabel((*iterLabel));
		 }
	}

	typedef itk::LabelMapToLabelImageFilter<typename BinaryImageToShapeLabelMapFilterType::OutputImageType, SegmentationImageType> LabelMapToLabelImageFilterType;
	typename LabelMapToLabelImageFilterType::Pointer labelMapToLabelImageFilter = LabelMapToLabelImageFilterType::New();
	labelMapToLabelImageFilter->SetInput(binaryImageToShapeLabelMapFilter->GetOutput());
	labelMapToLabelImageFilter->Update();

	// This image has still holes
	typename SegmentationImageType::Pointer mask = SegmentationImageType::New();
	mask = fillHolesByInvertFilter(labelMapToLabelImageFilter->GetOutput());
	mask = cleanBodySegmentedVolume(mask);
	return mask;
}


template<typename ImageType, typename SegmentationImageType>
typename SegmentationImageType::Pointer BodyClass<ImageType, SegmentationImageType>::cleanBodySegmentedVolume(SegmentationImageType* inputMask)
{
	// Clean the resulted image of noise, like the treatment table
	typedef itk::BinaryBallStructuringElement<unsigned char,SegmentationImageType::ImageDimension> StructuringElementType;
	StructuringElementType structuringElement;
	structuringElement.SetRadius(2);
	structuringElement.CreateStructuringElement();

	typedef itk::BinaryErodeImageFilter <SegmentationImageType, SegmentationImageType, StructuringElementType>
         BinaryErodeImageFilterType;
	typename BinaryErodeImageFilterType::Pointer erodeFilter
          = BinaryErodeImageFilterType::New();
	erodeFilter->SetInput(inputMask);
	erodeFilter->SetKernel(structuringElement);
	erodeFilter->Update();

	typedef itk::BinaryImageToShapeLabelMapFilter<SegmentationImageType> BinaryImageToShapeLabelMapFilterType;
	typename BinaryImageToShapeLabelMapFilterType::Pointer binaryImageToShapeLabelMapFilterCheck = BinaryImageToShapeLabelMapFilterType::New();
	binaryImageToShapeLabelMapFilterCheck->SetInput(erodeFilter->GetOutput());
	binaryImageToShapeLabelMapFilterCheck->Update();

	int numLabels	=	binaryImageToShapeLabelMapFilterCheck->GetOutput()->GetNumberOfLabelObjects();

	std::vector<int>	labelsArrayToRemove;
	int largestLabelSize	=	0;
	int largestLabelID		=	0;
	// Find the largest Label
	for(int i = 0; i < numLabels; i++)
	{
		typename BinaryImageToShapeLabelMapFilterType::OutputImageType::LabelObjectType* labelObject = binaryImageToShapeLabelMapFilterCheck->GetOutput()->GetNthLabelObject(i);
		 int numPixels =  labelObject->GetNumberOfPixels();
		 labelsArrayToRemove.push_back(labelObject->GetLabel());
		 if(numPixels > largestLabelSize)
		 {
			 largestLabelSize	=	numPixels;
			 largestLabelID		=	labelObject->GetLabel();
			 bodySize			=	labelObject->GetBoundingBox().GetSize();
		 }
	}

	// Remove all the other Labels
	std::vector<int>::iterator iterLabelClean  =  labelsArrayToRemove.begin();

	for(; iterLabelClean != labelsArrayToRemove.end() ; iterLabelClean++)
	{
		 if((*iterLabelClean) != largestLabelID)
		 {
			binaryImageToShapeLabelMapFilterCheck->GetOutput()->RemoveLabel((*iterLabelClean));
		 }
	}

	typedef itk::LabelMapToBinaryImageFilter<typename BinaryImageToShapeLabelMapFilterType::OutputImageType, SegmentationImageType>LabelMapToLabelImageFilterType;
	typename LabelMapToLabelImageFilterType::Pointer labelMapToLabelImageFilterMask = LabelMapToLabelImageFilterType::New();

	labelMapToLabelImageFilterMask->SetInput(binaryImageToShapeLabelMapFilterCheck->GetOutput());
	labelMapToLabelImageFilterMask->Update();

	// Clean now in 2D
	typedef  itk::Image<unsigned char, 2>			UCharImageType;
	typedef itk::ExtractImageFilter< SegmentationImageType,
									 UCharImageType >			 BinSliceTypeGray;

	typename BinSliceTypeGray::Pointer BinarySlice		    = BinSliceTypeGray::New();
	typename SegmentationImageType::RegionType regionToExtractBinarySlice;
	typename SegmentationImageType::SizeType   sizeBinarySlice ;
	typename SegmentationImageType::IndexType  startBinarySlice;

	int numSlices = labelMapToLabelImageFilterMask->GetOutput()->GetLargestPossibleRegion().GetSize()[2];
	startBinarySlice.Fill(0);
	startBinarySlice[0] = labelMapToLabelImageFilterMask->GetOutput()->GetLargestPossibleRegion().GetIndex()[0];
	startBinarySlice[1] = labelMapToLabelImageFilterMask->GetOutput()->GetLargestPossibleRegion().GetIndex()[1];

	sizeBinarySlice[0]     =  labelMapToLabelImageFilterMask->GetOutput()->GetLargestPossibleRegion().GetSize()[0];
	sizeBinarySlice[1]     =  labelMapToLabelImageFilterMask->GetOutput()->GetLargestPossibleRegion().GetSize()[1];
	sizeBinarySlice[2]  =  0;

	for(int sliceID = 0  ; sliceID < numSlices; sliceID++)
	{
		// Extract Slice
		startBinarySlice[2] =  sliceID;
		regionToExtractBinarySlice.SetIndex(startBinarySlice);
		regionToExtractBinarySlice.SetSize( sizeBinarySlice );
		BinarySlice->SetExtractionRegion( regionToExtractBinarySlice );
		BinarySlice->SetInput( labelMapToLabelImageFilterMask->GetOutput());

		#if ITK_VERSION_MAJOR >= 4
			BinarySlice->SetDirectionCollapseToIdentity();   // Mandatory in ITK Versions >= 4
		#endif

		BinarySlice->Update();

		typedef itk::RescaleIntensityImageFilter<UCharImageType,UCharImageType> RescaleIntensityImageFilterSliceType;
		RescaleIntensityImageFilterSliceType::Pointer bodySliceRescaler	=	RescaleIntensityImageFilterSliceType::New();
		bodySliceRescaler->SetInput(BinarySlice->GetOutput());
		bodySliceRescaler->SetOutputMaximum(255);
		bodySliceRescaler->SetOutputMinimum(0);

		typedef itk::BinaryBallStructuringElement<unsigned char,2> StructuringElementType2D;
		StructuringElementType2D structuringElement2D;
		structuringElement2D.SetRadius(3);
		structuringElement2D.CreateStructuringElement();

		typedef itk::BinaryErodeImageFilter <UCharImageType, UCharImageType, StructuringElementType2D>
			 BinaryErodeImageFilterType2D;
		BinaryErodeImageFilterType2D::Pointer erodeFilter2D
			  = BinaryErodeImageFilterType2D::New();
		erodeFilter2D->SetInput(bodySliceRescaler->GetOutput());
		erodeFilter2D->SetKernel(structuringElement2D);
		erodeFilter2D->Update();

		// Convert to label  in order to see whether some label exist
		typedef itk::BinaryImageToShapeLabelMapFilter<UCharImageType> BinaryImageToShapeLabelMapFilterType2D;
		BinaryImageToShapeLabelMapFilterType2D::Pointer binaryImageToShapeLabelMapFilter2D = BinaryImageToShapeLabelMapFilterType2D::New();
		binaryImageToShapeLabelMapFilter2D->SetInput(erodeFilter2D->GetOutput());
		binaryImageToShapeLabelMapFilter2D->Update();

		int num_Labels	=	binaryImageToShapeLabelMapFilter2D->GetOutput()->GetNumberOfLabelObjects();

		std::vector<int>	labelsArrayToRemove2D;
		for (int labelID = 0; labelID < num_Labels; labelID++)
		{
			BinaryImageToShapeLabelMapFilterType2D::OutputImageType::LabelObjectType* labelObject = binaryImageToShapeLabelMapFilter2D->GetOutput()->GetNthLabelObject(labelID);
			if( labelObject->GetNumberOfPixels() < 2500)
			{
				labelsArrayToRemove2D.push_back(labelObject->GetLabel());
			}
		}

		// Remove all the other Labels
		std::vector<int>::iterator iterLabelClean2D  =  labelsArrayToRemove2D.begin();
		for(; iterLabelClean2D != labelsArrayToRemove2D.end() ; iterLabelClean2D++)
		{
				binaryImageToShapeLabelMapFilter2D->GetOutput()->RemoveLabel((*iterLabelClean2D));
		}

		typedef itk::LabelMapToBinaryImageFilter<BinaryImageToShapeLabelMapFilterType2D::OutputImageType, UCharImageType>LabelMapToLabelImageFilterType;
		typename LabelMapToLabelImageFilterType::Pointer labelMapToLabelImageFilterSliceBack = LabelMapToLabelImageFilterType::New();

		labelMapToLabelImageFilterSliceBack->SetInput(binaryImageToShapeLabelMapFilter2D->GetOutput());
		labelMapToLabelImageFilterSliceBack->SetForegroundValue(255);
		labelMapToLabelImageFilterSliceBack->Update();

		// Dilate
		typedef itk::BinaryDilateImageFilter <UCharImageType, UCharImageType, StructuringElementType2D>
			 BinaryDilateImageFilterType2D;
		typename BinaryDilateImageFilterType2D::Pointer dilateFilter2D
			  = BinaryDilateImageFilterType2D::New();
		dilateFilter2D->SetInput(labelMapToLabelImageFilterSliceBack->GetOutput());
		dilateFilter2D->SetKernel(structuringElement2D);
		dilateFilter2D->Update();

		// Set the slice again
		typename SegmentationImageType::RegionType regionToPutBackBinarySlice;
		typename SegmentationImageType::SizeType   sizeBinaryBackSlice ;
		typename SegmentationImageType::IndexType  startBinaryBackSlice;

		startBinaryBackSlice.Fill(0);
		startBinaryBackSlice[0] = startBinarySlice[0];
		startBinaryBackSlice[1] = startBinarySlice[1];
		startBinaryBackSlice[2]  =  sliceID;
		sizeBinaryBackSlice[0]   =  sizeBinarySlice[0];
		sizeBinaryBackSlice[1]   =  sizeBinarySlice[1];
		sizeBinaryBackSlice[2]   =  1;

		regionToPutBackBinarySlice.SetIndex(startBinaryBackSlice);
		regionToPutBackBinarySlice.SetSize(sizeBinaryBackSlice);

		itk::ImageRegionIterator<SegmentationImageType> volumeIterator( labelMapToLabelImageFilterMask->GetOutput(), regionToPutBackBinarySlice);
		UCharImageType::RegionType regionSlice;
		UCharImageType::SizeType   sizeSlice  = BinarySlice->GetOutput()->GetLargestPossibleRegion().GetSize();
		UCharImageType::IndexType  startSlice = BinarySlice->GetOutput()->GetLargestPossibleRegion().GetIndex();

		regionSlice.SetIndex(startSlice);
		regionSlice.SetSize(sizeSlice);
		itk::ImageRegionIterator<UCharImageType>   sliceIterator(  dilateFilter2D->GetOutput(), regionSlice);

		while(!volumeIterator.IsAtEnd() && !sliceIterator.IsAtEnd())
		{
			volumeIterator.Set(sliceIterator.Get());

			++volumeIterator;
			++sliceIterator;
		}
	}

	typedef itk::BinaryDilateImageFilter <SegmentationImageType, SegmentationImageType, StructuringElementType>         BinaryDilateImageFilterType;
	typename BinaryDilateImageFilterType::Pointer dilateFilter	=	BinaryDilateImageFilterType::New();
	dilateFilter->SetInput(labelMapToLabelImageFilterMask->GetOutput());
	dilateFilter->SetKernel(structuringElement);
	dilateFilter->Update();

	labelMapToLabelImageFilterMask->Update();
	return dilateFilter->GetOutput();
}

template<typename ImageType, typename SegmentationImageType>
typename SegmentationImageType::Pointer BodyClass<ImageType, SegmentationImageType>::fillHolesByInvertFilter(SegmentationImageType* inputMask)
{
	typedef itk::RescaleIntensityImageFilter< SegmentationImageType, SegmentationImageType > RescaleFilterType;
	typename RescaleFilterType::Pointer rescaleFilter = RescaleFilterType::New();
	rescaleFilter->SetInput(inputMask);
	rescaleFilter->SetOutputMinimum(0);
	rescaleFilter->SetOutputMaximum(255);
	rescaleFilter->Update();

	typedef itk::InvertIntensityImageFilter <SegmentationImageType>
	InvertIntensityImageFilterType;
	typename InvertIntensityImageFilterType::Pointer invertIntensityFilter = InvertIntensityImageFilterType::New();
	invertIntensityFilter->SetInput(rescaleFilter->GetOutput());
	invertIntensityFilter->SetMaximum(255);
	invertIntensityFilter->Update();

	typedef itk::BinaryImageToShapeLabelMapFilter<SegmentationImageType> BinaryImageToShapeLabelMapFilterType;
	typename BinaryImageToShapeLabelMapFilterType::Pointer binaryImageToShapeLabelMapFilterMask = BinaryImageToShapeLabelMapFilterType::New();
	binaryImageToShapeLabelMapFilterMask->SetInput(invertIntensityFilter->GetOutput());
	binaryImageToShapeLabelMapFilterMask->Update();

	int numLabels			=	binaryImageToShapeLabelMapFilterMask->GetOutput()->GetNumberOfLabelObjects();
	std::cout <<" numLabels "  << numLabels <<  std::endl;
	int largestLabelID		=	-1;
	int largestLabelSize	=	-1;
	std::vector<int>	labelsArray;
	// Find the largest Label
	for(int i = 0; i < numLabels; i++)
	{
		typename BinaryImageToShapeLabelMapFilterType::OutputImageType::LabelObjectType* labelObject = binaryImageToShapeLabelMapFilterMask->GetOutput()->GetNthLabelObject(i);
		 int numPixels =  labelObject->GetNumberOfPixels();
		 labelsArray.push_back(labelObject->GetLabel());
		 if(numPixels > largestLabelSize)
		 {
			 largestLabelSize	=	numPixels;
			 largestLabelID		=	labelObject->GetLabel();
		 }
	}

	// Remove all the other Labels
	std::vector<int>::iterator iterLabel  =  labelsArray.begin();

	for(; iterLabel != labelsArray.end() ; iterLabel++)
	{
		 if((*iterLabel) != largestLabelID)
		 {
			binaryImageToShapeLabelMapFilterMask->GetOutput()->RemoveLabel((*iterLabel));
		 }
	}

	typedef itk::LabelMapToBinaryImageFilter<typename BinaryImageToShapeLabelMapFilterType::OutputImageType, SegmentationImageType>LabelMapToLabelImageFilterType;
	typename LabelMapToLabelImageFilterType::Pointer labelMapToLabelImageFilterMask = LabelMapToLabelImageFilterType::New();

	labelMapToLabelImageFilterMask->SetInput(binaryImageToShapeLabelMapFilterMask->GetOutput());
	labelMapToLabelImageFilterMask->Update();

	invertIntensityFilter->SetInput(labelMapToLabelImageFilterMask->GetOutput());
	invertIntensityFilter->SetMaximum(255);
	invertIntensityFilter->Update();

	return invertIntensityFilter->GetOutput();
}

#endif  // __TRACHEACLASS_H
