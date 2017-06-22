/*
 * BodyClassWrapper.cpp
 *
 *  Created on: 18.05.2017
 *      Author: fechter
 */

#include "BodyClassWrapper.h"

BodyClassWrapper::BodyClassWrapper(ImageWrapper* imgwrp)
{
	this->imgwrp = imgwrp;
}

BodyClassWrapper::~BodyClassWrapper()
{
	// TODO Auto-generated destructor stub
}

ImageWrapper* BodyClassWrapper::runSegmentation()
{
	int dimension = imgwrp->getDimensions();
	std::string type = imgwrp->getType();

	ImageWrapper* outputImg = 0;
	if (dimension == 2)
	{
		typedef itk::Image<unsigned char, 2> SegmentationImageType;
		typename SegmentationImageType::Pointer algoOutput;
		if (type == "short")
		{
			typedef itk::Image<short, 2> ImageType;
			BodyClass<ImageType, SegmentationImageType> bc;
			Image<ImageType>* inputImg = (Image<ImageType>*)imgwrp;
			algoOutput = bc.run_BodyExtraction(inputImg->getItkImage());
		}
		else if (type == "unsigned short")
		{
			typedef itk::Image<unsigned short, 2> ImageType;
			BodyClass<ImageType, SegmentationImageType> bc;
			Image<ImageType>* inputImg = (Image<ImageType>*)imgwrp;
			algoOutput = bc.run_BodyExtraction(inputImg->getItkImage());
		}
		else if (type == "int")
		{
			typedef itk::Image<int, 2> ImageType;
			BodyClass<ImageType, SegmentationImageType> bc;
			Image<ImageType>* inputImg = (Image<ImageType>*)imgwrp;
			algoOutput = bc.run_BodyExtraction(inputImg->getItkImage());
		}
		else if (type == "unsigned int")
		{
			typedef itk::Image<unsigned int, 2> ImageType;
			BodyClass<ImageType, SegmentationImageType> bc;
			Image<ImageType>* inputImg = (Image<ImageType>*)imgwrp;
			algoOutput = bc.run_BodyExtraction(inputImg->getItkImage());
		}
		else if (type == "float")
		{
			typedef itk::Image<float, 2> ImageType;
			BodyClass<ImageType, SegmentationImageType> bc;
			Image<ImageType>* inputImg = (Image<ImageType>*)imgwrp;
			algoOutput = bc.run_BodyExtraction(inputImg->getItkImage());
		}
		else if (type == "double")
		{
			typedef itk::Image<double, 2> ImageType;
			BodyClass<ImageType, SegmentationImageType> bc;
			Image<ImageType>* inputImg = (Image<ImageType>*)imgwrp;
			algoOutput = bc.run_BodyExtraction(inputImg->getItkImage());
		}
		else if (type == "char")
		{
			typedef itk::Image<char, 2> ImageType;
			BodyClass<ImageType, SegmentationImageType> bc;
			Image<ImageType>* inputImg = (Image<ImageType>*)imgwrp;
			algoOutput = bc.run_BodyExtraction(inputImg->getItkImage());
		}
		else if (type == "unsigned char")
		{
			typedef itk::Image<unsigned char, 2> ImageType;
			BodyClass<ImageType, SegmentationImageType> bc;
			Image<ImageType>* inputImg = (Image<ImageType>*)imgwrp;
			algoOutput = bc.run_BodyExtraction(inputImg->getItkImage());
		}
		outputImg = new Image<SegmentationImageType>(algoOutput);
	}
	else if (dimension == 3)
	{

		typedef itk::Image<unsigned char, 3> SegmentationImageType;
		typename SegmentationImageType::Pointer algoOutput;
		if (type == "short")
		{
			typedef itk::Image<short, 3> ImageType;
			BodyClass<ImageType, SegmentationImageType> bc;
			Image<ImageType>* inputImg = (Image<ImageType>*)imgwrp;
			algoOutput = bc.run_BodyExtraction(inputImg->getItkImage());
		}
		else if (type == "unsigned short")
		{
			typedef itk::Image<unsigned short, 3> ImageType;
			BodyClass<ImageType, SegmentationImageType> bc;
			Image<ImageType>* inputImg = (Image<ImageType>*)imgwrp;
			algoOutput = bc.run_BodyExtraction(inputImg->getItkImage());
		}
		else if (type == "int")
		{
			typedef itk::Image<int, 3> ImageType;
			BodyClass<ImageType, SegmentationImageType> bc;
			Image<ImageType>* inputImg = (Image<ImageType>*)imgwrp;
			algoOutput = bc.run_BodyExtraction(inputImg->getItkImage());
		}
		else if (type == "unsigned int")
		{
			typedef itk::Image<unsigned int, 3> ImageType;
			BodyClass<ImageType, SegmentationImageType> bc;
			Image<ImageType>* inputImg = (Image<ImageType>*)imgwrp;
			algoOutput = bc.run_BodyExtraction(inputImg->getItkImage());
		}
		else if (type == "float")
		{
			typedef itk::Image<float, 3> ImageType;
			BodyClass<ImageType, SegmentationImageType> bc;
			Image<ImageType>* inputImg = (Image<ImageType>*)imgwrp;
			algoOutput = bc.run_BodyExtraction(inputImg->getItkImage());
		}
		else if (type == "double")
		{
			typedef itk::Image<double, 3> ImageType;
			BodyClass<ImageType, SegmentationImageType> bc;
			Image<ImageType>* inputImg = (Image<ImageType>*)imgwrp;
			algoOutput = bc.run_BodyExtraction(inputImg->getItkImage());
		}
		else if (type == "char")
		{
			typedef itk::Image<char, 3> ImageType;
			BodyClass<ImageType, SegmentationImageType> bc;
			Image<ImageType>* inputImg = (Image<ImageType>*)imgwrp;
			algoOutput = bc.run_BodyExtraction(inputImg->getItkImage());
		}
		else if (type == "unsigned char")
		{
			typedef itk::Image<unsigned char, 3> ImageType;
			BodyClass<ImageType, SegmentationImageType> bc;
			Image<ImageType>* inputImg = (Image<ImageType>*)imgwrp;
			algoOutput = bc.run_BodyExtraction(inputImg->getItkImage());
		}
		outputImg = new Image<SegmentationImageType>(algoOutput);
	}
	return outputImg;
}
