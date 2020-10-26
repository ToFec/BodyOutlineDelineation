//============================================================================
// Name        : BodyOutlineDelineation.cpp
// Author      : Tobias Fechter
// Version     :
// Copyright   : 
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <fstream>
#include <Utils.h>
#include <ImageWrapper.h>
#include <ImageFactory.h>
#include "Body/BodyClassWrapper.h"
#include <NRRD/nrrd_image.hxx>

using namespace std;

int main(int argc, char ** argv)
{
	if (argc < 3)
	{
		std::cerr << "no enough input parameters!";
		return 1;
	}
	std::string filename = std::string(argv[1]);
	ImageWrapper* imgwrp = ImageFactory::instance()->createImage(filename);
	BodyClassWrapper bodyClassWrp(imgwrp);
	ImageWrapper* outputImg = bodyClassWrp.runSegmentation();

	std::string outputFileName = std::string(argv[2]);
	outputImg->saveImage(outputFileName);

	if (argc > 3)
	{
		double computationTime = bodyClassWrp.getComputationTime();
		std::string metaInfFile = std::string(argv[3]);
		ofstream myfile;
		myfile.open(metaInfFile.c_str());
		myfile << "computationTime:" << computationTime;
		myfile.close();
	}

	return 0;
}
