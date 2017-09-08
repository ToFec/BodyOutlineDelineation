/*
 * BodyClassWrapper.h
 *
 *  Created on: 18.05.2017
 *      Author: fechter
 */

#ifndef BODY_BODYCLASSWRAPPER_H_
#define BODY_BODYCLASSWRAPPER_H_

#include "BodyClass.h"
#include <Image.h>
#include <ImageWrapper.h>

class BodyClassWrapper
{
public:
	BodyClassWrapper(ImageWrapper* imgwrp);
	virtual ~BodyClassWrapper();

	ImageWrapper* runSegmentation();

	double getComputationTime();

private:
	ImageWrapper* imgwrp = 0;
	double computationTime = 0.0;
};

#endif /* BODY_BODYCLASSWRAPPER_H_ */
