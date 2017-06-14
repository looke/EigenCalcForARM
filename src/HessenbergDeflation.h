/*
 * HessenbergDeflation.h
 *
 *  Created on: 2017��6��1��
 *      Author: looke
 */

#ifndef TRANSFORMATION_BASIC_HESSENBERGDEFLATION_H_
#define TRANSFORMATION_BASIC_HESSENBERGDEFLATION_H_

#include "BasicMatrix.h"
class HessenbergDeflation
{
public:
	HessenbergDeflation();

	bool findNewDeflationPoint(BasicMatrix* p_HessenbergMatrix, int deflationStart, int deflationEnd);

	int getNewDeflationStart();
	int getNewDeflationEnd();
protected:
	int deflateStart_new;
	int deflateEnd_new;
};

#endif /* TRANSFORMATION_BASIC_HESSENBERGDEFLATION_H_ */
