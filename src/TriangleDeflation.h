/*
 * TriangleDeflation.h
 *
 *  Created on: 2017��6��22��
 *      Author: looke
 */

#ifndef TRIANGLEDEFLATION_H_
#define TRIANGLEDEFLATION_H_

#include "BasicMatrix.h"
class TriangleDeflation
{
public:
	TriangleDeflation();

	bool findNewDeflationPoint(BasicMatrix* p_TriangleMatrix, int deflationStart, int deflationEnd);

	int getNewDeflationStart();
	int getNewDeflationEnd();
protected:
	int deflateStart_new;
	int deflateEnd_new;
};



#endif /* TRIANGLEDEFLATION_H_ */
