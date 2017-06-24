/*
 * TriangleDeflation.h
 *
 *  Created on: 2017Äê6ÔÂ22ÈÕ
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
