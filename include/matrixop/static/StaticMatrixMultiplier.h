/*
 * StaticMatrixMultiplier.h
 *
 *  Created on: 2017Äê2ÔÂ25ÈÕ
 *      Author: looke
 */

#ifndef STATICMATRIXMULTIPLIER_H_
#define STATICMATRIXMULTIPLIER_H_
#include "..\include\matrixop\basic\MatrixMultiplier.h"

class StaticMatrixMultiplier:public MatrixMultiplier
{
public:
	StaticMatrixMultiplier(BasicMatrix* leftOp, BasicMatrix* rightOp, BasicMatrix* resultOp);

	void reload(BasicMatrix* leftOp, BasicMatrix* rightOp, BasicMatrix* resultOp);
	void init(BasicMatrix* leftOp, BasicMatrix* rightOp, BasicMatrix* resultOp) throw(length_error);
	~StaticMatrixMultiplier() {};

protected:

};

#endif /* STATICMATRIXMULTIPLIER_H_ */
