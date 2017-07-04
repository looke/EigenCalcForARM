/*
 * MatrixMultiplier.h
 *
 *  Created on: 2017年2月25日
 *      Author: looke
 */

#ifndef MATRIXMULTIPLIER_H_
#define MATRIXMULTIPLIER_H_
#include "BasicMatrix.h"

class MatrixMultiplier
{
public:

	MatrixMultiplier(BasicMatrix* leftOp, BasicMatrix* rightOp, BasicMatrix* resultOp);

	void reload(BasicMatrix* leftOp, BasicMatrix* rightOp, BasicMatrix* resultOp);
	void init(BasicMatrix* leftOp, BasicMatrix* rightOp, BasicMatrix* resultOp);

	bool multiplyCalc();

	//void printMultiplyResult();

	//string getMatrixLengthErrorMessage(BasicMatrix* leftOp, BasicMatrix* rightOp, BasicMatrix* resultOp);

	//BasicMatrix* getMultiplyResult();
	virtual ~MatrixMultiplier() {};
protected:

	BasicMatrix* p_leftOpMatrix;	//左操作数矩阵
	BasicMatrix* p_rightOpMatrix;//右操作数矩阵

	BasicMatrix* p_MultiResult;
};

#endif /* MATRIXMULTIPLIER_H_ */
