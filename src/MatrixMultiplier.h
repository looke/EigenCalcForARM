/*
 * MatrixMultiplier.h
 *
 *  Created on: 2017��2��25��
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

	BasicMatrix* p_leftOpMatrix;	//�����������
	BasicMatrix* p_rightOpMatrix;//�Ҳ���������

	BasicMatrix* p_MultiResult;
};

#endif /* MATRIXMULTIPLIER_H_ */
