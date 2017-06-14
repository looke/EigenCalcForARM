/*
 * MatrixTransposer.h
 *
 *  Created on: 2017��2��25��
 *      Author: looke
 */

#ifndef MATRIXTRANSPOSER_H_
#define MATRIXTRANSPOSER_H_
#include "BasicMatrix.h"

class MatrixTransposer
{
public:
	MatrixTransposer(){};

	//ת�÷���
	bool transposeSquareMatrix(BasicMatrix* p_opMatrix);

	//ת���������
	bool transposeMatrix(BasicMatrix* p_opMatrix, BasicMatrix* p_resultMatrix);

	//string getMatrixNotSquareErrorMessage(int rowNumber, int columnNumber);
	//string getMatrixNotSquareErrorMessage(int op_rowNumber, int op_columnNumber,int re_rowNumber, int re_columnNumber);
	//~MatrixTransposer() {};
protected:

};

#endif /* MATRIXTRANSPOSER_H_ */
