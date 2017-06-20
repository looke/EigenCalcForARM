/*
 * MatrixInverser.h
 *
 *  Created on: 2017Äê2ÔÂ25ÈÕ
 *      Author: looke
 */

#ifndef MATRIXINVERSER_H_
#define MATRIXINVERSER_H_
#include "BasicMatrix.h"

class MatrixInverser
{
public:

	MatrixInverser(BasicMatrix* p_input_opMatrix,BasicMatrix* p_input_inverseMatrix);
	void reload(BasicMatrix* p_input_opMatrix,BasicMatrix* p_input_inverseMatrix);
	void init(BasicMatrix* p_input_opMatrix,BasicMatrix* p_input_inverseMatrix);

	void unZeroPivotRow(int pivotRow);
	void normalizePivotRow(int pivotRow);
	void eliminateSubMatrix(int pivotRow);
	void eliminateUpperMatrix(int startRow);
	void generateInverseMatrix();


	bool isOpMatrixFullRank();
	//BasicMatrix* getInverseMatrix();


	//void printOpMatrix();
	//void printInverseMatrix();
	//virtual ~MatrixInverser() {};
protected:

	BasicMatrix* p_operateMatrix;
	BasicMatrix* p_inverseMatrix;


	bool isSquareMatrix;
	int rowNumber;
	int columnNumber;
	bool isFullRank;
};

#endif /* MATRIXINVERSER_H_ */
