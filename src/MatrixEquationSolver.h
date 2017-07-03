/*
 * MatrixEquationSolver.h
 *
 *  Created on: 2017Äê2ÔÂ19ÈÕ
 *      Author: looke
 */

#ifndef MATRIXEQUATIONSOLVER_H_
#define MATRIXEQUATIONSOLVER_H_

#include "BasicMatrix.h"

class MatrixEquationSolver
{

public:

	MatrixEquationSolver(BasicMatrix* input_opMatrix);

	void init(BasicMatrix* input_opMatrix);
	void reload(BasicMatrix* input_opMatrix);

	bool solveMatrixEquation();
	//double* getRoots();

	//virtual ~MatrixEquationSolver() {};

protected:

	BasicMatrix* p_opMatrix;

	//double *roots;
	int rowNumber;
	int columnNumber;

	void swapRow(int from, int to);
	int findPivotRow(int startRow);
	void normalizePivotRow(int pivotRow);
	void eliminateSubMatrix(int pivotRow);

	void gaussElim_ColmnPrin();

	int coefficientMatrixRank();
	int augmentedMatrixRank();

	bool isAllZeroForArgumentRow(int arguRowNo);
	bool isAllZeroForCoRow(int coRowNo);

	bool calcRoots();

};

#endif /* MATRIXEQUATIONSOLVER_H_ */
