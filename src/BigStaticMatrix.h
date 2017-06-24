/*
 * BigStaticMatrix.h
 *
 *  Created on: 2017年6月24日
 *      Author: looke
 */

#ifndef BIGSTATICMATRIX_H_
#define BIGSTATICMATRIX_H_
#include "BasicMatrix.h"
#include "BigStaticVector.h"

class BigStaticMatrix : public BasicMatrix
{
public:
	BigStaticMatrix();
	BigStaticMatrix(int inputRowNum, int inputColumnNum);

	//int rowNum;
	//int columnNum;

	bool setMatrixElement(int rNum, int cNum, double val);
	double getMatrixElement(int rNum, int cNum);


	void printMatrix();

	//交换行
	bool swapRow(int from, int to);

	//交换列
	bool swapColumn(int from, int to);

	//交换对角线主元
	bool swapDiagElement(int from, int to);

	//将矩阵重置为单位矩阵
	bool resetMatrixToI();

	//将矩阵重置为单位矩阵
	void resetMatrixToZero();

	//获取指定列向量
	BasicVector* getColumnVector(int columnNo);

	//获取指定行向量
	BasicVector* getRowVector(int rowNo);

	//获取指定列向量
	BasicVector* getSubMatrixColumnVector(int subMatrixIndex, int columnNo);

	//获取指定行向量
	BasicVector* getSubMatrixRowVector(int subMatrixIndex, int rowNo);

	//获取指定对角子矩阵hessenberg列向量
	BasicVector* getSubMatrixHessenColumnVector(int subMatrixIndex);

protected:

	double matrixNxN[20][20];
	BigStaticVector matrixVector;
	void initMatrix();

	//获取矩阵指定元素的值(对元素值进行整形，小于精度的值直接返回0)
	double getMatrixElementRegulared(int rNum, int cNum, double lowEdge);
};

#endif /* BIGSTATICMATRIX_H_ */
