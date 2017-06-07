/*
 * StaticMatrix.h
 *
 *  Created on: 2017年4月25日
 *      Author: looke
 */

#ifndef MATRIX_STATIC_STATICMATRIX_H_
#define MATRIX_STATIC_STATICMATRIX_H_

#include "..\include\matrix\basic\BasicMatrix.h"

class StaticMatrix:public BasicMatrix
{
public:
	StaticMatrix();
	StaticMatrix(int inputRowNum, int inputColumnNum);

	//int rowNum;
	//int columnNum;

	void setMatrixElement(int rNum, int cNum, double val);
	double getMatrixElement(int rNum, int cNum);
	//获取矩阵指定元素的值(对元素值进行整形，小于精度的值直接返回0)
	double getMatrixElementRegulared(int rNum, int cNum, double lowEdge);

	void printMatrix();

	//交换行
	void swapRow(int from, int to);

	//交换列
	void swapColumn(int from, int to);

	//交换对角线主元
	void swapDiagElement(int from, int to);

	//将矩阵重置为单位矩阵
	void resetMatrixToI();

	//将矩阵重置为单位矩阵
	void resetMatrixToZero();

	//获取指定列向量
	void getColumnVector(int columnNo, BasicVector* p_Vector);

	//获取指定行向量
	void getRowVector(int rowNo, BasicVector* p_Vector);

	//获取指定列向量
	void getSubMatrixColumnVector(int subMatrixIndex, int columnNo, BasicVector* p_Vector);

	//获取指定行向量
	void getSubMatrixRowVector(int subMatrixIndex, int rowNo, BasicVector* p_Vector);

	//获取指定对角子矩阵hessenberg列向量
	void getSubMatrixHessenColumnVector(int subMatrixIndex, BasicVector* p_Vector);

private:

	double matrixNxN[20][20];
	//StaticVector columnVector;
	//StaticVector rowVector;
	void initMatrix();

};

#endif /* MATRIX_STATIC_STATICMATRIX_H_ */
