/*
 * StaticMatrix.h
 *
 *  Created on: 2017��4��25��
 *      Author: looke
 */

#ifndef MATRIX_STATIC_STATICMATRIX_H_
#define MATRIX_STATIC_STATICMATRIX_H_

#include "BasicMatrix.h"
#include "StaticVector.h"

class StaticMatrix:public BasicMatrix
{
public:
	StaticMatrix();
	StaticMatrix(int inputRowNum, int inputColumnNum);

	//int rowNum;
	//int columnNum;

	bool setMatrixElement(int rNum, int cNum, double val);
	double getMatrixElement(int rNum, int cNum);


	//void printMatrix();

	//������
	bool swapRow(int from, int to);

	//������
	bool swapColumn(int from, int to);

	//�����Խ�����Ԫ
	bool swapDiagElement(int from, int to);

	//����������Ϊ��λ����
	bool resetMatrixToI();

	//����������Ϊ��λ����
	void resetMatrixToZero();

	//��ȡָ��������
	BasicVector* getColumnVector(int columnNo);

	//��ȡָ��������
	BasicVector* getRowVector(int rowNo);

	//��ȡָ��������
	BasicVector* getSubMatrixColumnVector(int subMatrixIndex, int columnNo);

	//��ȡָ��������
	BasicVector* getSubMatrixRowVector(int subMatrixIndex, int rowNo);

	//��ȡָ���Խ��Ӿ���hessenberg������
	BasicVector* getSubMatrixHessenColumnVector(int subMatrixIndex);

protected:

	double matrixNxN[10][10];
	StaticVector matrixVector;
	void initMatrix();

	//��ȡ����ָ��Ԫ�ص�ֵ(��Ԫ��ֵ�������Σ�С�ھ��ȵ�ֱֵ�ӷ���0)
	double getMatrixElementRegulared(int rNum, int cNum, double lowEdge);

};

#endif /* MATRIX_STATIC_STATICMATRIX_H_ */
