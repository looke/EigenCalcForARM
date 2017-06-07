/*
 * StaticMatrix.h
 *
 *  Created on: 2017��4��25��
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
	//��ȡ����ָ��Ԫ�ص�ֵ(��Ԫ��ֵ�������Σ�С�ھ��ȵ�ֱֵ�ӷ���0)
	double getMatrixElementRegulared(int rNum, int cNum, double lowEdge);

	void printMatrix();

	//������
	void swapRow(int from, int to);

	//������
	void swapColumn(int from, int to);

	//�����Խ�����Ԫ
	void swapDiagElement(int from, int to);

	//����������Ϊ��λ����
	void resetMatrixToI();

	//����������Ϊ��λ����
	void resetMatrixToZero();

	//��ȡָ��������
	void getColumnVector(int columnNo, BasicVector* p_Vector);

	//��ȡָ��������
	void getRowVector(int rowNo, BasicVector* p_Vector);

	//��ȡָ��������
	void getSubMatrixColumnVector(int subMatrixIndex, int columnNo, BasicVector* p_Vector);

	//��ȡָ��������
	void getSubMatrixRowVector(int subMatrixIndex, int rowNo, BasicVector* p_Vector);

	//��ȡָ���Խ��Ӿ���hessenberg������
	void getSubMatrixHessenColumnVector(int subMatrixIndex, BasicVector* p_Vector);

private:

	double matrixNxN[20][20];
	//StaticVector columnVector;
	//StaticVector rowVector;
	void initMatrix();

};

#endif /* MATRIX_STATIC_STATICMATRIX_H_ */
