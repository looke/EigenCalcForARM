/*
 * BasicMatrix.h
 *
 *  Created on: 2017��4��1��
 *      Author: looke
 */

#ifndef MATRIX_BASICMATRIX_H_
#define MATRIX_BASICMATRIX_H_

#include "BasicVector.h"
#include "math.h"

class BasicMatrix
{
public:
	BasicMatrix();
	BasicMatrix(int inputRowNum, int inputColumnNum);

	int rowNum;
	int columnNum;
	int space;

	//���þ���ָ��Ԫ�ص�ֵ
	virtual bool setMatrixElement(int rNum, int cNum, double val);

	//��ȡ����ָ��Ԫ�ص�ֵ
	virtual double getMatrixElement(int rNum, int cNum);

	//��ӡ����
	virtual void printMatrix();

	//������
	virtual bool swapRow(int from, int to);

	//������
	virtual bool swapColumn(int from, int to);

	//�����Խ�����Ԫ
	virtual bool swapDiagElement(int from, int to);

	//����������Ϊ��λ����
	virtual bool resetMatrixToI();

	//����������Ϊ��λ����
	virtual void resetMatrixToZero();

	//��ȡָ��������
	virtual BasicVector* getColumnVector(int columnNo);

	//��ȡָ��������
	virtual BasicVector* getRowVector(int rowNo);

	//��ȡָ���Խ��Ӿ���������
	virtual BasicVector* getSubMatrixColumnVector(int subMatrixIndex, int columnNo);

	//��ȡָ���Խ��Ӿ���������
	virtual BasicVector* getSubMatrixRowVector(int subMatrixIndex, int rowNo);

	//��ȡָ���Խ��Ӿ���hessenberg������
	virtual BasicVector* getSubMatrixHessenColumnVector(int subMatrixIndex);

	//���������������Ԫ�ص�����ֵ(���������)
	double calcMaxDifferentialWithCheck(BasicMatrix* targetMatrix);

	//���������������Ԫ�ص�����ֵ(�����������)
	double calcMaxDifferentialNoCheck(BasicMatrix* targetMatrix);

	//��������Ԫ��(���������)
	bool copyMatrixElementWithCheck(BasicMatrix* input_Matrix);

	//��������Ԫ��(�����������)
	void copyMatrixElementNoCheck(BasicMatrix* input_Matrix);

	//�ж��Ƿ�Ϊ�ԳƷ���
	bool isInputSymmetryMatrix();

	//��������Ԫ�س������
	void rMatrix(double r);

	//����ָ���е�����Ԫ�ط��ŷ�ת
	void reverseSignOfColumn(int columnIndex);

	//�Ծ�������Ԫ�ؽ�������,������ֵС�ھ��ȵ�Ԫ������Ϊ0
	void regularZeroElement();

	//Ѱ�����Խ�����Ϊ0��Ԫ��,��������λ��
	//����ж��0Ԫ,��������ֵ�����Ǹ�
	int maxIndexOfZeroOnDiagonal();

	//�Խ���ȫ��Ԫ�ؼ�ȥָ��ֵ
	void diagonalSubtraction(double subValue);
	//�Խ���ȫ��Ԫ�ؼ���ָ��ֵ
	void diagonalAddition(double addValue);

	//��ʼ��������Ϣ
	void initPrecision();

	//get������Ϣ
	double getPrecision();

	//��������Frobenious���� ||A||f
	double FrobeniousNorm();

	/*
	 * ����������ȥ����
	 */
	double getLowEdge();

	//�����趨����ά��
	bool resizeMatrix(int row, int column);

	//�жϾ����Ƿ�����upper hessenberg���� ���ڼ�Сֵ��0����
	bool isUpperHessenbergMatrix();
	//�жϾ����Ƿ�����upper Triangle���� ���ڼ�Сֵ��0����
	bool isUpperTriangleMatrix();

	//�ضԽ��� �����ƶ�ָ���Խ��Ӿ��� �ƶ�ָ������
	bool moveDiagonalSubMatrixDown(int headIndex, int tailIndex, int steps);

	//�ضԽ��� �����ƶ�ָ���Խ��Ӿ��� �ƶ�ָ������
	bool moveDiagonalSubMatrixUp(int headIndex, int tailIndex, int steps);
	//��ȡ����ָ��Ԫ�ص�ֵ(��Ԫ��ֵ�������Σ�С�ھ��ȵ�ֱֵ�ӷ���0)
	virtual double getMatrixElementRegulared(int rNum, int cNum, double lowEdge);
	virtual ~BasicMatrix() {};
protected:

	double precision;
	virtual void initMatrix();


};


#endif /* MATRIX_BASICMATRIX_H_ */
