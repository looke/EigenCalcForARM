/*
 * QRDecomposition.cpp
 *
 *  Created on: 2017��4��25��
 *      Author: looke
 */


#include "QRDecomposition.h"
#include <iostream>
using namespace std;

//QRDecomposition::QRDecomposition()
//{};

QRDecomposition::QRDecomposition(BasicMatrix* p_input_OpMatrix, BasicMatrix* p_input_QMatrix, BasicMatrix* p_input_householderMatrix, BasicMatrix* p_input_TempMatrix):m_Multiplier(p_input_OpMatrix,p_input_householderMatrix,p_input_TempMatrix),m_HouseholderTrans(p_input_OpMatrix->getColumnVector(0))
{
	this->init(p_input_OpMatrix, p_input_QMatrix, p_input_householderMatrix, p_input_TempMatrix);
};

void QRDecomposition::init(BasicMatrix* p_input_OpMatrix, BasicMatrix* p_input_QMatrix, BasicMatrix* p_input_householderMatrix, BasicMatrix* p_input_TempMatrix)
{
	this->p_OpMatrix = p_input_OpMatrix;
	this->p_QMatrix = p_input_QMatrix;
	this->p_householderMatrix = p_input_householderMatrix;
	this->p_TempMatrix = p_input_TempMatrix;
};

void QRDecomposition::reload(BasicMatrix* p_input_OpMatrix, BasicMatrix* p_input_QMatrix, BasicMatrix* p_input_householderMatrix, BasicMatrix* p_input_TempMatrix)
{
	this->init(p_input_OpMatrix, p_input_QMatrix, p_input_householderMatrix, p_input_TempMatrix);
};

BasicMatrix* QRDecomposition::getQMatrix()
{
	return this->p_QMatrix;
};

BasicMatrix* QRDecomposition::getQTMatrix()
{
	return this->p_QMatrix;
};

BasicMatrix* QRDecomposition::getRMatrix()
{
	return this->p_OpMatrix;
};

BasicMatrix* QRDecomposition::getHouseholderMatrix()
{
	return this->p_householderMatrix;
};
/*
 * ����OP����ĸ�������,����Householder�任����,�𲽽�OP����Խǻ�
 */
void QRDecomposition::calcQRMatrix()
{
	this->p_householderMatrix->resetMatrixToI();
	this->p_QMatrix->resetMatrixToI();
	//this->p_QTMatrix->resetMatrixToI();
	int iterationMax = this->p_OpMatrix->rowNum-1;

	for(int i=0; i<iterationMax; i++)
	{
		generateHouseholderMatrix(i);
		updateMatrix();
	}

};

//����HouseHolder�任����  index ��ʾ��ǰ�������еĶԽ���Ԫ�ص�����
void QRDecomposition::generateHouseholderMatrix(int index)
{
	BasicVector *p_currentColumnVector = this->p_OpMatrix->getSubMatrixColumnVector(index,0);
	int vectorSize = p_currentColumnVector->getDimension();
	m_HouseholderTrans.reload(p_currentColumnVector);

	//���µ���householderMatrixά�ȣ�ƥ���Ӧ����ת������
	p_householderMatrix->resizeMatrix(vectorSize, vectorSize);
	//p_householderMatrix->printMatrix();
	m_HouseholderTrans.getHouseholderMatrixToE1_ReverseElement(p_householderMatrix);
	//p_householderMatrix->printMatrix();
	//��ԭhouseholderMatrixά�ȣ�ƥ��ԭʼ����ά��
	p_householderMatrix->resizeMatrix(p_OpMatrix->rowNum, p_OpMatrix->columnNum);
	//p_householderMatrix->printMatrix();

	p_householderMatrix->moveDiagonalSubMatrixDown(0,p_OpMatrix->columnNum-1-index,index);
	//p_householderMatrix->printMatrix();
	//cout << "Sub Householder Matrix:" << endl;
	//p_householderMatrix->printMatrix();
	//����Householder����
	//this->p_householderMatrix->resetMatrixToI();
	//int m=0;
	//int n=0;
	//for(int i=index, m=0;i<this->p_householderMatrix->rowNum; i++,m++)
	//{
		//n=0;
	//	for(int j=index, n=0; j<this->p_householderMatrix->columnNum; j++,n++)
	//	{
	//		p_householderMatrix->setMatrixElement(i,j,p_subHouseholderMatrix->getMatrixElement(m,n));
			//n++;
	//	}
		//m++;
	//}

};

/*
 * ��������ε�����householder�����,����OP�����Q����
 */
void QRDecomposition::updateMatrix()
{
	//update op matrix
	this->m_Multiplier.reload(this->p_householderMatrix, this->p_OpMatrix, p_TempMatrix);
	this->m_Multiplier.multiplyCalc();
	this->p_OpMatrix->copyMatrixElementNoCheck(p_TempMatrix);
	//this->p_OpMatrix->regularZeroElement();
	//update QT Matrix
	//this->m_Multiplier.reload(this->p_householderMatrix, this->p_QTMatrix);
	//this->m_Multiplier.multiplyCalc();
	//this->p_QTMatrix->copyMatrixElementNoCheck(this->p_Multiplier->getMultiplyResult());
	//this->p_QTMatrix->regularZeroElement();
	//Update Q matrix
	this->m_Multiplier.reload(this->p_QMatrix, this->p_householderMatrix, p_TempMatrix);
	this->m_Multiplier.multiplyCalc();
	this->p_QMatrix->copyMatrixElementNoCheck(p_TempMatrix);
	//this->p_QMatrix->regularZeroElement();
	p_householderMatrix->resetMatrixToI();

	//cout << "OP Matrix After QR update:" << endl;
	//this->p_OpMatrix->printMatrix();
	//cout << "QT Matrix After QR update:" << endl;
	//this->p_QTMatrix->printMatrix();
	//cout << "Q Matrix After QR update:" << endl;
	//this->p_QMatrix->printMatrix();
};
