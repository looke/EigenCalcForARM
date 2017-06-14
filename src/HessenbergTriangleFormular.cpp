/*
 * HessenbergTriangleFormular.cpp
 *
 *  Created on: 2017��4��28��
 *      Author: looke
 */

#include "HessenbergTriangleFormular.h"
#include <iostream>
using namespace std;

//HessenbergTriangleFormular::HessenbergTriangleFormular()
//{};

HessenbergTriangleFormular::HessenbergTriangleFormular(BasicMatrix* input_Matrix_A, BasicMatrix* input_Matrix_B,BasicMatrix* input_QMatrix_Total,BasicMatrix* input_ZMatrix_Total,BasicMatrix* input_QZMatrix_Step,BasicMatrix* input_TempMatrix_Trans,BasicMatrix* input_TempMatrix):m_Transposer(),m_Multiplier(input_Matrix_A, input_Matrix_B,input_TempMatrix),m_GivensTrans(input_Matrix_A->getColumnVector(0)),m_QRDecomp(input_Matrix_B, input_QZMatrix_Step, input_TempMatrix_Trans, input_TempMatrix)
{
	this->init(input_Matrix_A, input_Matrix_B, input_QMatrix_Total, input_ZMatrix_Total, input_QZMatrix_Step, input_TempMatrix_Trans, input_TempMatrix);
};

void HessenbergTriangleFormular::init(BasicMatrix* input_Matrix_A, BasicMatrix* input_Matrix_B,BasicMatrix* input_QMatrix_Total,BasicMatrix* input_ZMatrix_Total,BasicMatrix* input_QZMatrix_Step,BasicMatrix* input_TempMatrix_Trans,BasicMatrix* input_TempMatrix)
{
	//ԭʼ��������A
	this->p_OpMatrix_A = input_Matrix_A;
	//ԭʼ��������B
	this->p_OpMatrix_B = input_Matrix_B;

	//����Q�任����,���ڻ�����˲�������
	this->p_QMatrix_Total = input_QMatrix_Total;

	//����Z�任����,���ڻ����ҳ˲�������
	this->p_ZMatrix_Total = input_ZMatrix_Total;

	//����QZ�任����,������˲�������
	this->p_QZMatrix_Step = input_QZMatrix_Step;

	//�м�ת������
	this->p_TempMatrix_Trans = input_TempMatrix_Trans;

	//�м���̾���
	this->p_TempMatrix = input_TempMatrix;
};

void HessenbergTriangleFormular::reload(BasicMatrix* input_Matrix_A, BasicMatrix* input_Matrix_B,BasicMatrix* input_QMatrix_Total,BasicMatrix* input_ZMatrix_Total,BasicMatrix* input_QZMatrix_Step,BasicMatrix* input_TempMatrix_Trans,BasicMatrix* input_TempMatrix)
{
	this->init(input_Matrix_A, input_Matrix_B, input_QMatrix_Total, input_ZMatrix_Total, input_QZMatrix_Step, input_TempMatrix_Trans, input_TempMatrix);
};

BasicMatrix* HessenbergTriangleFormular::getMatrixA()
{
	return this->p_OpMatrix_A;
};

BasicMatrix* HessenbergTriangleFormular::getMatrixB()
{
	return this->p_OpMatrix_B;
};
BasicMatrix* HessenbergTriangleFormular::getHessenbergMatrixA()
{
	return this->p_OpMatrix_A;
};
BasicMatrix* HessenbergTriangleFormular::getTriangleMatrixB()
{
	return this->p_OpMatrix_B;
};
BasicMatrix* HessenbergTriangleFormular::getMatrixQ_Total()
{
	return this->p_QMatrix_Total;
};

BasicMatrix* HessenbergTriangleFormular::getMatrixZ_Total()
{
	return this->p_ZMatrix_Total;
};

//BasicMatrix* HessenbergTriangleFormular::getMatrixQ_Step()
//{
//	return this->p_QMatrix_Step;
//};

//BasicMatrix* HessenbergTriangleFormular::getMatrixZ_Step()
//{
//	return this->p_ZMatrix_Step;
//};

BasicMatrix* HessenbergTriangleFormular::getMatrixQZ_Step()
{
	return this->p_QZMatrix_Step;
};

/*
 * �γ�A-Hessenberg B-Triangle �����
 */
void HessenbergTriangleFormular::formularABMatrix()
{
	p_QMatrix_Total->resetMatrixToI();
	p_ZMatrix_Total->resetMatrixToI();
	//p_QMatrix_Step->resetMatrixToI();
	p_QZMatrix_Step->resetMatrixToI();
	if(this->p_OpMatrix_A->isUpperHessenbergMatrix() && this->p_OpMatrix_B->isUpperTriangleMatrix())
	{
		//this->p_OpHessenbergMatrix_A->copyMatrixElementNoCheck(p_OpMatrix_A);
		//this->p_OpTriangleMatrix_B->copyMatrixElementNoCheck(p_OpMatrix_B);
		return;
	}

	this->initABMatrix();
	cout << "After HessenBergTriangle init--- A:" << endl;
	p_OpMatrix_A->printMatrix();
	cout << "After HessenBergTriangle init--- B:" << endl;
	p_OpMatrix_B->printMatrix();
	cout << "After HessenBergTriangle init--- Q Total:" << endl;
	p_QMatrix_Total->printMatrix();
	cout << "After HessenBergTriangle init--- Z Total:" << endl;
	p_ZMatrix_Total->printMatrix();

	for(int i=0; i<this->p_OpMatrix_A->columnNum; i++)
	{
		this->formularColumnVector(i);
		//p_OpHessenbergMatrix_A->regularZeroElement();
		//p_OpTriangleMatrix_B->regularZeroElement();
	}
};


/*
 * ��A�����ָ���н���hessenberg��ʽ��
 */
void HessenbergTriangleFormular::formularColumnVector(int columnIndex)
{
	double element,lowEdge;
	BasicVector* p_ColumnVector;
	//BasicMatrix* p_GivensMatrix;
	for(int i=this->p_OpMatrix_A->rowNum-1; i>columnIndex+1; i--)
	{
		p_QZMatrix_Step->resetMatrixToI();

		lowEdge = p_OpMatrix_A->getLowEdge();
		//�ȼ������A����
		element = p_OpMatrix_A->getMatrixElementRegulared(i,columnIndex,lowEdge);
		if(0 != element)
		{
			p_ColumnVector = p_OpMatrix_A->getColumnVector(columnIndex);

			this->m_GivensTrans.reload(p_ColumnVector);
			this->m_GivensTrans.getGivensMatrixPreMultiple(i, p_QZMatrix_Step);
			cout << "After HessenBergTriangle--- Left Givens Matrix:" << i << endl;
			p_QZMatrix_Step->printMatrix();;
		}


		updateOpMatrix_A_ByQ();
		updateOpMatrix_B_ByQ();
		updateQMatrix_Total();

		p_QZMatrix_Step->resetMatrixToI();

		//���ż������B����
		lowEdge = p_OpMatrix_B->getLowEdge();
		element = p_OpMatrix_B->getMatrixElementRegulared(i,i-1,lowEdge);
		if(0 != element)
		{
			p_ColumnVector = p_OpMatrix_B->getRowVector(i);
			cout << "p_OpTriangleMatrix_B:getRowVector:" << i << endl;
			p_ColumnVector->printVector();
			this->m_GivensTrans.reload(p_ColumnVector);
			//this->p_GivensTrans->setIsUsingPreElement(false);
			this->m_GivensTrans.getGivensMatrixAfterMultiple(i-1, p_QZMatrix_Step);
			cout << "After HessenBergTriangle--- Right Givens Matrix:" << i << endl;
			p_QZMatrix_Step->printMatrix();
		}

		updateOpMatrixByZ();
		updateZMatrix_Total();

		cout << "After HessenBergTriangle--- A:" << i << endl;
		p_OpMatrix_A->printMatrix();
		cout << "After HessenBergTriangle--- B:" << i<< endl;
		p_OpMatrix_B->printMatrix();
		cout << "After HessenBergTriangle--- Q:" << i << endl;
		p_QMatrix_Total->printMatrix();
		cout << "After HessenBergTriangle--- Z:" << i << endl;
		p_ZMatrix_Total->printMatrix();
	}
};


/*
 * ʹ��QT���� ���²�������
 *  QT * A
 *  QT * B
 */
void HessenbergTriangleFormular::updateOpMatrix_A_ByQ()
{
	//���²�������A
	m_Multiplier.reload(p_QZMatrix_Step, p_OpMatrix_A, p_TempMatrix);
	m_Multiplier.multiplyCalc();
	this->p_OpMatrix_A->copyMatrixElementNoCheck(p_TempMatrix);
	//this->p_OpMatrix_A->regularZeroElement();
};

void HessenbergTriangleFormular::updateOpMatrix_B_ByQ()
{
	//���²�������B
	m_Multiplier.reload(p_QZMatrix_Step, p_OpMatrix_B, p_TempMatrix);
	m_Multiplier.multiplyCalc();
	this->p_OpMatrix_B->copyMatrixElementNoCheck(p_TempMatrix);
	//this->p_OpTriangleMatrix_B->regularZeroElement();
};
/*
 * ʹ��Z���� ���²�������
 *  A * Z
 *  B * Z
 */
void HessenbergTriangleFormular::updateOpMatrixByZ()
{
	//���²�������A
	m_Multiplier.reload(p_OpMatrix_A, p_QZMatrix_Step, p_TempMatrix);
	m_Multiplier.multiplyCalc();
	this->p_OpMatrix_A->copyMatrixElementNoCheck(p_TempMatrix);
	//this->p_OpHessenbergMatrix_A->regularZeroElement();

	//���²�������B
	m_Multiplier.reload(p_OpMatrix_B, p_QZMatrix_Step, p_TempMatrix);
	m_Multiplier.multiplyCalc();
	this->p_OpMatrix_B->copyMatrixElementNoCheck(p_TempMatrix);
	//this->p_OpTriangleMatrix_B->regularZeroElement();

};

void HessenbergTriangleFormular::updateQMatrix_Total()
{
	//����Q����
	m_Multiplier.reload(p_QZMatrix_Step, p_QMatrix_Total, p_TempMatrix);
	m_Multiplier.multiplyCalc();
	this->p_QMatrix_Total->copyMatrixElementNoCheck(p_TempMatrix);
};

void HessenbergTriangleFormular::updateZMatrix_Total()
{
	//����Z����
	m_Multiplier.reload(p_ZMatrix_Total, p_QZMatrix_Step, p_TempMatrix);
	m_Multiplier.multiplyCalc();
	this->p_ZMatrix_Total->copyMatrixElementNoCheck(p_TempMatrix);
};

/*
 * ԭʼ��������B ����QR�ֽ� ��Ϊ�����Ǿ��󣬲���QT���A
 */
void HessenbergTriangleFormular::initABMatrix()
{
	//this->p_OpHessenbergMatrix_A->copyMatrixElementNoCheck(this->p_OpMatrix_A);
	//this->p_OpTriangleMatrix_B->copyMatrixElementNoCheck(this->p_OpMatrix_B);

	//��B�������QR�ֽ⣬
	this->m_QRDecomp.reload(p_OpMatrix_B, p_QZMatrix_Step, p_TempMatrix_Trans, p_TempMatrix);
	this->m_QRDecomp.calcQRMatrix();
	//Q�任����ΪQR�ֽ����Q�����ת��
	m_Transposer.transposeSquareMatrix(p_QZMatrix_Step);

	//this->p_QMatrix_Step->copyMatrixElementNoCheck(this->p_QRDecomp->getQTMatrix());
	updateOpMatrix_A_ByQ();
	updateQMatrix_Total();

	//Z�任����Ϊ��λ��
	//this->p_QZMatrix_Step->resetMatrixToI();
	//updateOpMatrixByZ();
	//updateZMatrix_Total();
};
