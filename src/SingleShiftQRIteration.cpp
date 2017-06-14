/*
 * SingelShiftQRIteration.cpp
 *
 *  Created on: 2017��5��13��
 *      Author: looke
 */

#include "SingleShiftQRIteration.h"
//#include <iostream>
//using namespace std;

//SingleShiftQRIteration::SingleShiftQRIteration()
//{};

SingleShiftQRIteration::SingleShiftQRIteration(BasicMatrix* p_input_OpMatrix, BasicMatrix* p_input_QTMatrix_Implicit_Total,BasicMatrix* p_input_Q_QT_Matrix_Implicit_Step,BasicMatrix* p_input_TempMatrix):m_Transposer(),m_Multiplier(p_input_OpMatrix,p_input_OpMatrix,p_input_TempMatrix),m_GivensTrans(p_input_OpMatrix->getColumnVector(0)),m_HessenbergForm(p_input_OpMatrix, p_input_QTMatrix_Implicit_Total, p_input_Q_QT_Matrix_Implicit_Step, p_input_TempMatrix),m_QRDecomp(p_input_OpMatrix, p_input_QTMatrix_Implicit_Total, p_input_Q_QT_Matrix_Implicit_Step, p_input_TempMatrix)
{
	this->init(p_input_OpMatrix, p_input_QTMatrix_Implicit_Total, p_input_Q_QT_Matrix_Implicit_Step, p_input_TempMatrix);
};

void SingleShiftQRIteration::init(BasicMatrix* p_input_OpMatrix,BasicMatrix* p_input_QTMatrix_Implicit_Total,BasicMatrix* p_input_Q_QT_Matrix_Implicit_Step,BasicMatrix* p_input_TempMatrix)
{
	//��������
	this->p_OpMatrix = p_input_OpMatrix;

	//Q ���� ��ʽ���� �ۺ� Q�����ҳ�OP����
	//this->p_QMatrix_Implicit_Total = p_input_QMatrix_Implicit_Total;
	//QT ���� ��ʽ���� �ۺ� QT�������OP����
	this->p_QTMatrix_Implicit_Total = p_input_QTMatrix_Implicit_Total;

	//Q/QT ���� ��ʽ���� �ֲ� Q�����ҳ�OP���� QT�������OP����
	this->p_Q_QT_Matrix_Implicit_Step = p_input_Q_QT_Matrix_Implicit_Step;

	//�м���̾���
	this->p_TempMatrix = p_input_TempMatrix;

	generateHessenbergOpMatrix();
};

void SingleShiftQRIteration::reload(BasicMatrix* p_input_OpMatrix, BasicMatrix* p_input_QTMatrix_Implicit_Total,BasicMatrix* p_input_Q_QT_Matrix_Implicit_Step,BasicMatrix* p_input_TempMatrix)
{
	this->init(p_input_OpMatrix, p_input_QTMatrix_Implicit_Total, p_input_Q_QT_Matrix_Implicit_Step, p_input_TempMatrix);
};

/*
 * ����hessenberg��������
 */
void SingleShiftQRIteration::generateHessenbergOpMatrix()
{
	m_HessenbergForm.reload(p_OpMatrix, p_QTMatrix_Implicit_Total, p_Q_QT_Matrix_Implicit_Step, p_TempMatrix);
	m_HessenbergForm.formularUpperHessnbergMatrix();
	//this->p_OpHessenbergMatrix->copyMatrixElementNoCheck(p_HessenbergForm->getOpMatrix());
	//p_QMatrix_Implicit_Total->copyMatrixElementNoCheck(p_QTMatrix_Implicit_Total);
	//m_Transposer.transposeSquareMatrix(p_QMatrix_Implicit_Total);
};


/*
 * ��ֵλ��QR���� ��ʽ
 */
/*
void SingleShiftQRIteration::explicit_QRIteration(double input_shiftValue)
{
	int iterationNumber = 10;

	//this->generateHessenbergOpMatrix();
	cout << "SingelShiftQRIteration--explicit_QRIteration----OP Hessenberg Matrix" << endl;
	this->p_OpHessenbergMatrix->printMatrix();

	int i=0;
	while (i<iterationNumber)
	{
		explicit_QRIteration_Step(input_shiftValue);
		cout << "SingelShiftQRIteration--explicit_QRIteration----OP Hessenberg Matrix after:" << i << "iteration" <<endl;
		p_OpHessenbergMatrix->printMatrix();
		i++;
	}
};
*/
/*
//��ֵλ��QR���� ��ʽ����
void SingleShiftQRIteration::explicit_QRIteration_Step(double input_shiftValue)
{
	//�Խ�Ԫ��ȥ�ƶ�ֵ
	p_OpHessenbergMatrix->diagonalSubtraction(input_shiftValue);
	cout << "SingelShiftQRIteration--explicit_QRIteration_Step----OP Hessenberg Matrix after Subtraction" <<endl;
	p_OpHessenbergMatrix->printMatrix();

	this->p_QRDecomp->reload(p_OpHessenbergMatrix);
	this->p_QRDecomp->calcQRMatrix();
	this->p_QMatrix_Explicit->copyMatrixElementNoCheck(this->p_QRDecomp->getQMatrix());
	p_OpHessenbergMatrix->copyMatrixElementNoCheck(this->p_QRDecomp->getRMatrix());

	cout << "SingelShiftQRIteration--explicit_QRIteration_Step----OP Hessenberg Matrix after QRDecomp" <<endl;
	p_OpHessenbergMatrix->printMatrix();

	this->p_Multiplier->reload(p_OpHessenbergMatrix, p_QMatrix_Explicit);
	this->p_Multiplier->multiplyCalc();
	p_OpHessenbergMatrix->copyMatrixElementNoCheck(p_Multiplier->getMultiplyResult());

	cout << "SingelShiftQRIteration--explicit_QRIteration_Step----OP Hessenberg Matrix after R*Q" <<endl;
	p_OpHessenbergMatrix->printMatrix();
	//�Խ�Ԫ�����ƶ�ֵ
	p_OpHessenbergMatrix->diagonalAddition(input_shiftValue);
	cout << "SingelShiftQRIteration--explicit_QRIteration_Step----OP Hessenberg Matrix after Addition" <<endl;
	p_OpHessenbergMatrix->printMatrix();
};
*/

//��ֵλ��QR���� ��ʽ ��ʼ��
void SingleShiftQRIteration::initForImplicitQR(double input_shiftValue)
{
	this->p_Q_QT_Matrix_Implicit_Step->resetMatrixToI();
	//������λ
	double temp_old = this->p_OpMatrix->getMatrixElement(0,0);
	double temp_new = temp_old - input_shiftValue;
	this->p_OpMatrix->setMatrixElement(0,0,temp_new);


	BasicVector* p_firstColumnVector = this->p_OpMatrix->getColumnVector(0);
	//��ԭOpMatrix
	this->p_OpMatrix->setMatrixElement(0,0,temp_old);

	//������ʽQ��QT
	//������ʽQT ��˾���
	m_GivensTrans.reload(p_firstColumnVector);
	m_GivensTrans.getGivensMatrixPreMultiple(1,p_Q_QT_Matrix_Implicit_Step);

	//���²������� A-Hess
	this->updateOpMatrix_By_QT_IM_QRIteration();
	//��������ת������QT
	this->updateQTMatrix_Total_IM_QRIteration();

	//cout << "SingelShiftQRIteration--initForImplicitQR----QTMatrix_Implicit_Step" <<endl;
	//p_Q_QT_Matrix_Implicit_Step->printMatrix();

	//this->p_QTMatrix_Implicit_Step->copyMatrixElementNoCheck(p_GivensTrans->getGivensMatrixPreMultiple(1));

	//������ʽQ �ҳ˾���
	this->m_Transposer.transposeSquareMatrix(p_Q_QT_Matrix_Implicit_Step);

	//���²������� A-Hess
	this->updateOpMatrix_By_Q_IM_QRIteration();
	//��������ת������Q
	//this->updateQMatrix_Total_IM_QRIteration();

	//cout << "SingelShiftQRIteration--initForImplicitQR----QMatrix_Implicit_Step" <<endl;
	//p_Q_QT_Matrix_Implicit_Step->printMatrix();

	//cout << "SingelShiftQRIteration--initForImplicitQR----QTMatrix_Implicit_Total" <<endl;
	//p_QTMatrix_Implicit_Total->printMatrix();
	//cout << "SingelShiftQRIteration--initForImplicitQR----QMatrix_Implicit_Total" <<endl;
	//p_QMatrix_Implicit_Total->printMatrix();
	//�˴��Ѿ��γ���ʽ����QR�����ĳ�ʼ����������
	//cout << "SingelShiftQRIteration--initForImplicitQR----OP Hessenberg Matrix after init"<<endl;
	//p_OpMatrix->printMatrix();
};

//��ʽQR���� ����hessenberg�������� QT*A
void SingleShiftQRIteration::updateOpMatrix_By_QT_IM_QRIteration()
{
	//����QT * OP-Hessenberg
	this->m_Multiplier.reload(p_Q_QT_Matrix_Implicit_Step, p_OpMatrix,p_TempMatrix);
	this->m_Multiplier.multiplyCalc();
	p_OpMatrix->copyMatrixElementNoCheck(p_TempMatrix);

};

//��ʽQR���� ����hessenberg�������� A*Q
void SingleShiftQRIteration::updateOpMatrix_By_Q_IM_QRIteration()
{
	//���� OP-Hessenberg * Q
	this->m_Multiplier.reload(p_OpMatrix, p_Q_QT_Matrix_Implicit_Step, p_TempMatrix);
	this->m_Multiplier.multiplyCalc();
	p_OpMatrix->copyMatrixElementNoCheck(p_TempMatrix);
};

//��ֵλ��QR���� ��ʽ
void SingleShiftQRIteration::implicit_QRIteration(double input_shiftValue)
{
	//��������ת��ΪHessenberg����
	//this->generateHessenbergOpMatrix();
	//cout << "SingelShiftQRIteration--initForImplicitQR----OP Hessenberg Matrix" << endl;
	//this->p_OpMatrix->printMatrix();

	int iteratNumber = 10;
	int i = 0;
	while(i < iteratNumber)
	{
		implicit_QRIteration_Step(input_shiftValue);
		//cout << "SingelShiftQRIteration--implicit_QRIteration----OP Hessenberg Matrix after: " << i << " iteration" <<endl;
		//p_OpMatrix->printMatrix();
		i++;
	}
};

//��ֵλ��QR���� ��ʽ
void SingleShiftQRIteration::implicit_QRIteration_Step(double input_shiftValue)
{
	//p_QTMatrix_Implicit_Total->resetMatrixToI();
	//p_QMatrix_Implicit_Total->resetMatrixToI();
	initForImplicitQR(input_shiftValue);

	int bulgeNumber = this->p_OpMatrix->columnNum-2;
	for(int j=0; j<bulgeNumber; j++)
	{
		//������ʽQ��QT
		BasicVector* p_columnVector = this->p_OpMatrix->getColumnVector(j);
		m_GivensTrans.reload(p_columnVector);
		//����QT
		m_GivensTrans.getGivensMatrixPreMultiple(j+2,p_Q_QT_Matrix_Implicit_Step);

		//����A-Hessenberg����
		this->updateOpMatrix_By_QT_IM_QRIteration();
		//��������ת������ QT
		this->updateQTMatrix_Total_IM_QRIteration();

		//����Q
		m_Transposer.transposeSquareMatrix(p_Q_QT_Matrix_Implicit_Step);
		//this->p_Transposer->transposeMatrix();
		//this->->copyMatrixElementNoCheck(this->p_Transposer->getTransposeMatrix());
		//����A-Hessenberg����
		this->updateOpMatrix_By_Q_IM_QRIteration();

		//��������ת������Q
		//this->updateQMatrix_Total_IM_QRIteration();
	}
};

/*
 * ������λ��QR���� ��ʽ
 */
/*
void SingleShiftQRIteration::rayleigh_Quotient_EX_QRIteration()
{
	int iterationNumber = 10;

	//this->generateHessenbergOpMatrix();
	cout << "SingelShiftQRIteration--rayleigh_Quotient_EX_QRIteration----OP Hessenberg Matrix" << endl;
	this->p_OpHessenbergMatrix->printMatrix();

	int rayleighValueIndex = p_OpHessenbergMatrix->rowNum - 1;
	double rayleighValue;
	int i=0;
	while (i<iterationNumber)
	{
		rayleighValue = p_OpHessenbergMatrix->getMatrixElement(rayleighValueIndex,rayleighValueIndex);
		explicit_QRIteration_Step(rayleighValue);
		cout << "SingelShiftQRIteration--rayleigh_Quotient_EX_QRIteration----OP Hessenberg Matrix after:" << i << "iteration" <<endl;
		p_OpHessenbergMatrix->printMatrix();
		i++;
	}
};
*/

//������λ��QR���� ��ʽ 10�ε���
void SingleShiftQRIteration::rayleigh_Quotient_IM_QRIteration()
{
	this->rayleigh_Quotient_IM_QRIteration(10);
}

//������λ��QR���� ��ʽ ָ����������
void SingleShiftQRIteration::rayleigh_Quotient_IM_QRIteration(int iterateNum)
{
	int iterationNumber = iterateNum;

	//this->generateHessenbergOpMatrix();
	//cout << "SingelShiftQRIteration--rayleigh_Quotient_IM_QRIteration----OP Hessenberg Matrix" << endl;
	//this->p_OpMatrix->printMatrix();

	int rayleighValueIndex = p_OpMatrix->rowNum - 1;
	double rayleighValue;
	int i=0;
	while (i<iterationNumber)
	{
		rayleighValue = p_OpMatrix->getMatrixElement(rayleighValueIndex,rayleighValueIndex);
		implicit_QRIteration_Step(rayleighValue);
		//cout << "SingelShiftQRIteration--rayleigh_Quotient_IM_QRIteration----OP Hessenberg Matrix after:" << i << "iteration" <<endl;
		//p_OpMatrix->printMatrix();
		i++;
	}
};

//������λ��QR���� ��ʽ ����
void SingleShiftQRIteration::rayleigh_Quotient_IM_QRIteration_Step()
{
	//this->generateHessenbergOpMatrix();
	//cout << "SingelShiftQRIteration--rayleigh_Quotient_IM_QRIteration Step----OP Hessenberg Matrix" << endl;
	//this->p_OpMatrix->printMatrix();

	int rayleighValueIndex = p_OpMatrix->rowNum - 1;
	double rayleighValue= p_OpMatrix->getMatrixElement(rayleighValueIndex,rayleighValueIndex);

	implicit_QRIteration_Step(rayleighValue);
	//cout << "SingelShiftQRIteration--rayleigh_Quotient_IM_QRIteration Step----OP Hessenberg Matrix" <<endl;
	//p_OpMatrix->printMatrix();
};

/*
 * �����ۺ�ת������QT Total
 */
void SingleShiftQRIteration::updateQTMatrix_Total_IM_QRIteration()
{
	//��������ۺϾ��� QT_Step * QT_Total
	this->m_Multiplier.reload(p_Q_QT_Matrix_Implicit_Step, p_QTMatrix_Implicit_Total,p_TempMatrix);
	this->m_Multiplier.multiplyCalc();
	p_QTMatrix_Implicit_Total->copyMatrixElementNoCheck(p_TempMatrix);
};

/*
 * �����ۺ�ת������Q Total
 */
//void SingleShiftQRIteration::updateQMatrix_Total_IM_QRIteration()
//{
	//�����ҳ��ۺϾ��� Q_Total * Q_Step
//	this->m_Multiplier.reload(p_QMatrix_Implicit_Total, p_Q_QT_Matrix_Implicit_Step, p_TempMatrix);
//	this->m_Multiplier.multiplyCalc();
//	p_QMatrix_Implicit_Total->copyMatrixElementNoCheck(p_TempMatrix);
//};

//BasicMatrix* SingleShiftQRIteration::getOpHessenbergMatrix()
//{
//	return this->p_OpHessenbergMatrix;
//};

/*
 * ��ȡ����ת������ QT
 */
BasicMatrix* SingleShiftQRIteration::getQTMatrix_Total()
{
	return this->p_QTMatrix_Implicit_Total;
};

/*
 * ��ȡ����ת������ Q
 */
//BasicMatrix* SingleShiftQRIteration::getQMatrix_Total()
//{
//	return this->p_QMatrix_Implicit_Total;
//};
