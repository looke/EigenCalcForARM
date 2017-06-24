/*
 * DoubleShiftQRIteration.cpp
 *
 *  Created on: 2017��5��17��
 *      Author: looke
 */

#include "DoubleShiftQRIteration.h"
#include <iostream>
using namespace std;

//DoubleShiftQRIteration::DoubleShiftQRIteration()
//{};

DoubleShiftQRIteration::DoubleShiftQRIteration(BasicMatrix* p_input_OpMatrix,BasicVector* p_input_TransVector,BasicMatrix* p_input_QTMatrix_Total,BasicMatrix* p_input_QQTMatrix_Step,BasicMatrix* p_input_TempMatrix)
:m_Transposer(),m_Multiplier(p_input_OpMatrix,p_input_OpMatrix,p_input_TempMatrix),m_GivensTrans(p_input_TransVector),m_HouseholderTrans(p_input_TransVector),m_HessenbergForm(p_input_OpMatrix,p_input_QTMatrix_Total,p_input_QQTMatrix_Step,p_input_TempMatrix)
{
	this->init(p_input_OpMatrix,p_input_TransVector,p_input_QTMatrix_Total,p_input_QQTMatrix_Step,p_input_TempMatrix);
};

void DoubleShiftQRIteration::init(BasicMatrix* p_input_OpMatrix,BasicVector* p_input_TransVector,BasicMatrix* p_input_QTMatrix_Total,BasicMatrix* p_input_QQTMatrix_Step,BasicMatrix* p_input_TempMatrix)
{
	//��������
	this->p_OpMatrix = p_input_OpMatrix;

	//���ڼ�����ʽ˫��QR��������ת������Qi/Qn-1 ������---3ά/2ά
	this->p_TransVectorForQStep = p_input_TransVector;

	//Q ���� ��ʽ����  ����Q�����ҳ�OP����
	//this->p_QMatrix_Implicit_Total = p_input_QMatrix_Total;
	//QT ���� ��ʽ���� ����QT�������OP����
	this->p_QTMatrix_Implicit_Total = p_input_QTMatrix_Total;

	//Q/QT ���� ��ʽ���� �ֲ� QT�����ҳ�/���OP����
	this->p_QQTMatrix_Implicit_Step = p_input_QQTMatrix_Step;

	//�м���̾���
	this->p_TempMatrix = p_input_TempMatrix;

	//this->generateHessenbergOpMatrix();
};

void DoubleShiftQRIteration::reload(BasicMatrix* p_input_OpMatrix,BasicVector* p_input_TransVector,BasicMatrix* p_input_QTMatrix_Total,BasicMatrix* p_input_QQTMatrix_Step,BasicMatrix* p_input_TempMatrix)
{
	this->init(p_input_OpMatrix,p_input_TransVector,p_input_QTMatrix_Total,p_input_QQTMatrix_Step,p_input_TempMatrix);
};


/*
 * ����hessenberg��������
 */
void DoubleShiftQRIteration::generateHessenbergOpMatrix()
{
	m_HessenbergForm.reload(this->p_OpMatrix, this->p_QTMatrix_Implicit_Total,this->p_QQTMatrix_Implicit_Step,this->p_TempMatrix);
	m_HessenbergForm.formularUpperHessnbergMatrix();
	//this->p_OpHessenbergMatrix->copyMatrixElementNoCheck(m_HessenbergForm.getOpMatrix());
	//��ʱp_OpMatrix�ѱ����Hessenberg����
	//p_QMatrix_Implicit_Total->copyMatrixElementNoCheck(p_QTMatrix_Implicit_Total);
	//this->m_Transposer.transposeSquareMatrix(p_QMatrix_Implicit_Total);
};

/*
 * ����hessenberg�������2x2�Ӿ��� ������������ֵ�ĺ��Լ��˻�
 */
void DoubleShiftQRIteration::generateWilkinsonShift()
{
	int lastRowIndex = this->p_OpMatrix->rowNum - 1;

	double An_1n_1 = this->p_OpMatrix->getMatrixElement(lastRowIndex-1,lastRowIndex-1);
	double Ann = this->p_OpMatrix->getMatrixElement(lastRowIndex,lastRowIndex);

	double An_1n = this->p_OpMatrix->getMatrixElement(lastRowIndex-1,lastRowIndex);
	double Ann_1 = this->p_OpMatrix->getMatrixElement(lastRowIndex,lastRowIndex-1);

	this->trace = An_1n_1 + Ann;
	this->determinant = An_1n_1*Ann - An_1n*Ann_1;
};

/*
 * B = (A-pI)*(A-tI); ����hessenberg���� ����˫����λ��ľ����һ��
 */
void DoubleShiftQRIteration::generateBMatrixFirstColumn()
{
	generateWilkinsonShift();
	double a11, a12, a21, a22, a32;
	a11 = this->p_OpMatrix->getMatrixElement(0,0);
	a12 = this->p_OpMatrix->getMatrixElement(0,1);
	a21 = this->p_OpMatrix->getMatrixElement(1,0);
	a22 = this->p_OpMatrix->getMatrixElement(1,1);
	a32 = this->p_OpMatrix->getMatrixElement(2,1);

	this->b11 = a11*a11 - a11*this->trace + this->determinant + a12*a21;
	this->b21 = a21*a11 + a21*a22 -a21*this->trace;
	this->b31 = a21*a32;
};

/*
 * Wilkinsonλ��QR���� ��ʽ ��ʼ��
 */
void DoubleShiftQRIteration::initForWilkinsonImplicitQR()
{
	//this->p_QQTSubMatrix_Implicit_Step->resetMatrixToI();
	//this->p_QSubMatrix_Implicit_Step->resetMatrixToI();
	//����ת������
	this->p_QQTMatrix_Implicit_Step->resetMatrixToI();
	//����ת�����󣬴˴�ӦΪ3x3ת���Ӿ���
	this->p_QQTMatrix_Implicit_Step->resizeMatrix(3,3);

	//����wilkinsonƫ�� trace�Լ�determinant
	//generateWilkinsonShift();

	//hessenberg���� ����˫����λ��ľ����һ��
	generateBMatrixFirstColumn();

	//˫����λ��ľ�����������Ԫ�ع�����Ԫ����
	p_TransVectorForQStep->resetDimension(3);
	p_TransVectorForQStep->setElement(0,this->b11);
	p_TransVectorForQStep->setElement(1,this->b21);
	p_TransVectorForQStep->setElement(2,this->b31);

	//����Q1T
	this->m_HouseholderTrans.reload(this->p_TransVectorForQStep);
	//p_QTSubMatrix_Implicit_Step->copyMatrixElementNoCheck(this->p_HouseholderTrans->getHouseholderMatrixToE1(true));
	this->m_HouseholderTrans.getHouseholderMatrixToE1_ReverseElement(p_QQTMatrix_Implicit_Step);

	cout << "initForWilkinsonImplicitQR---Q1T sub" << endl;
	p_QQTMatrix_Implicit_Step->printMatrix();
	//���Ӿ���QT ����Ϊȫ�ߴ�QT����
	upgradeQQTSubMatrix(0);
	cout << "initForWilkinsonImplicitQR---Q1T" << endl;
	p_QQTMatrix_Implicit_Step->printMatrix();
	//��������ת������QT total
	updateQT_Total();
	//����hessenberg���� By QT
	updateHessenbergOpMatrix_By_QT_IM_QRIteration();

	//����Q1
	this->m_Transposer.transposeSquareMatrix(p_QQTMatrix_Implicit_Step);
	//this->p_Transposer->transposeMatrix();
	//this->p_QSubMatrix_Implicit_Step->copyMatrixElementNoCheck(this->p_Transposer->getTransposeMatrix());

	//cout << "initForWilkinsonImplicitQR---Q1 sub" << endl;
	//p_QQTMatrix_Implicit_Step->printMatrix();
	//���Ӿ���Q ����Ϊȫ�ߴ�Q����
	//upgradeQQTSubMatrix(0);
	cout << "initForWilkinsonImplicitQR---Q1" << endl;
	p_QQTMatrix_Implicit_Step->printMatrix();
	//��������ת������Q total
	//updateQ_Total();
	//����hessenberg���� By Q
	updateHessenbergOpMatrix_By_Q_IM_QRIteration();

	cout << "initForWilkinsonImplicitQR---Q1T * Hessenberg * Q1" << endl;
	this->p_OpMatrix->printMatrix();
};

/*
 * ��ת���Ӿ���Q/QT ���ݵ�ǰ��������,����Ϊȫ�ߴ�Q/QT����
 */
void DoubleShiftQRIteration::upgradeQQTSubMatrix(int iterateNum)
{
	this->p_QQTMatrix_Implicit_Step->resizeMatrix(this->p_OpMatrix->rowNum, this->p_OpMatrix->columnNum);

	//���϶Խ���3���Ӿ�������
	p_QQTMatrix_Implicit_Step->moveDiagonalSubMatrixDown(0,2,iterateNum);

	/*
	for(int i=0,m=iterateNum; i<p_QSubMatrix_Implicit_Step->rowNum; i++,m++)
	{
		for(int j=0,n=iterateNum; j<p_QSubMatrix_Implicit_Step->columnNum; j++,n++)
		{
			p_QTMatrix_Implicit_Step->setMatrixElement(m,n,p_QTSubMatrix_Implicit_Step->getMatrixElement(i,j));
			p_QMatrix_Implicit_Step->setMatrixElement(m,n,p_QSubMatrix_Implicit_Step->getMatrixElement(i,j));
		}
	}
	*/
};

/*
 * ��Qn-1�Ӿ���,����Ϊȫ�ߴ�Q����
 */
void DoubleShiftQRIteration::upgradeQQTLastSubMatrix()
{
	int lastQIndex = p_OpMatrix->rowNum-1;

	this->p_QQTMatrix_Implicit_Step->resizeMatrix(this->p_OpMatrix->rowNum, this->p_OpMatrix->columnNum);
	this->p_QQTMatrix_Implicit_Step->moveDiagonalSubMatrixDown(0,1,lastQIndex-1);
	//���϶Խ���2���Ӿ�������������


	//this->p_QTMatrix_Implicit_Step->resetMatrixToI();
	//this->p_QMatrix_Implicit_Step->resetMatrixToI();
	/*
	for(int i=0,m=lastQIndex-1; i<p_QSubMatrix_Implicit_LastStep->rowNum; i++,m++)
	{
		for(int j=0,n=lastQIndex-1; j<p_QSubMatrix_Implicit_LastStep->columnNum; j++,n++)
		{
			p_QTMatrix_Implicit_Step->setMatrixElement(m,n,p_QTSubMatrix_Implicit_LastStep->getMatrixElement(i,j));
			p_QMatrix_Implicit_Step->setMatrixElement(m,n,p_QSubMatrix_Implicit_LastStep->getMatrixElement(i,j));
		}
	}
	*/
};


/*
 * ��ʽQR���� ����hessenberg�������� ����QT * OP-Hessenberg
 */
void DoubleShiftQRIteration::updateHessenbergOpMatrix_By_QT_IM_QRIteration()
{
	//����QT * OP-Hessenberg
	this->m_Multiplier.reload(p_QQTMatrix_Implicit_Step, p_OpMatrix, p_TempMatrix);
	this->m_Multiplier.multiplyCalc();
	p_OpMatrix->copyMatrixElementNoCheck(p_TempMatrix);
};

/*
 * ��ʽQR���� ����hessenberg�������� ���� OP-Hessenberg * Q
 */
void DoubleShiftQRIteration::updateHessenbergOpMatrix_By_Q_IM_QRIteration()
{
	//����QT * OP-Hessenberg * Q
	this->m_Multiplier.reload(p_OpMatrix, p_QQTMatrix_Implicit_Step, p_TempMatrix);
	this->m_Multiplier.multiplyCalc();
	p_OpMatrix->copyMatrixElementNoCheck(p_TempMatrix);
};

/*
 * Wilkinson˫λ��QR���� ��ʽ
 */
void DoubleShiftQRIteration::wilkinson_IM_QRIteration()
{
	p_QTMatrix_Implicit_Total->resetMatrixToI();
	//��������ת��ΪHessenberg����
	this->generateHessenbergOpMatrix();
	cout << "DoubleShiftQRIteration--initForImplicitQR----OP Hessenberg Matrix" << endl;
	this->p_OpMatrix->printMatrix();

	int iteratNumber = 10;
	int i = 0;
	while(i < iteratNumber)
	{
		wilkinson_IM_QRIteration_Step();
		i++;
	}
};

/*
 * Wilkinson˫λ��QR���� ��ʽ ��Ϊ�����ӿڵ���ʹ��
 */
void DoubleShiftQRIteration::wilkinson_IM_QRIteration_Single()
{

	cout << "DoubleShiftQRIteration--initForImplicitQR----OP Hessenberg Matrix" << endl;
	this->p_OpMatrix->printMatrix();

	this->wilkinson_IM_QRIteration_Step();
};
/*
 * Wilkinson˫λ��QR���� ��ʽ -����
 */
void DoubleShiftQRIteration::wilkinson_IM_QRIteration_Step()
{
	//this->p_QTMatrix_Implicit_Total->resetMatrixToI();
	//this->p_QMatrix_Implicit_Total->resetMatrixToI();

	initForWilkinsonImplicitQR();

	//chasing bulge
	int bulgeNumber = this->p_OpMatrix->columnNum-3;
	for(int j=0; j<bulgeNumber; j++)
	{
		//������ʽQ��QT
		double a1,a2,a3;
		a1 = this->p_OpMatrix->getMatrixElement(j+1,j);
		a2 = this->p_OpMatrix->getMatrixElement(j+2,j);
		a3 = this->p_OpMatrix->getMatrixElement(j+3,j);

		//˫����λ ����Ԫ�ع�����Ԫ����
		this->p_TransVectorForQStep->resetDimension(3);
		this->p_TransVectorForQStep->setElement(0,a1);
		this->p_TransVectorForQStep->setElement(1,a2);
		this->p_TransVectorForQStep->setElement(2,a3);

		//����ת������
		this->p_QQTMatrix_Implicit_Step->resetMatrixToI();
		//����ת�����󣬴˴�ӦΪ3x3ת���Ӿ���
		this->p_QQTMatrix_Implicit_Step->resizeMatrix(3,3);

		//����Q1T
		this->m_HouseholderTrans.reload(this->p_TransVectorForQStep);
		//p_QTSubMatrix_Implicit_Step->copyMatrixElementNoCheck(this->p_HouseholderTrans->getHouseholderMatrixToE1(true));
		this->m_HouseholderTrans.getHouseholderMatrixToE1_ReverseElement(p_QQTMatrix_Implicit_Step);

		cout << "wilkinson_IM_QRIteration_Step---QiT Sub:" << j <<endl;
		p_QQTMatrix_Implicit_Step->printMatrix();

		//���Ӿ���QiT ����Ϊȫ�ߴ�QiT����
		upgradeQQTSubMatrix(j+1);

		cout << "wilkinson_IM_QRIteration_Step---QiT Full:" << j <<endl;
		p_QQTMatrix_Implicit_Step->printMatrix();

		//��������ת������QT total
		updateQT_Total();
		//����hessenberg���� By QT
		updateHessenbergOpMatrix_By_QT_IM_QRIteration();

		//����ȫ�ߴ�Qi
		this->m_Transposer.transposeSquareMatrix(p_QQTMatrix_Implicit_Step);
		//this->p_Transposer->transposeMatrix();
		//this->p_QSubMatrix_Implicit_Step->copyMatrixElementNoCheck(this->p_Transposer->getTransposeMatrix());
		cout << "wilkinson_IM_QRIteration_Step---Qi Full:" << j <<endl;
		p_QQTMatrix_Implicit_Step->printMatrix();

		//��������ת������Q total
		//updateQ_Total();
		//����hessenberg����
		updateHessenbergOpMatrix_By_Q_IM_QRIteration();
		cout << "wilkinson_IM_QRIteration_Step---QiT * Hessenberg * Qi" << endl;
		p_OpMatrix->printMatrix();
	}

	endForWilkinsonImplicitQR();
};

/*
 * Wilkinsonλ��QR���� ��ʽ ��β���һ�ε���
 */
void DoubleShiftQRIteration::endForWilkinsonImplicitQR()
{
	int lastQIndex = p_OpMatrix->rowNum-1;
	//������ʽQ��QT
	double an1,an2;
	an1 = this->p_OpMatrix->getMatrixElement(lastQIndex-1,lastQIndex-2);
	an2 = this->p_OpMatrix->getMatrixElement(lastQIndex,lastQIndex-2);

	//˫����λ���һ�� 2��Ԫ�ع���2Ԫ����
	this->p_TransVectorForQStep->resetDimension(2);
	this->p_TransVectorForQStep->setElement(0,an1);
	this->p_TransVectorForQStep->setElement(1,an2);

	//����ת������
	this->p_QQTMatrix_Implicit_Step->resetMatrixToI();
	//����ת�����󣬴˴�ӦΪ2x2ת���Ӿ���
	this->p_QQTMatrix_Implicit_Step->resizeMatrix(2,2);

	this->m_GivensTrans.reload(p_TransVectorForQStep);
	//this->p_QTSubMatrix_Implicit_LastStep->copyMatrixElementNoCheck(p_GivensTrans->getGivensMatrixPreMultiple(1));
	this->m_GivensTrans.getGivensMatrixPreMultiple(1,p_QQTMatrix_Implicit_Step);

	cout << "endForWilkinsonImplicitQR---Qn_1T Sub:" <<endl;
	p_QQTMatrix_Implicit_Step->printMatrix();
	//���Ӿ�������Ϊȫ�ߴ�ת������
	upgradeQQTLastSubMatrix();
	cout << "endForWilkinsonImplicitQR---Qn_1T Full:" <<endl;
	p_QQTMatrix_Implicit_Step->printMatrix();

	//��������ת������QT total
	updateQT_Total();

	//����hessenberg���� By QT
	updateHessenbergOpMatrix_By_QT_IM_QRIteration();

	this->m_Transposer.transposeSquareMatrix(p_QQTMatrix_Implicit_Step);
	//this->p_Transposer->transposeMatrix();
	//this->p_QSubMatrix_Implicit_LastStep->copyMatrixElementNoCheck(this->p_Transposer->getTransposeMatrix());


	cout << "endForWilkinsonImplicitQR---Qn_1 Full:" << endl;
	p_QQTMatrix_Implicit_Step->printMatrix();

	//���Ӿ�������Ϊȫ�ߴ�ת������
	//upgradeQQTLastSubMatrix();
	//cout << "endForWilkinsonImplicitQR---Qn_1T" <<endl;
	//p_QTMatrix_Implicit_Step->printMatrix();
	//cout << "endForWilkinsonImplicitQR---Qn_1" << endl;
	//p_QMatrix_Implicit_Step->printMatrix();

	//��������ת������Q total
	//updateQ_Total();

	//����hessenberg���� By Q
	updateHessenbergOpMatrix_By_Q_IM_QRIteration();

	cout << "endForWilkinsonImplicitQR---Qn_1T * Hessenberg * Qn_1" << endl;
	p_OpMatrix->printMatrix();
};

/*
 * ��������ת������Q QT
 */
void DoubleShiftQRIteration::updateQT_Total()
{
	//����QT_Step * QT_Total
	this->m_Multiplier.reload(p_QQTMatrix_Implicit_Step, p_QTMatrix_Implicit_Total, p_TempMatrix);
	this->m_Multiplier.multiplyCalc();
	this->p_QTMatrix_Implicit_Total->copyMatrixElementNoCheck(p_TempMatrix);

};

/*
 * ��������ת������Q QT
 */
//void DoubleShiftQRIteration::updateQ_Total()
//{
	//����Q_Total * Q_Step
//	this->m_Multiplier.reload(p_QMatrix_Implicit_Total, p_QQTMatrix_Implicit_Step, p_TempMatrix);
//	this->m_Multiplier.multiplyCalc();
//	this->p_QMatrix_Implicit_Total->copyMatrixElementNoCheck(p_TempMatrix);
//};

BasicMatrix* DoubleShiftQRIteration::getOpHessenbergMatrix()
{
	return this->p_OpMatrix;
};

//��ȡ����ת������
BasicMatrix* DoubleShiftQRIteration::getQTMatrix_Total()
{
	return this->p_QTMatrix_Implicit_Total;
};

//BasicMatrix* DoubleShiftQRIteration::getQMatrix_Total()
//{
//	return this->p_QMatrix_Implicit_Total;
//};
