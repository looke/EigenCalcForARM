/*
 * SingleShiftQZIteration.cpp
 *
 *  Created on: 2017��5��27��
 *      Author: looke
 */

#include "SingleShiftQZIteration.h"
//#include <iostream>
using namespace std;

//SingleShiftQZIteration::SingleShiftQZIteration()
//{};

SingleShiftQZIteration::SingleShiftQZIteration(BasicMatrix* p_input_OpMatrix_A, BasicMatrix* p_input_OpMatrix_B,BasicMatrix* p_input_QMatrix_Total,BasicMatrix* p_input_ZMatrix_Total,BasicMatrix* p_input_QZMatrix_Step,BasicMatrix* p_input_TempMatrix_Trans,BasicMatrix* p_input_TempMatrix)
:m_ABInvCalc(),m_GivensTrans(p_input_OpMatrix_A->getColumnVector(0)),m_Multiplier(p_input_OpMatrix_A,p_input_OpMatrix_B,p_input_TempMatrix),m_HessenbergTriangleFormular(p_input_OpMatrix_A,p_input_OpMatrix_B,p_input_QMatrix_Total,p_input_ZMatrix_Total,p_input_QZMatrix_Step,p_input_TempMatrix_Trans,p_input_TempMatrix)
{
	this->init(p_input_OpMatrix_A,p_input_OpMatrix_B,p_input_QMatrix_Total,p_input_ZMatrix_Total,p_input_QZMatrix_Step,p_input_TempMatrix_Trans,p_input_TempMatrix);
};

void SingleShiftQZIteration::init(BasicMatrix* p_input_OpMatrix_A, BasicMatrix* p_input_OpMatrix_B,BasicMatrix* p_input_QMatrix_Total,BasicMatrix* p_input_ZMatrix_Total,BasicMatrix* p_input_QZMatrix_Step,BasicMatrix* p_input_TempMatrix_Trans,BasicMatrix* p_input_TempMatrix)
{
	//ԭʼ��������A Hessenberg
	this->p_OpMatrix_A = p_input_OpMatrix_A;
	//ԭʼ��������B Triangle
	this->p_OpMatrix_B = p_input_OpMatrix_B;

	//Q ���� ��ʽ���� �ۺ� Z�����ҳ�OP����
	this->p_QMatrix_Implicit_Total = p_input_QMatrix_Total;
	//Z ���� ��ʽ���� �ۺ� Q�������OP����
	this->p_ZMatrix_Implicit_Total = p_input_ZMatrix_Total;

	//Q/Z ���� ��ʽ���� �ֲ� Q�������OP���� Z�����ҳ�OP����
	this->p_QZMatrix_Implicit_Step = p_input_QZMatrix_Step;

	//�м���̾���
	this->p_TempMatrix_Trans = p_input_TempMatrix_Trans;
	//�м���̾���
	this->p_TempMatrix = p_input_TempMatrix;

	//this->generateHessenTriangleOpMatrix();
};

void SingleShiftQZIteration::reload(BasicMatrix* p_intput_OpMatrix_A, BasicMatrix* p_intput_OpMatrix_B,BasicMatrix* p_input_QMatrix_Total,BasicMatrix* p_input_ZMatrix_Total,BasicMatrix* p_input_QZMatrix_Step,BasicMatrix* p_input_TempMatrix_Trans,BasicMatrix* p_input_TempMatrix)
{
	this->init(p_intput_OpMatrix_A,p_intput_OpMatrix_B,p_input_QMatrix_Total,p_input_ZMatrix_Total,p_input_QZMatrix_Step,p_input_TempMatrix_Trans,p_input_TempMatrix);
};


/*
 * ����H-T���������
 */
void SingleShiftQZIteration::generateHessenTriangleOpMatrix()
{
	this->m_HessenbergTriangleFormular.reload(p_OpMatrix_A, p_OpMatrix_B,p_QMatrix_Implicit_Total,p_ZMatrix_Implicit_Total,p_QZMatrix_Implicit_Step,p_TempMatrix_Trans,p_TempMatrix);
	this->m_HessenbergTriangleFormular.formularABMatrix();
	//this->p_OpMatrix_Hessenberg->copyMatrixElementNoCheck(p_HessenbergTriangleFormular->getHessenbergMatrixA());
	//this->p_OpMatrix_Triangle->copyMatrixElementNoCheck(p_HessenbergTriangleFormular->getTriangleMatrixB());

	//cout << "SingleShiftQZIteration--generateHessenTriangleOpMatrix----OP Hessenberg Matrix" << endl;
	//this->p_OpMatrix_A->printMatrix();
	//cout << "SingleShiftQZIteration--generateHessenTriangleOpMatrix----OP Triangle Matrix" << endl;
	//this->p_OpMatrix_B->printMatrix();
};

/*
 * ʹ��Q�������H-T�����
 */
void SingleShiftQZIteration::updateHTMatrixByQ()
{
	//cout << "SingleShiftQZIteration--updateHTMatrixByQ----OP Hessenberg Matrix Before" << endl;
	//this->p_OpMatrix_A->printMatrix();
	//cout << "SingleShiftQZIteration--updateHTMatrixByQ----OP Triangle Matrix Before" << endl;
	//this->p_OpMatrix_B->printMatrix();

	//���¾���Hessenberg A
	this->m_Multiplier.reload(p_QZMatrix_Implicit_Step, p_OpMatrix_A, p_TempMatrix);
	this->m_Multiplier.multiplyCalc();
	p_OpMatrix_A->copyMatrixElementNoCheck(p_TempMatrix);

	//���¾���Triangle B
	this->m_Multiplier.reload(p_QZMatrix_Implicit_Step, p_OpMatrix_B,p_TempMatrix);
	this->m_Multiplier.multiplyCalc();
	p_OpMatrix_B->copyMatrixElementNoCheck(p_TempMatrix);

	//cout << "SingleShiftQZIteration--updateHTMatrixByQ----OP Hessenberg Matrix After" << endl;
	//this->p_OpMatrix_A->printMatrix();
	//cout << "SingleShiftQZIteration--updateHTMatrixByQ----OP Triangle Matrix After" << endl;
	//this->p_OpMatrix_B->printMatrix();
};

/*
 * ʹ��Z�����ҳ�H-T�����
 */
void SingleShiftQZIteration::updateHTMatrixByZ()
{
	//cout << "SingleShiftQZIteration--updateHTMatrixByZ----OP Hessenberg Matrix Before" << endl;
	//this->p_OpMatrix_A->printMatrix();
	//cout << "SingleShiftQZIteration--updateHTMatrixByZ----OP Triangle Matrix Before" << endl;
	//this->p_OpMatrix_B->printMatrix();

	//���¾���Hessenberg A
	this->m_Multiplier.reload(p_OpMatrix_A, p_QZMatrix_Implicit_Step, p_TempMatrix);
	this->m_Multiplier.multiplyCalc();
	p_OpMatrix_A->copyMatrixElementNoCheck(p_TempMatrix);

	//���¾���Triangle B
	this->m_Multiplier.reload(p_OpMatrix_B, p_QZMatrix_Implicit_Step, p_TempMatrix);
	this->m_Multiplier.multiplyCalc();
	p_OpMatrix_B->copyMatrixElementNoCheck(p_TempMatrix);

	//cout << "SingleShiftQZIteration--updateHTMatrixByZ----OP Hessenberg Matrix After" << endl;
	//this->p_OpMatrix_A->printMatrix();
	//cout << "SingleShiftQZIteration--updateHTMatrixByZ----OP Triangle Matrix After" << endl;
	//this->p_OpMatrix_B->printMatrix();
};

/*
 * �����ۺ�ת������Q/Z Total(��Step�ϲ���Total)
 */
void SingleShiftQZIteration::updateQMatrix_Total()
{
	//���¾���Q total
	this->m_Multiplier.reload(p_QZMatrix_Implicit_Step, p_QMatrix_Implicit_Total, p_TempMatrix);
	this->m_Multiplier.multiplyCalc();
	p_QMatrix_Implicit_Total->copyMatrixElementNoCheck(p_TempMatrix);

};

/*
 * �����ۺ�ת������Q/Z Total(��Step�ϲ���Total)
 */
void SingleShiftQZIteration::updateZMatrix_Total()
{
	//���¾���Z total
	this->m_Multiplier.reload(p_ZMatrix_Implicit_Total, p_QZMatrix_Implicit_Step, p_TempMatrix);
	this->m_Multiplier.multiplyCalc();
	p_ZMatrix_Implicit_Total->copyMatrixElementNoCheck(p_TempMatrix);
};

/*
 * ��ʼ����ʽ����
 * ִ����ʽ������һ������H-T����ԵĴζԽ������γ�һ��ͻ����λ
 */
void SingleShiftQZIteration::initForImplicitQZ(double input_ShiftValue)
{
	//cout << "SingleShiftQZIteration--initForImplicitQZ" << endl;
	//p_QMatrix_Implicit_Step->resetMatrixToI();
	p_QZMatrix_Implicit_Step->resetMatrixToI();

	//QZ������λ A*B^-1 ��һ��
	double temp_old_a00 = this->p_OpMatrix_A->getMatrixElement(0,0);
	double temp_old_a10 = this->p_OpMatrix_A->getMatrixElement(1,0);

	double temp_old_b00 = this->p_OpMatrix_B->getMatrixElement(0,0);


	double temp_new_a00 = temp_old_a00/temp_old_b00 - input_ShiftValue;
	double temp_new_a10 = temp_old_a10/temp_old_b00;

	//���µ�һ��
	this->p_OpMatrix_A->setMatrixElement(0,0,temp_new_a00);
	this->p_OpMatrix_A->setMatrixElement(1,0,temp_new_a10);

	//������ʽQ
	BasicVector* p_firstColumnVector = this->p_OpMatrix_A->getColumnVector(0);
	m_GivensTrans.reload(p_firstColumnVector);
	m_GivensTrans.getGivensMatrixPreMultiple(1,p_QZMatrix_Implicit_Step);

	//cout << "SingleShiftQZIteration--initForImplicitQZ-----Q-Step" << endl;
	//p_QZMatrix_Implicit_Step->printMatrix();
	//��ԭHessenberg����
	this->p_OpMatrix_A->setMatrixElement(0,0,temp_old_a00);
	this->p_OpMatrix_A->setMatrixElement(1,0,temp_old_a10);

	//����H-T�����
	this->updateHTMatrixByQ();
	//�����ۺ�ת������Q Total
	updateQMatrix_Total();

	//������ʽZ
	BasicVector* p_secondRowVector = this->p_OpMatrix_B->getRowVector(1);
	m_GivensTrans.reload(p_secondRowVector);
	m_GivensTrans.getGivensMatrixAfterMultiple(0,p_QZMatrix_Implicit_Step);

	//cout << "SingleShiftQZIteration--initForImplicitQZ-----Z-Step" << endl;
	//p_QZMatrix_Implicit_Step->printMatrix();
	//����H-T�����
	this->updateHTMatrixByZ();
	//�����ۺ�ת������Z Total
	updateZMatrix_Total();

	//cout << "SingleShiftQZIteration--initForImplicitQZ-----END" << endl;

};

//��ֵλ��QZ���� ��ʽ����
void SingleShiftQZIteration::implicit_QZIteration_Step(double input_shiftValue)
{
	//p_QMatrix_Implicit_Total->resetMatrixToI();
	//p_ZMatrix_Implicit_Total->resetMatrixToI();

	this->initForImplicitQZ(input_shiftValue);

	int bulgeNumber = this->p_OpMatrix_A->columnNum-2;

	for(int j=0; j<bulgeNumber; j++)
	{
		//������ʽQ
		BasicVector* p_columnVector = this->p_OpMatrix_A->getColumnVector(j);
		m_GivensTrans.reload(p_columnVector);
		//this->p_QMatrix_Implicit_Step->copyMatrixElementNoCheck(p_GivensTrans->getGivensMatrixPreMultiple(j+2));
		m_GivensTrans.getGivensMatrixPreMultiple(j+2,p_QZMatrix_Implicit_Step);
		//����H-T�����
		this->updateHTMatrixByQ();
		//�����ۺ�ת������Q Total
		updateQMatrix_Total();

		//������ʽZ
		BasicVector* p_secondRowVector = this->p_OpMatrix_B->getRowVector(j+2);
		m_GivensTrans.reload(p_secondRowVector);
		//this->p_ZMatrix_Implicit_Step->copyMatrixElementNoCheck(p_GivensTrans->getGivensMatrixAfterMultiple(j+1));
		m_GivensTrans.getGivensMatrixAfterMultiple(j+1,p_QZMatrix_Implicit_Step);
		//����H-T�����
		this->updateHTMatrixByZ();
		//�����ۺ�ת������Q\Z Total
		updateZMatrix_Total();
	}

};

//��ֵλ��QZ���� ��ʽ
void SingleShiftQZIteration::implicit_QZIteration(double input_shiftValue)
{
	//generateHessenTriangleOpMatrix();
	int i=0;
	while (i<10)
	{
		implicit_QZIteration_Step(input_shiftValue);
		//cout << "SingleShiftQZIteration--implicit_QZIteration----OP Hessenberg Matrix after:" << i << "iteration" <<endl;
		//p_OpMatrix_A->printMatrix();
		i++;
	}
};
//��ֵrayleigh��λ��QZ���� ��ʽ
void SingleShiftQZIteration::rayleigh_Quotient_IM_QZIteration(int iterateNum)
{
	p_QMatrix_Implicit_Total->resetMatrixToI();
	p_ZMatrix_Implicit_Total->resetMatrixToI();

	this->generateHessenTriangleOpMatrix();

	//int rayleighValueIndex = p_OpMatrix_A->rowNum - 1;
	//double rayleighValue;
	int i=0;
	while (i<iterateNum)
	{
		this->rayleigh_Quotient_IM_QZIteration_Step();
		//rayleighValue = p_OpMatrix_A->getMatrixElement(rayleighValueIndex,rayleighValueIndex);
		//implicit_QZIteration_Step(rayleighValue);
		//cout << "SingleShiftQZIteration--rayleigh_Quotient_IM_QRIteration----OP Hessenberg Matrix after:" << i << "iteration" <<endl;
		//p_OpMatrix_Hessenberg->printMatrix();
		i++;
	}
};

//��ֵrayleigh��λ��QZ���� ��ʽ ����
void SingleShiftQZIteration::rayleigh_Quotient_IM_QZIteration_Step()
{
	//����A*B^-1 ���½�Ԫ�� ��Ϊ����λ��ֵ
	m_ABInvCalc.generateABinvLastOne(p_OpMatrix_A, p_OpMatrix_B);
	double rayleighValue = m_ABInvCalc.getABinv_N_N();

	implicit_QZIteration_Step(rayleighValue);
};

//��ȡQ ����ת������
BasicMatrix* SingleShiftQZIteration::getQMatrix_Total()
{
	return this->p_QMatrix_Implicit_Total;
};

//��ȡZ ����ת������
BasicMatrix* SingleShiftQZIteration::getZMatrix_Total()
{
	return this->p_ZMatrix_Implicit_Total;
};
