/*
 * SingleShiftQZIteration.h
 *
 *  Created on: 2017��5��26��
 *      Author: looke
 */

#ifndef EIGEN_BASIC_SINGLESHIFTQZITERATION_H_
#define EIGEN_BASIC_SINGLESHIFTQZITERATION_H_

#include "BasicMatrix.h"
#include "MatrixMultiplier.h"
#include "HessenbergTriangleFormular.h"
#include "GivensTransformation.h"
#include "ABInverseCalculator.h"

class SingleShiftQZIteration
{
public:
	//SingleShiftQZIteration();
	SingleShiftQZIteration(BasicMatrix* p_input_OpMatrix_A, BasicMatrix* p_input_OpMatrix_B,BasicMatrix* p_input_QMatrix_Total,BasicMatrix* p_input_ZMatrix_Total,BasicMatrix* p_input_QZMatrix_Step,BasicMatrix* p_input_TempMatrix_Trans,BasicMatrix* p_input_TempMatrix);

	void init(BasicMatrix* p_input_OpMatrix_A, BasicMatrix* p_input_OpMatrix_B,BasicMatrix* p_input_QMatrix_Total,BasicMatrix* p_input_ZMatrix_Total,BasicMatrix* p_input_QZMatrix_Step,BasicMatrix* p_input_TempMatrix_Trans,BasicMatrix* p_input_TempMatrix);
	void reload(BasicMatrix* p_input_OpMatrix_A, BasicMatrix* p_input_OpMatrix_B,BasicMatrix* p_input_QMatrix_Total,BasicMatrix* p_input_ZMatrix_Total,BasicMatrix* p_input_QZMatrix_Step,BasicMatrix* p_input_TempMatrix_Trans,BasicMatrix* p_input_TempMatrix);
	//virtual ~SingleShiftQZIteration(){};

	//����H-T���������
	void generateHessenTriangleOpMatrix();

	//����B^-1���2�� ���ڼ��㵥��λ��ֵ
	//void generateBinvLastTwoRow();

	//����AB^-1 ���½�Ԫ��ֵ ���ڵ���λ��
	//double generateABinvLastOne();

	//��ʼ����ʽ����
	void initForImplicitQZ(double input_ShiftValue);

	//��ֵλ��QZ���� ��ʽ����
	void implicit_QZIteration_Step(double input_shiftValue);

	//��ֵλ��QZ���� ��ʽ
	void implicit_QZIteration(double input_shiftValue);

	//��ֵrayleigh��λ��QZ���� ��ʽ ���
	void rayleigh_Quotient_IM_QZIteration(int iterateNum);

	//��ֵrayleigh��λ��QZ���� ��ʽ ����
	void rayleigh_Quotient_IM_QZIteration_Step();

	//ʹ��Q�������H-T�����
	void updateHTMatrixByQ();
	//ʹ��Z�����ҳ�H-T�����
	void updateHTMatrixByZ();

	//�����ۺ�ת������Q/Z Total(��Step�ϲ���Total)
	void updateQMatrix_Total();
	void updateZMatrix_Total();
	//��ȡQ/Z ����ת������
	BasicMatrix* getQMatrix_Total();
	BasicMatrix* getZMatrix_Total();
protected:
	//B^-1 �������Ԫ��
	//double Binv_n_1_n_1;
	//double Binv_n_1_n;
	//double Binv_n_n;

	//ԭʼ��������A Hessenberg
	BasicMatrix* p_OpMatrix_A;
	//ԭʼ��������B Triangle
	BasicMatrix* p_OpMatrix_B;

	//Z ���� ��ʽ���� �ۺ� Z�����ҳ�OP����
	BasicMatrix* p_ZMatrix_Implicit_Total;
	//Q ���� ��ʽ���� �ۺ� Q�������OP����
	BasicMatrix* p_QMatrix_Implicit_Total;

	//Q/Z ���� ��ʽ���� �ֲ� Q�������OP���� Z�����ҳ�OP����
	BasicMatrix* p_QZMatrix_Implicit_Step;

	//�м���̾���
	BasicMatrix* p_TempMatrix_Trans;
	//�м���̾���
	BasicMatrix* p_TempMatrix;

	//Hessenberg ��������A
	//BasicMatrix* p_OpMatrix_Hessenberg;
	//Triangle ��������B
	//BasicMatrix* p_OpMatrix_Triangle;

	//Z ���� ��ʽ���� �ֲ� Z�����ҳ�OP����
	//BasicMatrix* p_ZMatrix_Implicit_Step;
	//Q ���� ��ʽ���� �ֲ� Q�������OP����
	//BasicMatrix* p_QMatrix_Implicit_Step;



	//hessenberg-triangle��ʽ��
	HessenbergTriangleFormular m_HessenbergTriangleFormular;

	//QZ-triangle 0Ԫ��λ
	//QZTriangleZeroChasing* p_QZTriangleZeroChasing;

	//Givens�任��
	GivensTransformation m_GivensTrans;

	//�˷���
	MatrixMultiplier m_Multiplier;

	//AB^-1 Ԫ�ؼ�����
	ABInverseCalculator m_ABInvCalc;
};



#endif /* EIGEN_BASIC_SINGLESHIFTQZITERATION_H_ */
