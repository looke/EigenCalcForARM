/*
 * DoubleShiftQZIteration.h
 *
 *  Created on: 2017��5��29��
 *      Author: looke
 */

#ifndef EIGEN_BASIC_DOUBLESHIFTQZITERATION_H_
#define EIGEN_BASIC_DOUBLESHIFTQZITERATION_H_

#include "BasicMatrix.h"
#include "BasicVector.h"
#include "ABInverseCalculator.h"
#include "QRDecomposition.h"
#include "HessenbergTriangleFormular.h"
#include "HouseholderTransformation.h"

class DoubleShiftQZIteration
{
public:
	//DoubleShiftQZIteration();
	DoubleShiftQZIteration(BasicMatrix* p_input_OpMatrix_A, BasicMatrix* p_input_OpMatrix_B,BasicVector* p_input_TransVectorForQZStep,BasicMatrix* p_input_QMatrix_Total,BasicMatrix* p_input_ZMatrix_Total,BasicMatrix* p_input_QZMatrix_Step,BasicMatrix* p_input_TempMatrix_Trans,BasicMatrix* p_input_TempMatrix);

	void init(BasicMatrix* p_input_OpMatrix_A, BasicMatrix* p_input_OpMatrix_B,BasicVector* p_input_TransVectorForQZStep,BasicMatrix* p_input_QMatrix_Total,BasicMatrix* p_input_ZMatrix_Total,BasicMatrix* p_input_QZMatrix_Step,BasicMatrix* p_input_TempMatrix_Trans,BasicMatrix* p_input_TempMatrix);
	void reload(BasicMatrix* p_input_OpMatrix_A, BasicMatrix* p_input_OpMatrix_B,BasicVector* p_input_TransVectorForQZStep,BasicMatrix* p_input_QMatrix_Total,BasicMatrix* p_input_ZMatrix_Total,BasicMatrix* p_input_QZMatrix_Step,BasicMatrix* p_input_TempMatrix_Trans,BasicMatrix* p_input_TempMatrix);
	//virtual ~DoubleShiftQZIteration(){};

	//����H-T���������
	void generateHessenTriangleOpMatrix();

	//Wilkinsonλ��QZ���� ��ʽ ��ʼ��
	void initForWilkinsonImplicitQZ();

	//Wilkinson˫λ��QZ���� ��ʽ -��β
	void endForWilkinsonImplicitQZ();

	//Wilkinsonλ��QZ���� ��ʽ ���ֵ���
	void wilkinson_IM_QZIteration_Step();
	//Wilkinsonλ��QZ���� ��ʽ ���ֵ�������ӿ�
	void wilkinson_IM_QZIteration_Single();
	//Wilkinsonλ��QZ���� ��ʽ 10����������
	void wilkinson_IM_QZIteration();

	//����A*B^-1 �������2x2�Ӿ��� ��������Wilkinson��λ��trace & determinant
	void generateWilkinsonShift();

	//���� A*B^-1 ǰ����Ԫ��
	void generateABinvFirst3x2();

	//���� A*B^-1 ˫����λ�� ��һ��Ԫ��
	void generateABinvShiftedFirstColumn();

	//���ݵ�ǰ�������ȣ��������Ǿ������Ͳ�λ��Ԫ��
	void generateTriangleBulgeElement(int iterateNum);

	//��QZ�Ӿ�����ݵ�ǰ��������,����Ϊȫ�ߴ�QZ����
	void upgradeQSubMatrix(int iterateNum);
	void upgradeZSubMatrix(int iterateNum);
	//��2x2 Q�Ӿ�������Ϊ 3x3 Q�Ӿ��� 2x2λ������
	void upgradeQMiniToSubMatrix_RightEnd();
	//��2x2 Z�Ӿ�������Ϊ 3x3 Z�Ӿ��� 2x2λ������
	void upgradeZMiniToSubMatrix_RightEnd();
	//��2x2 Z�Ӿ�������Ϊ 3x3 Z�Ӿ��� 2x2λ������
	void upgradeZMiniToSubMatrix_LeftTop();

	//ʹ��Q�������H-T�����
	void updateHTMatrixByQ();
	//ʹ��Z�����ҳ�H-T�����
	void updateHTMatrixByZ();
	//�����ۺ�ת������Q/Z Total(��Step�ϲ���Total)
	void updateQMatrix_Total();
	void updateZMatrix_Total();

	//��ȡQ/Z �ۺ�ת������
	BasicMatrix* getQMatrix_Total();
	BasicMatrix* getZMatrix_Total();
protected:
	//p+t ����wilkinsonλ��ֵ�ĺ�
	double trace;

	//p*t ����wilkinsonλ��ֵ�ĳ˻�
	double determinant;

	//A*B^-1 ǰ����Ԫ��
	double ABinv_11;
	double ABinv_12;
	double ABinv_21;
	double ABinv_22;
	double ABinv_32;

	// C = (AB^-1 -pI)*(AB^-1 -tI); ˫����λ��ľ����һ��
	double C11;
	double C21;
	double C31;

    //���Ǿ������Ͳ��ֵ�Ԫ��
	double B_21;
	double B_22;
	double B_31;
	double B_32;
	double B_33;
	//��������
	BasicMatrix* p_OpMatrix_A;
	BasicMatrix* p_OpMatrix_B;

	//Hessenberg����
	//BasicMatrix* p_OpMatrix_Hessenberg;
	//Triangle����
	//BasicMatrix* p_OpMatrix_Triangle;

	//���ڼ�����ʽ˫��QZ��������ת������Qi������---3ά
	BasicVector* p_TransVectorForQZStep;

	//Q ȫ�ߴ���� ��ʽ����  ����Q�������OP����
	BasicMatrix* p_QMatrix_Implicit_Total;
	//Z ȫ�ߴ���� ��ʽ���� ����Z�����ҳ�OP����
	BasicMatrix* p_ZMatrix_Implicit_Total;

	//���ڼ�����ʽ˫��QZ��������ת������Qi������---3ά
	//BasicVector* p_TransVectorForQStep_3;
	//���ڼ�����ʽ˫��QZ��������ת������Qn-1������---2ά
	//BasicVector* p_TransVectorForQStep_2;

	//���ڼ�����ʽ˫��QZ��������ת������Zi������---3ά
	//BasicVector* p_TransVectorForZStep_3;
	//���ڼ�����ʽ˫��QZ��������ת������Zi������---2ά
	//BasicVector* p_TransVectorForZStep_2;

	//QZ ȫ�ߴ���� ��ʽ���� �ֲ� Q�������OP����
	BasicMatrix* p_QZMatrix_Implicit_Step;

	//Q ȫ�ߴ���� ��ʽ���� �ֲ� Q�������OP����
	//BasicMatrix* p_QMatrix_Implicit_Step;
	//Z ȫ�ߴ���� ��ʽ���� �ֲ� Z�����ҳ�OP����
	//BasicMatrix* p_ZMatrix_Implicit_Step;


	//�м���̾���
	BasicMatrix* p_TempMatrix_Trans;
	//�м���̾���
	BasicMatrix* p_TempMatrix;

	//Qi �Ӿ��� ��ʽ���� �ֲ� Qi�������OP���� 3x3
	//BasicMatrix* p_QSubMatrix_Implicit_Step_3;
	//Qi �Ӿ��� ��ʽ���� �ֲ� Qn-1�������OP���� 2x2
	//BasicMatrix* p_QSubMatrix_Implicit_Step_2;

	//Zi �Ӿ��� ��ʽ���� �ֲ� QTi�����ҳ�OP���� 3x3
	//BasicMatrix* p_ZSubMatrix_Implicit_Step_3;
	//Zi �Ӿ��� ��ʽ���� �ֲ� QT�����ҳ�OP���� 2x2
	//BasicMatrix* p_ZSubMatrix_Implicit_Step_2;

	//hessenberg-triangle��ʽ��
	HessenbergTriangleFormular m_HessenbergTriangleFormular;


	//householder�任
	HouseholderTransformation m_HouseholderTrans;

	//�˷���
	MatrixMultiplier m_Multiplier;

	//AB^-1 Ԫ�ؼ�����
	ABInverseCalculator m_ABInvCalc;
};



#endif /* EIGEN_BASIC_DOUBLESHIFTQZITERATION_H_ */
