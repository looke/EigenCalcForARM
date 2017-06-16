/*
 * DoubleShiftQRIteration.h
 *
 *  Created on: 2017��5��17��
 *      Author: looke
 */

#ifndef EIGEN_BASIC_DOUBLESHIFTQRITERATION_H_
#define EIGEN_BASIC_DOUBLESHIFTQRITERATION_H_

#include "BasicMatrix.h"
#include "MatrixTransposer.h"
#include "BasicVector.h"
#include "HessenbergFormular.h"
#include "GivensTransformation.h"
#include "HouseholderTransformation.h"

class DoubleShiftQRIteration
{
public:
	//DoubleShiftQRIteration();
	DoubleShiftQRIteration(BasicMatrix* p_input_OpMatrix,BasicVector* p_input_TransVector,BasicMatrix* p_input_QTMatrix_Total,BasicMatrix* p_input_QQTMatrix_Step,BasicMatrix* p_input_TempMatrix);

	void init(BasicMatrix* p_input_OpMatrix,BasicVector* p_input_TransVector,BasicMatrix* p_input_QTMatrix_Total,BasicMatrix* p_input_QQTMatrix_Step,BasicMatrix* p_input_TempMatrix);
	void reload(BasicMatrix* p_input_OpMatrix,BasicVector* p_input_TransVector,BasicMatrix* p_input_QTMatrix_Total,BasicMatrix* p_input_QQTMatrix_Step,BasicMatrix* p_input_TempMatrix);
	//virtual ~DoubleShiftQRIteration(){};

	//Wilkinson˫λ��QR����  ��ʽ
	//void wilkinson_EX_QRIteration();

	//Wilkinson˫λ��QR���� ��ʽ
	void wilkinson_IM_QRIteration();
	//��Ϊ�����ӿڵ���ʹ��
	void wilkinson_IM_QRIteration_Single();

	BasicMatrix* getOpHessenbergMatrix();

	//��ȡ����ת������
	BasicMatrix* getQTMatrix_Total();
	//BasicMatrix* getQMatrix_Total();

	//Wilkinson˫λ��QR���� ��ʽ -����
	void wilkinson_IM_QRIteration_Step();

	//Wilkinsonλ��QR���� ��ʽ ��ʼ��
	void initForWilkinsonImplicitQR();

	//Wilkinson˫λ��QR���� ��ʽ -��β
	void endForWilkinsonImplicitQR();

	//����hessenberg��������
	void generateHessenbergOpMatrix();

	//����hessenberg�������2x2�Ӿ��� ��������
	void generateWilkinsonShift();

	//B = (A-pI)*(A-tI); ����hessenberg���� ����˫����λ��ľ����һ��
	void generateBMatrixFirstColumn();

	//��Q�Ӿ�����ݵ�ǰ��������,����Ϊȫ�ߴ�Q����
	void upgradeQQTSubMatrix(int iterateNum);

	//��Qn-1�Ӿ���,����Ϊȫ�ߴ�Q����
	void upgradeQQTLastSubMatrix();

	//��ʽQR���� ����hessenberg��������
	//void updateHessenbergOpMatrix_IM_QRIteration();
	void updateHessenbergOpMatrix_By_QT_IM_QRIteration();
	void updateHessenbergOpMatrix_By_Q_IM_QRIteration();
	//��ʽQR���� ��������ת������Q QT
	//void updateTotalQQT();
	void updateQT_Total();
	//void updateQ_Total();
protected:

	//p+t ����wilkinsonλ��ֵ�ĺ�
	double trace;

	//p*t ����wilkinsonλ��ֵ�ĳ˻�
	double determinant;

	// B = (A-pI)*(A-tI); ˫����λ��ľ����һ��
	double b11;
	double b21;
	double b31;

	//��������
	BasicMatrix* p_OpMatrix;

	//���ڼ�����ʽ˫��QR��������ת������Qi/Qn-1 ������---3ά/2ά
	BasicVector* p_TransVectorForQStep;

	//Q ���� ��ʽ����  ����Q�����ҳ�OP����
	//BasicMatrix* p_QMatrix_Implicit_Total;
	//QT ���� ��ʽ���� ����QT�������OP����
	BasicMatrix* p_QTMatrix_Implicit_Total;

	//Q/QT ���� ��ʽ���� �ֲ� QT�����ҳ�/���OP����
	BasicMatrix* p_QQTMatrix_Implicit_Step;

	//�м���̾���
	BasicMatrix* p_TempMatrix;


	//Hessenberg����
	//BasicMatrix* p_OpHessenbergMatrix;

	//���ڼ�����ʽ˫��QR��������ת������Qi������---3ά
	//BasicVector* p_TransVectorForQStep;
	//���ڼ�����ʽ˫��QR��������ת������Qn-1������---2ά
	//BasicVector* p_TransVectorForQ_LastStep;



	//Q ���� ��ʽ���� �ֲ� Q�����ҳ�OP����
	//BasicMatrix* p_QMatrix_Implicit_Step;
	//QT ���� ��ʽ���� �ֲ� QT�������OP����
	//BasicMatrix* p_QTMatrix_Implicit_Step;


	//Qi �Ӿ��� ��ʽ���� �ֲ� Qi�����ҳ�OP���� 3x3
	//BasicMatrix* p_QSubMatrix_Implicit_Step;

	//Qn-1 �Ӿ��� ��ʽ���� �ֲ� Qn-1�����ҳ�OP���� 2x2
	//BasicMatrix* p_QSubMatrix_Implicit_LastStep;

	//QTi �Ӿ��� ��ʽ���� �ֲ� QTi�������OP���� 3x3
	//BasicMatrix* p_QTSubMatrix_Implicit_Step;

	//QTn-1 �Ӿ��� ��ʽ���� �ֲ� QT�������OP���� 2x2
	//BasicMatrix* p_QTSubMatrix_Implicit_LastStep;

	//QR�ֽ�
	//QRDecomposition* p_QRDecomp;

	//hessen��ʽ��
	HessenbergFormular m_HessenbergForm;

	//Givens�任��
	GivensTransformation m_GivensTrans;

	//householder�任
	HouseholderTransformation m_HouseholderTrans;

	//�˷���
	MatrixMultiplier m_Multiplier;

	//ת����
	MatrixTransposer m_Transposer;

};



#endif /* EIGEN_BASIC_DOUBLESHIFTQRITERATION_H_ */
