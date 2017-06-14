/*
 * SingelShiftQRIteration.h
 *
 *  Created on: 2017��5��13��
 *      Author: looke
 */

#ifndef EIGEN_BASIC_SINGLESHIFTQRITERATION_H_
#define EIGEN_BASIC_SINGLESHIFTQRITERATION_H_

#include "BasicMatrix.h"
#include "MatrixTransposer.h"
#include "QRDecomposition.h"
#include "HessenbergFormular.h"
#include "GivensTransformation.h"


class SingleShiftQRIteration
{
public:
	//SingleShiftQRIteration();
	SingleShiftQRIteration(BasicMatrix* p_input_OpMatrix, BasicMatrix* p_input_QTMatrix_Implicit_Total,BasicMatrix* p_input_Q_QT_Matrix_Implicit_Step,BasicMatrix* p_input_TempMatrix);
	void init(BasicMatrix* p_input_OpMatrix,BasicMatrix* p_input_QTMatrix_Implicit_Total,BasicMatrix* p_input_Q_QT_Matrix_Implicit_Step,BasicMatrix* p_input_TempMatrix);
	void reload(BasicMatrix* p_input_OpMatrix, BasicMatrix* p_input_QTMatrix_Implicit_Total,BasicMatrix* p_input_Q_QT_Matrix_Implicit_Step,BasicMatrix* p_input_TempMatrix);
	//virtual ~SingleShiftQRIteration(){};

	//������λ��QR���� ��ʽ
	//void rayleigh_Quotient_EX_QRIteration();

	//������λ��QR���� ��ʽ 10�ε���
	void rayleigh_Quotient_IM_QRIteration();
	//������λ��QR���� ��ʽ ָ����������
	void rayleigh_Quotient_IM_QRIteration(int iterateNum);
	//������λ��QR���� ��ʽ ����
	void rayleigh_Quotient_IM_QRIteration_Step();

	//��ֵλ��QR���� ��ʽ
	//void explicit_QRIteration(double input_shiftValue);
	//��ֵλ��QR���� ��ʽ����
	//void explicit_QRIteration_Step(double input_shiftValue);

	//��ֵλ��QR���� ��ʽ
	void implicit_QRIteration(double input_shiftValue);
	//��ֵλ��QR���� ��ʽ����
	void implicit_QRIteration_Step(double input_shiftValue);

	//��ֵλ��QR���� ��ʽ ��ʼ��
	void initForImplicitQR(double input_shiftValue);

	//�����ۺ�ת������QT Total
	void updateQTMatrix_Total_IM_QRIteration();
	//�����ۺ�ת������Q Total
	//void updateQMatrix_Total_IM_QRIteration();

	//����hessenberg��������
	void generateHessenbergOpMatrix();

	BasicMatrix* getOpHessenbergMatrix();

	//��ȡ����ת������
	BasicMatrix* getQTMatrix_Total();
	//BasicMatrix* getQMatrix_Total();

protected:

	//��ʽQR���� ����hessenberg��������
	//void updateHessenbergOpMatrix_IM_QRIteration();
	//��ʽQR���� ���²�������
	void updateOpMatrix_By_QT_IM_QRIteration();
	void updateOpMatrix_By_Q_IM_QRIteration();
	//��������
	BasicMatrix* p_OpMatrix;

	//Q ���� ��ʽ���� �ۺ� Q�����ҳ�OP����
	//BasicMatrix* p_QMatrix_Implicit_Total;
	//QT ���� ��ʽ���� �ۺ� QT�������OP����
	BasicMatrix* p_QTMatrix_Implicit_Total;

	//Hessenberg����
	//BasicMatrix* p_OpHessenbergMatrix;

	//Q ���� ��ʽ����
	//BasicMatrix* p_QMatrix_Explicit;

	//Q ���� ��ʽ���� �ֲ� Q�����ҳ�OP����
	//BasicMatrix* p_QMatrix_Implicit_Step;
	//QT ���� ��ʽ���� �ֲ� QT�������OP����
	//BasicMatrix* p_QTMatrix_Implicit_Step;

	//Q/QT ���� ��ʽ���� �ֲ� Q�����ҳ�OP���� QT�������OP����
	BasicMatrix* p_Q_QT_Matrix_Implicit_Step;

	//�м���̾���
	BasicMatrix* p_TempMatrix;


	//QR�ֽ�
	QRDecomposition m_QRDecomp;

	//hessen��ʽ��
	HessenbergFormular m_HessenbergForm;

	//Givens�任��
	GivensTransformation m_GivensTrans;

	//�˷���
	MatrixMultiplier m_Multiplier;

	//ת����
	MatrixTransposer m_Transposer;
};



#endif /* EIGEN_BASIC_SINGLESHIFTQRITERATION_H_ */
