/*
 * NormalEigenSolverForReal.h
 *
 *  Created on: 2017��5��20��
 *      Author: looke
 */

#ifndef EIGEN_BASIC_NORMALEIGENSOLVERFORREAL_H_
#define EIGEN_BASIC_NORMALEIGENSOLVERFORREAL_H_

#include "BasicMatrix.h"
#include "DoubleShiftQRIteration.h"
#include "SingleShiftQRIteration.h"
#include "HessenbergDeflation.h"

class NormalEigenSolverForReal
{
public:
	//NormalEigenSolverForReal();
	NormalEigenSolverForReal(BasicMatrix* p_input_OpMatrix,BasicVector* p_input_Vector,BasicMatrix* p_input_QT_Total,BasicMatrix* p_input_Q_Total,BasicMatrix* p_input_QQTMatrix_It,BasicMatrix* p_input_OpMatrix_deflated,BasicMatrix* p_TempMatrix_Trans,BasicMatrix* p_TempMatrix);

	void init(BasicMatrix* p_input_OpMatrix,BasicVector* p_input_Vector,BasicMatrix* p_input_QT_Total,BasicMatrix* p_input_Q_Total,BasicMatrix* p_input_QQTMatrix_It,BasicMatrix* p_input_OpMatrix_deflated,BasicMatrix* p_TempMatrix_Trans,BasicMatrix* p_TempMatrix);
	void reload(BasicMatrix* p_input_OpMatrix,BasicVector* p_input_Vector,BasicMatrix* p_input_QT_Total,BasicMatrix* p_input_Q_Total,BasicMatrix* p_input_QQTMatrix_It,BasicMatrix* p_input_OpMatrix_deflated,BasicMatrix* p_TempMatrix_Trans,BasicMatrix* p_TempMatrix);
	//virtual ~NormalEigenSolverForReal(){};

	//���Ҳ����½��׵�
	//bool findNewDeflationPoint();

	//���ݽ��׵� ���ɽ��׵�Hessenberg����
	void generateDeflatedHessenbergMatrix();

	//���ѽ��׵ı任���� ������Ϊȫ�ߴ�任����
	void upgradeDeflatedQQTMatrix();

	//���������任���� �ϲ���Ϊ����任����
	void updateQTMatrixTotal();
	void updateQMatrixTotal();

	//����ԭʼhessenberg��������
	void generateHessenbergOpMatrix();

	//����ԭʼhessenberg��������
	void updateHessenbergOpMatrix_By_QT();
	void updateHessenbergOpMatrix_By_Q();
	//��������ֵ
	void calcEigenValue();

	//��ʼ������ֵ������ؾ���
	void initEigenCalcMatrix();

	//��������ж��߼�
	bool hasFinishedIteration();

	//2x2�Խǿ��Ƿ�Ϊ��������ֵ�ж�
	bool isDiagonalBlockComplexEigen(BasicMatrix* p_Input_OpMatrix);

	//������Ϊ�Խǿ��Ժ����һ�����������Խ����ϵ�2x2�Խǿ���������ǻ�
	void lastStepIteration(int startIndex);

	//������������С
	void resizeMatrixForDeflation();
	//���ѽ��׵��м���� ������Ϊȫ�ߴ�任����
	void upgradeDeflatedTempMatrix();
	//BasicMatrix* getEigenValueMatrix();
	BasicMatrix* getOpHessenbergMatrix();

	//BasicMatrix* getQTMatrix_Iteration();
	//BasicMatrix* getQMatrix_Iteration();

	BasicMatrix* getQTMatrix_Total();
	BasicMatrix* getQMatrix_Total();

	BasicMatrix* getOpHessenbergMatrix_deflated();

	//BasicMatrix* getQTMatrix_Deflated_Iteration();
	//BasicMatrix* getQMatrix_Deflated_Iteration();

	//���Դ�ӡ��QT_Total * OP * Q
	//void showQTOPxQ();
protected:

	//�����������ָʾ����ʼ��Ϊ0
	int deflationStart;
	//�����յ�����ָʾ����ʼ��Ϊn-1
	int deflationEnd;

	//ԭʼ��������
	BasicMatrix* p_OpMatrix;

	BasicVector* p_TransVector;

	//ԭʼ���������任���� QTΪ��˾��� QΪ�ҳ˾���
	BasicMatrix* p_QTMatrix_Total;
	BasicMatrix* p_QMatrix_Total;

	//ԭʼ����ĵ������任���� QTΪ��˾��� QΪ�ҳ˾���
	BasicMatrix* p_QQTMatrix_Iteration;

	//�ѽ��׵� ����Hessenberg����
	BasicMatrix* p_OpHessenbergMatrix_deflated;

	//�м���̾���
	BasicMatrix* p_TempMatrix_Trans;
	//�м���̾���
	BasicMatrix* p_TempMatrix;


	//ԭʼ���������Ӧ������ֵ����
	//BasicMatrix* p_EigenValueMatrix;

	//ԭʼ����Hessenberg����
	//BasicMatrix* p_OpHessenbergMatrix;


	//ԭʼ����ĵ������任���� QTΪ��˾��� QΪ�ҳ˾���
	//BasicMatrix* p_QTMatrix_Iteration;
	//BasicMatrix* p_QMatrix_Iteration;





	//�ѽ��׾���ĵ������任����Q\QT
	//BasicMatrix* p_QTMatrix_Deflated_Iteration;
	//BasicMatrix* p_QMatrix_Deflated_Iteration;

	//���һ����ԶԽ���2x2�����Ĳ���
	//BasicMatrix* p_LastStepMatrix_2x2;


	//˫�ز�QR������
	DoubleShiftQRIteration m_DoubleShifeQR;

	//����QR������
	SingleShiftQRIteration m_SingleShifeQR;

	//hessen��ʽ��
	HessenbergFormular m_HessenbergForm;

	//�˷���
	MatrixMultiplier m_Multiplier;

	//ת����
	MatrixTransposer m_Transposer;


	//hessenberg���׵������
	HessenbergDeflation m_HessenbergDeflation;

	//��������ת�þ������ʱ���Զ���
	//BasicMatrix* p_testForTemp_nxn;
	//MatrixMultiplier* p_testMulti;
};



#endif /* EIGEN_BASIC_NORMALEIGENSOLVERFORREAL_H_ */
