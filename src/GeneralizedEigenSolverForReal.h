/*
 * GeneralizedEigenSolverForReal.h
 *
 *  Created on: 2017��6��1��
 *      Author: looke
 */

#ifndef EIGEN_BASIC_GENERALIZEDEIGENSOLVERFORREAL_H_
#define EIGEN_BASIC_GENERALIZEDEIGENSOLVERFORREAL_H_

#include "BasicMatrix.h"
#include "StaticMatrix.h"
#include "BasicVector.h"
#include "DoubleShiftQZIteration.h"
#include "SingleShiftQZIteration.h"
#include "HessenbergDeflation.h"
#include "QZTriangleZeroChasing.h"

class GeneralizedEigenSolverForReal
{
public:
	//GeneralizedEigenSolverForReal();
	GeneralizedEigenSolverForReal(BasicMatrix* p_input_OpMatrix_A, BasicMatrix* p_input_OpMatrix_B, BasicVector* p_input_Vector, BasicMatrix* p_input_A_deflated, BasicMatrix* p_input_B_deflated, BasicMatrix* p_input_Q_Total, BasicMatrix* p_input_Z_Total, BasicMatrix* p_input_Q_Step, BasicMatrix* p_input_Z_Step, BasicMatrix* p_input_QZ_Step, BasicMatrix* p_input_TempMatrix_Trans,BasicMatrix* p_input_TempMatrix);

	void init(BasicMatrix* p_input_OpMatrix_A, BasicMatrix* p_input_OpMatrix_B, BasicVector* p_input_Vector, BasicMatrix* p_input_A_deflated, BasicMatrix* p_input_B_deflated, BasicMatrix* p_input_Q_Total, BasicMatrix* p_input_Z_Total, BasicMatrix* p_input_Q_Step, BasicMatrix* p_input_Z_Step, BasicMatrix* p_input_QZ_Step, BasicMatrix* p_input_TempMatrix_Trans,BasicMatrix* p_input_TempMatrix);
	void reload(BasicMatrix* p_input_OpMatrix_A, BasicMatrix* p_input_OpMatrix_B, BasicVector* p_input_Vector, BasicMatrix* p_input_A_deflated, BasicMatrix* p_input_B_deflated, BasicMatrix* p_input_Q_Total, BasicMatrix* p_input_Z_Total, BasicMatrix* p_input_Q_Step, BasicMatrix* p_input_Z_Step, BasicMatrix* p_input_QZ_Step, BasicMatrix* p_input_TempMatrix_Trans,BasicMatrix* p_input_TempMatrix);
	//virtual ~GeneralizedEigenSolverForReal(){};

	//����ԭʼhessenberg-Triangle��������
	void generateHTOpMatrix();

	//���ݽ��׵� ���ɽ��׵�Hessenberg-Tirangle����
	void generateDeflatedHTMatrixPair();

	//���ѽ��׵ı任���� ������Ϊȫ�ߴ�任����
	void upgradeDeflatedQMatrix();
	void upgradeDeflatedZMatrix();

	//���������任���� �ϲ���Ϊ����任����
	void updateQMatrixTotal();
	void updateZMatrixTotal();

	//ʹ��Q�������H-T�����
	void updateHTMatrixByQ();
	//ʹ��Z�����ҳ�H-T�����
	void updateHTMatrixByZ();

	//��������ֵ
	void calcEigenValue();

	//��ʼ������ֵ������ؾ���
	void initEigenCalcMatrix();

	//��������ж��߼�
	bool hasFinishedIteration();

	//H-T����� ����2x2�Խǿ��Ƿ�Ϊ��������ֵ�ж�
	bool isDiagonalBlockComplexEigen(BasicMatrix* p_Input_OpMatrix_A, BasicMatrix* p_Input_OpMatrix_B);

	//������Ϊ�Խǿ��Ժ����һ�����������Խ����ϵ�2x2�Խǿ���������ǻ�
	void lastStepIteration(int startIndex);

	//����resize ��С������ؾ���
	void resizeTransMatrix();

	//�����ת����������Ϊȫά��
	void upgradTransMatrix();

	//����������ֵ�������Խ���������������ڶ��������ֵ��ȡ������С��
	int findPositiveEigenValue();
	//��ȡH-T�����
	//BasicMatrix* getHessenbergMatrix();
	//BasicMatrix* getTriangleMatrix();

	//��ȡQ\Z �ۺ�ת������
	//BasicMatrix* getQMatrix_Total();
	//BasicMatrix* getZMatrix_Total();

	//���Դ�ӡ��Q_Total * OP * Z_Total
	void showQxOPxZ();

protected:
	//�����������ָʾ����ʼ��Ϊ0
	int deflationStart;
	//�����յ�����ָʾ����ʼ��Ϊn-1
	int deflationEnd;
	//����������ޣ����������г������޶���Ĵ�����Ȼû�н���ʱ�����������������
	int noDeflateLimite;
	//ԭʼ��������
	BasicMatrix* p_OpMatrix_A;
	BasicMatrix* p_OpMatrix_B;

	BasicVector* p_OpTransVector;

	//ԭʼ�����ȫά������任����
	BasicMatrix* p_QMatrix_Total;
	BasicMatrix* p_ZMatrix_Total;

	//ԭʼ����Hessenberg����
	//BasicMatrix* p_OpMatrix_Hessenberg;
	//ԭʼ����Triangle����
	//BasicMatrix* p_OpMatrix_Triangle;

	//ԭʼ�����ȫά�ȵ������任����
	//BasicMatrix* p_QMatrix_Iteration;
	//BasicMatrix* p_ZMatrix_Iteration;

	//ԭʼ����� �����任���� Q\Z
	BasicMatrix* p_QMatrix_Step;
	BasicMatrix* p_ZMatrix_Step;

	//�ѽ��׵� ����Hessenberg����
	BasicMatrix* p_OpMatrix_Hessenberg_deflated;
	//�ѽ��׵� ����Triangle����
	BasicMatrix* p_OpMatrix_Triangle_deflated;

	//�ѽ��׾���ĵ������任����Q\Z
	//BasicMatrix* p_QMatrix_Deflated_Iteration;
	//BasicMatrix* p_ZMatrix_Deflated_Iteration;

	//���һ����ԶԽ���2x2�����Ĳ���
	//BasicMatrix* p_LastStepMatrix_2x2_Hessenberg;
	//BasicMatrix* p_LastStepMatrix_2x2_Triangle;


	//�м���̾���
	BasicMatrix* p_QZMatrix_Step;
	//�м���̾���
	BasicMatrix* p_TempMatrix_Trans;
	//�м���̾���
	BasicMatrix* p_TempMatrix;

	//˫�ز�QZ������
	DoubleShiftQZIteration m_DoubleShifeQZ;

	//����QR������
	SingleShiftQZIteration m_SingleShifeQZ;

	//hessenberg-triangle��ʽ��
	HessenbergTriangleFormular m_HessenbergTriangleFormular;

	//�˷���
	MatrixMultiplier m_Multiplier;

	//hessenberg���׵������
	HessenbergDeflation m_HessenbergDeflation;

	//QZ-triangle 0Ԫ��λ
	QZTriangleZeroChasing m_QZTriangleZeroChasing;

	//AB^-1 Ԫ�ؼ�����
	ABInverseCalculator m_ABInvCalc;

	//��������ת�þ������ʱ���Զ���
	StaticMatrix testForTemp_A_nxn;
	StaticMatrix testForTemp_B_nxn;
	StaticMatrix testTemp_nxn;
	MatrixMultiplier m_testMulti;

};



#endif /* EIGEN_BASIC_GENERALIZEDEIGENSOLVERFORREAL_H_ */
