/*
 * MagCalibration.h
 *
 *  Created on: 2017��6��28��
 *      Author: looke
 */

#ifndef MAGCALIBRATION_H_
#define MAGCALIBRATION_H_

#include "GeneralizedEigenSolverForReal.h"
#include "GeneralizedEigenVectorCalcForReal.h"
#include "MatrixSquareRootSolver.h"
#include "MatrixEquationSolver.h"

class MagCalibration
{
public:
	MagCalibration(BasicMatrix* p_input_OpMatrix);

	void init(BasicMatrix* p_input_OpMatrix);
	void reload(BasicMatrix* p_input_OpMatrix);

	//����У������
	bool generateCaliInfo();

	//Mag������Բ����
	//v=a*x.^2+b*x.*y+c*y.^2+d*x.*z+e*y.*z+j*z.^2+p*x+q*y+r*z+s;
	double a,b,c,d,e,j,p,q,r,s;
	//Mag������Բ��������ϵ��
	//double scale;
	//���У������
	StaticMatrix m_Mag_SoftIron_CaliMatrix;

	//Ӳ��У��ϵ��
	double mag_X_HardIron_Cali_Para;
	double mag_Y_HardIron_Cali_Para;
	double mag_Z_HardIron_Cali_Para;

protected:

	//10x10 op matrix
	BasicMatrix* p_OpMatrix_Original;

	StaticMatrix m_OpMatrix_A;
	StaticMatrix m_OpMatrix_B;

	StaticVector m_TransVector;


	//�ѽ��׵� ����Hessenberg����
	StaticMatrix m_OpMatrix_Hessenberg_deflated;
	//�ѽ��׵� ����Triangle����
	StaticMatrix m_OpMatrix_Triangle_deflated;

	//����任���� Q\Z Total
	StaticMatrix m_QMatrix_Total;
	StaticMatrix m_ZMatrix_Total;


	//ԭʼ����� �����任���� Q\Z
	StaticMatrix m_QMatrix_Step;
	StaticMatrix m_ZMatrix_Step;

	//�м���̾���
	StaticMatrix m_QZMatrix_Step;
	//�м���̾���
	StaticMatrix m_TempMatrix_Trans;
	//�м���̾���
	StaticMatrix m_TempMatrix;

	//��������ֵ����
	GeneralizedEigenSolverForReal m_GeneralizedEigenValueCalc;

	//����������������
	GeneralizedEigenVectorCalcForReal m_GeneralizedEigenVectorCalc;

	//���󿪷�����
	MatrixSquareRootSolver m_MatrixSquareRoot;

	//�������Է��������
	MatrixEquationSolver m_MatrixEqSolver;

	//����˻�
	MatrixMultiplier m_MatrixMulti;

	void updateEllipseParaByScale(double scale);

	//Ӳ��У��ϵ�����
	bool generateHardIronCaliInfo();
	//���У��ϵ�����
	bool generateSoftIronCaliInfo();

};

#endif /* MAGCALIBRATION_H_ */
