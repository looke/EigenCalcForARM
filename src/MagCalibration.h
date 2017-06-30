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
class MagCalibration
{
public:
	MagCalibration(BasicMatrix* p_input_OpMatrix);

	void init(BasicMatrix* p_input_OpMatrix);
	void reload(BasicMatrix* p_input_OpMatrix);

	//����У������
	bool generateCaliInfo();

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

	//ԭʼ����� �����任���� Q\Z
	StaticMatrix m_QMatrix_Step;
	StaticMatrix m_ZMatrix_Step;

	//�ѽ��׵� ����Hessenberg����
	StaticMatrix m_OpMatrix_Hessenberg_deflated;
	//�ѽ��׵� ����Triangle����
	StaticMatrix m_OpMatrix_Triangle_deflated;

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
};

#endif /* MAGCALIBRATION_H_ */
