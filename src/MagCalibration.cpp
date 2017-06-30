/*
 * MagCalibration.cpp
 *
 *  Created on: 2017��6��30��
 *      Author: looke
 */
#include "MagCalibration.h"

MagCalibration::MagCalibration(BasicMatrix* p_input_OpMatrix):
m_GeneralizedEigenValueCalc(p_input_OpMatrix,p_input_OpMatrix,p_input_OpMatrix->getColumnVector(0),p_input_OpMatrix,p_input_OpMatrix,p_input_OpMatrix,p_input_OpMatrix,p_input_OpMatrix,p_input_OpMatrix,p_input_OpMatrix,p_input_OpMatrix,p_input_OpMatrix),
m_GeneralizedEigenVectorCalc(p_input_OpMatrix,p_input_OpMatrix,p_input_OpMatrix,p_input_OpMatrix,p_input_OpMatrix,p_input_OpMatrix),
m_MatrixSquareRoot(p_input_OpMatrix,p_input_OpMatrix,p_input_OpMatrix,p_input_OpMatrix,p_input_OpMatrix)
{
	//10x10 op matrix
	this->p_OpMatrix_Original = p_input_OpMatrix;

	this->m_OpMatrix_A.resizeMatrix(p_OpMatrix_Original->rowNum,p_OpMatrix_Original->columnNum);
	m_OpMatrix_A.copyMatrixElementNoCheck(p_OpMatrix_Original);

	this->m_OpMatrix_B.resizeMatrix(p_OpMatrix_Original->rowNum,p_OpMatrix_Original->columnNum);
	m_OpMatrix_B.resetMatrixToZero();
	m_OpMatrix_B.setMatrixElement(0,2,2);
	m_OpMatrix_B.setMatrixElement(1,1,-1);
	m_OpMatrix_B.setMatrixElement(2,0,2);

	//ת������
	this->m_TransVector.resetDimension(p_OpMatrix_Original->rowNum);
	this->m_TransVector.resetVectorElementToZero();

	//ԭʼ����� �����任���� Q\Z
	this->m_QMatrix_Step.resizeMatrix(p_OpMatrix_Original->rowNum,p_OpMatrix_Original->columnNum);
	m_QMatrix_Step.resetMatrixToI();

	this->m_ZMatrix_Step.resizeMatrix(p_OpMatrix_Original->rowNum,p_OpMatrix_Original->columnNum);
	m_ZMatrix_Step.resetMatrixToI();

	//�ѽ��׵� ����Hessenberg����
	this->m_OpMatrix_Hessenberg_deflated.resizeMatrix(p_OpMatrix_Original->rowNum,p_OpMatrix_Original->columnNum);
	m_OpMatrix_Hessenberg_deflated.resetMatrixToI();
	//�ѽ��׵� ����Triangle����
	this->m_OpMatrix_Triangle_deflated.resizeMatrix(p_OpMatrix_Original->rowNum,p_OpMatrix_Original->columnNum);
	m_OpMatrix_Triangle_deflated.resetMatrixToI();

	//�м���̾���
	this->m_QZMatrix_Step.resizeMatrix(p_OpMatrix_Original->rowNum,p_OpMatrix_Original->columnNum);
	m_QZMatrix_Step.resetMatrixToI();

	//�м���̾���
	this->m_TempMatrix_Trans.resizeMatrix(p_OpMatrix_Original->rowNum,p_OpMatrix_Original->columnNum);
	m_TempMatrix_Trans.resetMatrixToI();

	//�м���̾���
	this->m_TempMatrix.resizeMatrix(p_OpMatrix_Original->rowNum,p_OpMatrix_Original->columnNum);
	m_TempMatrix.resetMatrixToI();


	//���У������
	this->m_Mag_SoftIron_CaliMatrix.resizeMatrix(3,3);
	m_Mag_SoftIron_CaliMatrix.resetMatrixToI();

	//Ӳ��У��ϵ��
	this->mag_X_HardIron_Cali_Para = 0;
	this->mag_Y_HardIron_Cali_Para = 0;
	this->mag_Z_HardIron_Cali_Para = 0;
};



