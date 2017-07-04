/*
 * MagCalibration.cpp
 *
 *  Created on: 2017��6��30��
 *      Author: looke
 */
#include "MagCalibration.h"
//#include <iostream>
using namespace std;

MagCalibration::MagCalibration(BasicMatrix* p_input_OpMatrix):
m_GeneralizedEigenValueCalc(p_input_OpMatrix,p_input_OpMatrix,p_input_OpMatrix->getColumnVector(0),p_input_OpMatrix,p_input_OpMatrix,p_input_OpMatrix,p_input_OpMatrix,p_input_OpMatrix,p_input_OpMatrix,p_input_OpMatrix,p_input_OpMatrix,p_input_OpMatrix),
m_GeneralizedEigenVectorCalc(p_input_OpMatrix,p_input_OpMatrix,p_input_OpMatrix,p_input_OpMatrix,p_input_OpMatrix,p_input_OpMatrix),
m_MatrixSquareRoot(p_input_OpMatrix,p_input_OpMatrix,p_input_OpMatrix,p_input_OpMatrix,p_input_OpMatrix),
m_MatrixEqSolver(p_input_OpMatrix),
m_MatrixMulti(p_input_OpMatrix,p_input_OpMatrix,p_input_OpMatrix)
{
	this->init(p_input_OpMatrix);
};

void MagCalibration::init(BasicMatrix* p_input_OpMatrix)
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

	//����任����
	this->m_QMatrix_Total.resizeMatrix(p_OpMatrix_Original->rowNum,p_OpMatrix_Original->columnNum);
	m_QMatrix_Total.resetMatrixToI();

	this->m_ZMatrix_Total.resizeMatrix(p_OpMatrix_Original->rowNum,p_OpMatrix_Original->columnNum);
	m_ZMatrix_Total.resetMatrixToI();

	//�ѽ��׵� ����Hessenberg����
	this->m_OpMatrix_Hessenberg_deflated.resizeMatrix(p_OpMatrix_Original->rowNum,p_OpMatrix_Original->columnNum);
	m_OpMatrix_Hessenberg_deflated.resetMatrixToI();
	//�ѽ��׵� ����Triangle����
	this->m_OpMatrix_Triangle_deflated.resizeMatrix(p_OpMatrix_Original->rowNum,p_OpMatrix_Original->columnNum);
	m_OpMatrix_Triangle_deflated.resetMatrixToI();

	//ԭʼ����� �����任���� Q\Z
	this->m_QMatrix_Step.resizeMatrix(p_OpMatrix_Original->rowNum,p_OpMatrix_Original->columnNum);
	m_QMatrix_Step.resetMatrixToI();

	this->m_ZMatrix_Step.resizeMatrix(p_OpMatrix_Original->rowNum,p_OpMatrix_Original->columnNum);
	m_ZMatrix_Step.resetMatrixToI();

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

void MagCalibration::reload(BasicMatrix* p_input_OpMatrix)
{
	if(this->p_OpMatrix_Original->rowNum == p_input_OpMatrix->rowNum && this->p_OpMatrix_Original->columnNum == p_input_OpMatrix->columnNum)
	{
		//10x10 op matrix
		this->p_OpMatrix_Original = p_input_OpMatrix;

		m_OpMatrix_A.copyMatrixElementNoCheck(p_OpMatrix_Original);

		m_OpMatrix_B.resetMatrixToZero();
		m_OpMatrix_B.setMatrixElement(0,2,2);
		m_OpMatrix_B.setMatrixElement(1,1,-1);
		m_OpMatrix_B.setMatrixElement(2,0,2);

		//ת������
		this->m_TransVector.resetVectorElementToZero();

		//����任����
		m_QMatrix_Total.resetMatrixToI();
		m_ZMatrix_Total.resetMatrixToI();

		//�ѽ��׵� ����Hessenberg����
		m_OpMatrix_Hessenberg_deflated.resetMatrixToI();
		//�ѽ��׵� ����Triangle����
		m_OpMatrix_Triangle_deflated.resetMatrixToI();

		//ԭʼ����� �����任���� Q\Z
		m_QMatrix_Step.resetMatrixToI();
		m_ZMatrix_Step.resetMatrixToI();

		//�м���̾���
		m_QZMatrix_Step.resetMatrixToI();

		//�м���̾���
		m_TempMatrix_Trans.resetMatrixToI();

		//�м���̾���
		m_TempMatrix.resetMatrixToI();

		//���У������
		m_Mag_SoftIron_CaliMatrix.resetMatrixToI();

		//Ӳ��У��ϵ��
		this->mag_X_HardIron_Cali_Para = 0;
		this->mag_Y_HardIron_Cali_Para = 0;
		this->mag_Z_HardIron_Cali_Para = 0;
	}
	else
	{
		this->init(p_input_OpMatrix);
	}

};

//����У������
bool MagCalibration::generateCaliInfo()
{
	//�����������ֵ
	m_GeneralizedEigenValueCalc.reload(
			&m_OpMatrix_A,
			&m_OpMatrix_B,
			&m_TransVector,
			&m_OpMatrix_Hessenberg_deflated,
			&m_OpMatrix_Triangle_deflated,
			&m_QMatrix_Total,
			&m_ZMatrix_Total,
			&m_QMatrix_Step,
			&m_ZMatrix_Step,
			&m_QZMatrix_Step,
			&m_TempMatrix_Trans,
			&m_TempMatrix);
	m_GeneralizedEigenValueCalc.calcEigenValue();

	//����������ֵ�������Խ�������
	int indexOfPositiveEigenValue = m_GeneralizedEigenValueCalc.findPositiveEigenValue();

	//��������������
	m_GeneralizedEigenVectorCalc.reload(
			&m_OpMatrix_A,
			&m_OpMatrix_B,
			&m_QZMatrix_Step,
			&m_OpMatrix_Hessenberg_deflated,
			&m_OpMatrix_Triangle_deflated,
			&m_TempMatrix);

	bool result = m_GeneralizedEigenVectorCalc.getEigenVector(indexOfPositiveEigenValue, &m_TransVector);

	//����������������Ӳ�š����У������
	//v=a*x.^2+b*x.*y+c*y.^2+d*x.*z+e*y.*z+j*z.^2+p*x+q*y+r*z+s;
	if(result)
	{
		//����Z_total*��������
		StaticMatrix vectorForABInv = StaticMatrix(p_OpMatrix_Original->columnNum,1);
		for(int i=0;i<p_OpMatrix_Original->columnNum; i++)
		{
			vectorForABInv.setMatrixElement(i,0,m_TransVector.getElement(i));
		}

		m_TempMatrix.resizeMatrix(p_OpMatrix_Original->rowNum,1);
		//cout << "Z:" << endl;
		//m_ZMatrix_Total.printMatrix();
		//cout << "Vector for ABinv" << endl;
		//vectorForABInv.printMatrix();

		m_MatrixMulti.reload(&m_ZMatrix_Total,&vectorForABInv,&m_TempMatrix);
		m_MatrixMulti.multiplyCalc();

		for(int i=0;i<p_OpMatrix_Original->columnNum; i++)
		{
			m_TransVector.setElement(i,m_TempMatrix.getMatrixElement(i,0));
		}
		m_TransVector.normalizationVector();

		//��ʼ��10����Բ����
		/*
		a = m_TempMatrix.getMatrixElement(0,0);
		b = m_TempMatrix.getMatrixElement(1,0);
		c = m_TempMatrix.getMatrixElement(2,0);
		d = m_TempMatrix.getMatrixElement(3,0);
		e = m_TempMatrix.getMatrixElement(4,0);
		j = m_TempMatrix.getMatrixElement(5,0);
		p = m_TempMatrix.getMatrixElement(6,0);
		q = m_TempMatrix.getMatrixElement(7,0);
		r = m_TempMatrix.getMatrixElement(8,0);
		s = m_TempMatrix.getMatrixElement(9,0);
		*/
		a = m_TransVector.getElement(0);
		b = m_TransVector.getElement(1);
		c = m_TransVector.getElement(2);
		d = m_TransVector.getElement(3);
		e = m_TransVector.getElement(4);
		j = m_TransVector.getElement(5);
		p = m_TransVector.getElement(6);
		q = m_TransVector.getElement(7);
		r = m_TransVector.getElement(8);
		s = m_TransVector.getElement(9);

		//��ʼ������ϵ��
		double scale;
		scale = 4*a*c - b*b;
		scale = 1/scale;
		scale = sqrt(scale);

		updateEllipseParaByScale(scale);

		//Ӳ��У������
		bool hardResult = generateHardIronCaliInfo();

		//���У������
		bool softResult = generateSoftIronCaliInfo();
		if(!softResult)
		{
			updateEllipseParaByScale(-1);
			softResult = generateSoftIronCaliInfo();
		}
		result = hardResult&&softResult;
	}

	return result;
};

void MagCalibration::updateEllipseParaByScale(double scale)
{
	a = scale * a;
	b = scale * b;
	c = scale * c;
	d = scale * d;
	e = scale * e;
	j = scale * j;
	p = scale * p;
	q = scale * q;
	r = scale * r;
	s = scale * s;
};

//Ӳ��У��ϵ�����
bool MagCalibration::generateHardIronCaliInfo()
{
	//Ӳ��У�������ʼ��
	StaticMatrix hardIronMatrix_Cali_in = StaticMatrix(3,4);
	hardIronMatrix_Cali_in.setMatrixElement(0,0,0-2*a);
	hardIronMatrix_Cali_in.setMatrixElement(0,1,0-b);
	hardIronMatrix_Cali_in.setMatrixElement(0,2,0-d);
	hardIronMatrix_Cali_in.setMatrixElement(0,3,p);

	hardIronMatrix_Cali_in.setMatrixElement(1,0,0-b);
	hardIronMatrix_Cali_in.setMatrixElement(1,1,0-2*c);
	hardIronMatrix_Cali_in.setMatrixElement(1,2,0-e);
	hardIronMatrix_Cali_in.setMatrixElement(1,3,q);

	hardIronMatrix_Cali_in.setMatrixElement(2,0,0-d);
	hardIronMatrix_Cali_in.setMatrixElement(2,1,0-e);
	hardIronMatrix_Cali_in.setMatrixElement(2,2,0-2*j);
	hardIronMatrix_Cali_in.setMatrixElement(2,3,r);

	m_MatrixEqSolver.reload(&hardIronMatrix_Cali_in);
	bool result = m_MatrixEqSolver.solveMatrixEquation();
	if(result)
	{
		this->mag_X_HardIron_Cali_Para = hardIronMatrix_Cali_in.getMatrixElement(0,3);
		this->mag_Y_HardIron_Cali_Para = hardIronMatrix_Cali_in.getMatrixElement(1,3);
		this->mag_Z_HardIron_Cali_Para = hardIronMatrix_Cali_in.getMatrixElement(2,3);
	}
	return result;
};

//���У��ϵ�����
bool MagCalibration::generateSoftIronCaliInfo()
{

	m_Mag_SoftIron_CaliMatrix.setMatrixElement(0,0,a);
	m_Mag_SoftIron_CaliMatrix.setMatrixElement(0,1,0.5*b);
	m_Mag_SoftIron_CaliMatrix.setMatrixElement(0,2,0.5*d);

	m_Mag_SoftIron_CaliMatrix.setMatrixElement(1,0,0.5*b);
	m_Mag_SoftIron_CaliMatrix.setMatrixElement(1,1,c);
	m_Mag_SoftIron_CaliMatrix.setMatrixElement(1,2,0.5*e);

	m_Mag_SoftIron_CaliMatrix.setMatrixElement(2,0,0.5*d);
	m_Mag_SoftIron_CaliMatrix.setMatrixElement(2,1,0.5*e);
	m_Mag_SoftIron_CaliMatrix.setMatrixElement(2,2,j);

	bool result;
	m_ZMatrix_Step.resizeMatrix(3,3);
	m_QZMatrix_Step.resizeMatrix(3,3);
	m_TempMatrix_Trans.resizeMatrix(3,3);
	m_TempMatrix.resizeMatrix(3,3);

	m_MatrixSquareRoot.reload(&m_Mag_SoftIron_CaliMatrix, &m_ZMatrix_Step, &m_QZMatrix_Step, &m_TempMatrix_Trans, &m_TempMatrix);

	result = m_MatrixSquareRoot.generateSquareRootMatrix();

	//������ʧ�ܣ������У����������Ϊ��λ����
	if(!result)
	{
		m_Mag_SoftIron_CaliMatrix.resetMatrixToI();
	}

	m_ZMatrix_Step.resizeMatrix(p_OpMatrix_Original->rowNum,p_OpMatrix_Original->columnNum);
	m_QZMatrix_Step.resizeMatrix(p_OpMatrix_Original->rowNum,p_OpMatrix_Original->columnNum);
	m_TempMatrix_Trans.resizeMatrix(p_OpMatrix_Original->rowNum,p_OpMatrix_Original->columnNum);
	m_TempMatrix.resizeMatrix(p_OpMatrix_Original->rowNum,p_OpMatrix_Original->columnNum);

	return result;
};
