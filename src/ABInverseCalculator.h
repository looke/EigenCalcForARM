/*
 * ABInverseCalculator.h
 *
 *  Created on: 2017��6��3��
 *      Author: looke
 */

#ifndef EIGEN_UTIL_ABINVERSECALCULATOR_H_
#define EIGEN_UTIL_ABINVERSECALCULATOR_H_

#include "BasicMatrix.h"
class ABInverseCalculator
{
public:
	ABInverseCalculator();

	//����B^-1�������
	void generateBinvLastThreeRow(BasicMatrix* p_OpMatrix_Tirangle);

	//����B^-1������
	void generateBinvLastTwoRow(BasicMatrix* p_OpMatrix_Triangle);

	//���� A*B^-1 ����2x2�Ӿ���Ԫ��
	void generateABinvLast2x2(BasicMatrix* p_OpMatrix_Hessenberg, BasicMatrix* p_OpMatrix_Triangle);

	//���� A*B^-1 ǰ����Ԫ��
	void generateABinvFirst3x2(BasicMatrix* p_OpMatrix_Hessenberg, BasicMatrix* p_OpMatrix_Triangle);

	//���� A*B^-1 ����2x2�Ӿ���
	void generateABinvFirst2x2(BasicMatrix* p_OpMatrix_Hessenberg, BasicMatrix* p_OpMatrix_Triangle);

	//����AB^-1 ���½�Ԫ��ֵ ���ڵ���λ��
	void generateABinvLastOne(BasicMatrix* p_OpMatrix_Hessenberg, BasicMatrix* p_OpMatrix_Triangle);

	double getABinv11() const
	{
		return ABinv_11;
	}

	double getABinv12() const
	{
		return ABinv_12;
	}

	double getABinv21() const
	{
		return ABinv_21;
	}

	double getABinv22() const
	{
		return ABinv_22;
	}

	double getABinv32() const
	{
		return ABinv_32;
	}

	double getABinv_N_1_N() const
	{
		return ABinv_n_1_n;
	}

	double getABinv_N_1_N_1() const
	{
		return ABinv_n_1_n_1;
	}

	double getABinv_N_N() const
	{
		return ABinv_n_n;
	}

	double getABinv_N_N_1() const
	{
		return ABinv_n_n_1;
	}

	double getBinv_N_1_N() const
	{
		return Binv_n_1_n;
	}

	double getBinv_N_1_N_1() const
	{
		return Binv_n_1_n_1;
	}

	double getBinv_N_2_N() const
	{
		return Binv_n_2_n;
	}

	double getBinv_N_2_N_1() const
	{
		return Binv_n_2_n_1;
	}

	double getBinv_N_2_N_2() const
	{
		return Binv_n_2_n_2;
	}

	double getBinv_N_N() const
	{
		return Binv_n_n;
	}

protected:

	//B^-1 ���3��
	double Binv_n_2_n_2;
	double Binv_n_2_n_1;
	double Binv_n_2_n;
	double Binv_n_1_n_1;
	double Binv_n_1_n;
	double Binv_n_n;

	//A*B^-1 ����2x2�Ӿ���Ԫ��
	double ABinv_n_1_n_1;
	double ABinv_n_1_n;
	double ABinv_n_n_1;
	double ABinv_n_n;

	//A*B^-1 ǰ����Ԫ��
	double ABinv_11;
	double ABinv_12;
	double ABinv_21;
	double ABinv_22;
	double ABinv_32;
};

#endif /* EIGEN_UTIL_ABINVERSECALCULATOR_H_ */
