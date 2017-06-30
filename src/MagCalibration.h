/*
 * MagCalibration.h
 *
 *  Created on: 2017年6月28日
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

	//生成校正参数
	bool generateCaliInfo();

	//软磁校正矩阵
	StaticMatrix m_Mag_SoftIron_CaliMatrix;

	//硬磁校正系数
	double mag_X_HardIron_Cali_Para;
	double mag_Y_HardIron_Cali_Para;
	double mag_Z_HardIron_Cali_Para;

protected:

	//10x10 op matrix
	BasicMatrix* p_OpMatrix_Original;

	StaticMatrix m_OpMatrix_A;
	StaticMatrix m_OpMatrix_B;

	StaticVector m_TransVector;

	//原始矩阵的 迭代变换矩阵 Q\Z
	StaticMatrix m_QMatrix_Step;
	StaticMatrix m_ZMatrix_Step;

	//已降阶的 操作Hessenberg矩阵
	StaticMatrix m_OpMatrix_Hessenberg_deflated;
	//已降阶的 操作Triangle矩阵
	StaticMatrix m_OpMatrix_Triangle_deflated;

	//中间过程矩阵
	StaticMatrix m_QZMatrix_Step;
	//中间过程矩阵
	StaticMatrix m_TempMatrix_Trans;
	//中间过程矩阵
	StaticMatrix m_TempMatrix;


	//广义特征值计算
	GeneralizedEigenSolverForReal m_GeneralizedEigenValueCalc;

	//广义特征向量计算
	GeneralizedEigenVectorCalcForReal m_GeneralizedEigenVectorCalc;

	//矩阵开方计算
	MatrixSquareRootSolver m_MatrixSquareRoot;
};

#endif /* MAGCALIBRATION_H_ */
