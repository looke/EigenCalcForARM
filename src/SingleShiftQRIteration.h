/*
 * SingelShiftQRIteration.h
 *
 *  Created on: 2017年5月13日
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

	//瑞利商位移QR迭代 显式
	//void rayleigh_Quotient_EX_QRIteration();

	//瑞利商位移QR迭代 隐式 10次迭代
	//void rayleigh_Quotient_IM_QRIteration();
	//瑞利商位移QR迭代 隐式 指定迭代次数
	void rayleigh_Quotient_IM_QRIteration(int iterateNum);
	//瑞利商位移QR迭代 隐式 单步
	void rayleigh_Quotient_IM_QRIteration_Step();

	//单值位移QR迭代 显式
	//void explicit_QRIteration(double input_shiftValue);
	//单值位移QR迭代 显式单步
	//void explicit_QRIteration_Step(double input_shiftValue);

	//单值位移QR迭代 隐式
	void implicit_QRIteration(double input_shiftValue);
	//单值位移QR迭代 隐式单步
	void implicit_QRIteration_Step(double input_shiftValue);

	//单值位移QR迭代 隐式 初始化
	void initForImplicitQR(double input_shiftValue);

	//更新综合转换矩阵QT Total
	void updateQTMatrix_Total_IM_QRIteration();
	//更新综合转换矩阵Q Total
	//void updateQMatrix_Total_IM_QRIteration();

	//生成hessenberg操作矩阵
	void generateHessenbergOpMatrix();

	BasicMatrix* getOpHessenbergMatrix();

	//获取总体转换矩阵
	BasicMatrix* getQTMatrix_Total();
	//BasicMatrix* getQMatrix_Total();

protected:

	//隐式QR迭代 更新hessenberg操作矩阵
	//void updateHessenbergOpMatrix_IM_QRIteration();
	//隐式QR迭代 更新操作矩阵
	void updateOpMatrix_By_QT_IM_QRIteration();
	void updateOpMatrix_By_Q_IM_QRIteration();
	//操作矩阵
	BasicMatrix* p_OpMatrix;

	//Q 矩阵 隐式迭代 综合 Q用于右乘OP矩阵
	//BasicMatrix* p_QMatrix_Implicit_Total;
	//QT 矩阵 隐式迭代 综合 QT用于左乘OP矩阵
	BasicMatrix* p_QTMatrix_Implicit_Total;

	//Hessenberg矩阵
	//BasicMatrix* p_OpHessenbergMatrix;

	//Q 矩阵 显式迭代
	//BasicMatrix* p_QMatrix_Explicit;

	//Q 矩阵 隐式迭代 分步 Q用于右乘OP矩阵
	//BasicMatrix* p_QMatrix_Implicit_Step;
	//QT 矩阵 隐式迭代 分步 QT用于左乘OP矩阵
	//BasicMatrix* p_QTMatrix_Implicit_Step;

	//Q/QT 矩阵 隐式迭代 分步 Q用于右乘OP矩阵 QT用于左乘OP矩阵
	BasicMatrix* p_Q_QT_Matrix_Implicit_Step;

	//中间过程矩阵
	BasicMatrix* p_TempMatrix;


	//QR分解
	QRDecomposition m_QRDecomp;

	//hessen格式化
	HessenbergFormular m_HessenbergForm;

	//Givens变换器
	GivensTransformation m_GivensTrans;

	//乘法器
	MatrixMultiplier m_Multiplier;

	//转置器
	MatrixTransposer m_Transposer;
};



#endif /* EIGEN_BASIC_SINGLESHIFTQRITERATION_H_ */
