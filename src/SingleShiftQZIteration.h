/*
 * SingleShiftQZIteration.h
 *
 *  Created on: 2017年5月26日
 *      Author: looke
 */

#ifndef EIGEN_BASIC_SINGLESHIFTQZITERATION_H_
#define EIGEN_BASIC_SINGLESHIFTQZITERATION_H_

#include "BasicMatrix.h"
#include "MatrixMultiplier.h"
#include "HessenbergTriangleFormular.h"
#include "GivensTransformation.h"
#include "ABInverseCalculator.h"

class SingleShiftQZIteration
{
public:
	//SingleShiftQZIteration();
	SingleShiftQZIteration(BasicMatrix* p_input_OpMatrix_A, BasicMatrix* p_input_OpMatrix_B,BasicMatrix* p_input_QMatrix_Total,BasicMatrix* p_input_ZMatrix_Total,BasicMatrix* p_input_QZMatrix_Step,BasicMatrix* p_input_TempMatrix_Trans,BasicMatrix* p_input_TempMatrix);

	void init(BasicMatrix* p_input_OpMatrix_A, BasicMatrix* p_input_OpMatrix_B,BasicMatrix* p_input_QMatrix_Total,BasicMatrix* p_input_ZMatrix_Total,BasicMatrix* p_input_QZMatrix_Step,BasicMatrix* p_input_TempMatrix_Trans,BasicMatrix* p_input_TempMatrix);
	void reload(BasicMatrix* p_input_OpMatrix_A, BasicMatrix* p_input_OpMatrix_B,BasicMatrix* p_input_QMatrix_Total,BasicMatrix* p_input_ZMatrix_Total,BasicMatrix* p_input_QZMatrix_Step,BasicMatrix* p_input_TempMatrix_Trans,BasicMatrix* p_input_TempMatrix);
	//virtual ~SingleShiftQZIteration(){};

	//生成H-T操作矩阵对
	void generateHessenTriangleOpMatrix();

	//计算B^-1最后2行 用于计算单步位移值
	//void generateBinvLastTwoRow();

	//计算AB^-1 右下角元素值 用于单步位移
	//double generateABinvLastOne();

	//初始化隐式迭代
	void initForImplicitQZ(double input_ShiftValue);

	//单值位移QZ迭代 隐式单步
	void implicit_QZIteration_Step(double input_shiftValue);

	//单值位移QZ迭代 隐式
	void implicit_QZIteration(double input_shiftValue);

	//单值rayleigh商位移QZ迭代 隐式 多次
	void rayleigh_Quotient_IM_QZIteration(int iterateNum);

	//单值rayleigh商位移QZ迭代 隐式 单步
	void rayleigh_Quotient_IM_QZIteration_Step();

	//使用Q矩阵左乘H-T矩阵对
	void updateHTMatrixByQ();
	//使用Z矩阵右乘H-T矩阵对
	void updateHTMatrixByZ();

	//更新综合转换矩阵Q/Z Total(将Step合并入Total)
	void updateQMatrix_Total();
	void updateZMatrix_Total();
	//获取Q/Z 总体转换矩阵
	BasicMatrix* getQMatrix_Total();
	BasicMatrix* getZMatrix_Total();
protected:
	//B^-1 最后两行元素
	//double Binv_n_1_n_1;
	//double Binv_n_1_n;
	//double Binv_n_n;

	//原始操作矩阵A Hessenberg
	BasicMatrix* p_OpMatrix_A;
	//原始操作矩阵B Triangle
	BasicMatrix* p_OpMatrix_B;

	//Z 矩阵 隐式迭代 综合 Z用于右乘OP矩阵
	BasicMatrix* p_ZMatrix_Implicit_Total;
	//Q 矩阵 隐式迭代 综合 Q用于左乘OP矩阵
	BasicMatrix* p_QMatrix_Implicit_Total;

	//Q/Z 矩阵 隐式迭代 分步 Q用于左乘OP矩阵 Z用于右乘OP矩阵
	BasicMatrix* p_QZMatrix_Implicit_Step;

	//中间过程矩阵
	BasicMatrix* p_TempMatrix_Trans;
	//中间过程矩阵
	BasicMatrix* p_TempMatrix;

	//Hessenberg 操作矩阵A
	//BasicMatrix* p_OpMatrix_Hessenberg;
	//Triangle 操作矩阵B
	//BasicMatrix* p_OpMatrix_Triangle;

	//Z 矩阵 隐式迭代 分步 Z用于右乘OP矩阵
	//BasicMatrix* p_ZMatrix_Implicit_Step;
	//Q 矩阵 隐式迭代 分步 Q用于左乘OP矩阵
	//BasicMatrix* p_QMatrix_Implicit_Step;



	//hessenberg-triangle格式化
	HessenbergTriangleFormular m_HessenbergTriangleFormular;

	//QZ-triangle 0元移位
	//QZTriangleZeroChasing* p_QZTriangleZeroChasing;

	//Givens变换器
	GivensTransformation m_GivensTrans;

	//乘法器
	MatrixMultiplier m_Multiplier;

	//AB^-1 元素计算器
	ABInverseCalculator m_ABInvCalc;
};



#endif /* EIGEN_BASIC_SINGLESHIFTQZITERATION_H_ */
