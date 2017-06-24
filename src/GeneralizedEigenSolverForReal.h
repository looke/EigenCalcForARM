/*
 * GeneralizedEigenSolverForReal.h
 *
 *  Created on: 2017年6月1日
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

	//生成原始hessenberg-Triangle操作矩阵
	void generateHTOpMatrix();

	//根据降阶点 生成降阶的Hessenberg-Tirangle矩阵
	void generateDeflatedHTMatrixPair();

	//将已降阶的变换矩阵 升级成为全尺寸变换矩阵
	void upgradeDeflatedQMatrix();
	void upgradeDeflatedZMatrix();

	//将单迭代变换矩阵 合并成为总体变换矩阵
	void updateQMatrixTotal();
	void updateZMatrixTotal();

	//使用Q矩阵左乘H-T矩阵对
	void updateHTMatrixByQ();
	//使用Z矩阵右乘H-T矩阵对
	void updateHTMatrixByZ();

	//计算特征值
	void calcEigenValue();

	//初始化特征值计算相关矩阵
	void initEigenCalcMatrix();

	//解算结束判断逻辑
	bool hasFinishedIteration();

	//H-T矩阵对 左上2x2对角块是否为复数特征值判断
	bool isDiagonalBlockComplexEigen(BasicMatrix* p_Input_OpMatrix_A, BasicMatrix* p_Input_OpMatrix_B);

	//初步化为对角块以后，最后一步迭代，将对角线上的2x2对角块进行上三角化
	void lastStepIteration(int startIndex);

	//重新resize 缩小计算相关矩阵
	void resizeTransMatrix();

	//将相关转换矩阵升级为全维度
	void upgradTransMatrix();
	//获取H-T矩阵对
	BasicMatrix* getHessenbergMatrix();
	BasicMatrix* getTriangleMatrix();

	//获取Q\Z 综合转换矩阵
	BasicMatrix* getQMatrix_Total();
	BasicMatrix* getZMatrix_Total();

	//测试打印，Q_Total * OP * Z_Total
	void showQxOPxZ();

protected:
	//降阶起点索引指示，初始化为0
	int deflationStart;
	//降阶终点索引指示，初始化为n-1
	int deflationEnd;

	//原始操作矩阵
	BasicMatrix* p_OpMatrix_A;
	BasicMatrix* p_OpMatrix_B;

	BasicVector* p_OpTransVector;

	//原始矩阵的全维度总体变换矩阵
	BasicMatrix* p_QMatrix_Total;
	BasicMatrix* p_ZMatrix_Total;

	//原始操作Hessenberg矩阵
	//BasicMatrix* p_OpMatrix_Hessenberg;
	//原始操作Triangle矩阵
	//BasicMatrix* p_OpMatrix_Triangle;

	//原始矩阵的全维度单迭代变换矩阵
	//BasicMatrix* p_QMatrix_Iteration;
	//BasicMatrix* p_ZMatrix_Iteration;

	//原始矩阵的 迭代变换矩阵 Q\Z
	BasicMatrix* p_QMatrix_Step;
	BasicMatrix* p_ZMatrix_Step;

	//已降阶的 操作Hessenberg矩阵
	BasicMatrix* p_OpMatrix_Hessenberg_deflated;
	//已降阶的 操作Triangle矩阵
	BasicMatrix* p_OpMatrix_Triangle_deflated;

	//已降阶矩阵的单迭代变换矩阵Q\Z
	//BasicMatrix* p_QMatrix_Deflated_Iteration;
	//BasicMatrix* p_ZMatrix_Deflated_Iteration;

	//最后一步针对对角线2x2矩阵块的操作
	//BasicMatrix* p_LastStepMatrix_2x2_Hessenberg;
	//BasicMatrix* p_LastStepMatrix_2x2_Triangle;


	//中间过程矩阵
	BasicMatrix* p_QZMatrix_Step;
	//中间过程矩阵
	BasicMatrix* p_TempMatrix_Trans;
	//中间过程矩阵
	BasicMatrix* p_TempMatrix;

	//双重步QZ迭代器
	DoubleShiftQZIteration m_DoubleShifeQZ;

	//单步QR迭代器
	SingleShiftQZIteration m_SingleShifeQZ;

	//hessenberg-triangle格式化
	HessenbergTriangleFormular m_HessenbergTriangleFormular;

	//乘法器
	MatrixMultiplier m_Multiplier;

	//hessenberg降阶点查找器
	HessenbergDeflation m_HessenbergDeflation;

	//QZ-triangle 0元移位
	QZTriangleZeroChasing m_QZTriangleZeroChasing;

	//AB^-1 元素计算器
	ABInverseCalculator m_ABInvCalc;

	//测试总体转置矩阵的临时测试对象
	StaticMatrix testForTemp_A_nxn;
	StaticMatrix testForTemp_B_nxn;
	StaticMatrix testTemp_nxn;
	MatrixMultiplier m_testMulti;

};



#endif /* EIGEN_BASIC_GENERALIZEDEIGENSOLVERFORREAL_H_ */
