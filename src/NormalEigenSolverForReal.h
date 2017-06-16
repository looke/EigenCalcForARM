/*
 * NormalEigenSolverForReal.h
 *
 *  Created on: 2017年5月20日
 *      Author: looke
 */

#ifndef EIGEN_BASIC_NORMALEIGENSOLVERFORREAL_H_
#define EIGEN_BASIC_NORMALEIGENSOLVERFORREAL_H_

#include "BasicMatrix.h"
#include "DoubleShiftQRIteration.h"
#include "SingleShiftQRIteration.h"
#include "HessenbergDeflation.h"

class NormalEigenSolverForReal
{
public:
	//NormalEigenSolverForReal();
	NormalEigenSolverForReal(BasicMatrix* p_input_OpMatrix,BasicVector* p_input_Vector,BasicMatrix* p_input_QT_Total,BasicMatrix* p_input_Q_Total,BasicMatrix* p_input_QQTMatrix_It,BasicMatrix* p_input_OpMatrix_deflated,BasicMatrix* p_TempMatrix_Trans,BasicMatrix* p_TempMatrix);

	void init(BasicMatrix* p_input_OpMatrix,BasicVector* p_input_Vector,BasicMatrix* p_input_QT_Total,BasicMatrix* p_input_Q_Total,BasicMatrix* p_input_QQTMatrix_It,BasicMatrix* p_input_OpMatrix_deflated,BasicMatrix* p_TempMatrix_Trans,BasicMatrix* p_TempMatrix);
	void reload(BasicMatrix* p_input_OpMatrix,BasicVector* p_input_Vector,BasicMatrix* p_input_QT_Total,BasicMatrix* p_input_Q_Total,BasicMatrix* p_input_QQTMatrix_It,BasicMatrix* p_input_OpMatrix_deflated,BasicMatrix* p_TempMatrix_Trans,BasicMatrix* p_TempMatrix);
	//virtual ~NormalEigenSolverForReal(){};

	//查找并更新降阶点
	//bool findNewDeflationPoint();

	//根据降阶点 生成降阶的Hessenberg矩阵
	void generateDeflatedHessenbergMatrix();

	//将已降阶的变换矩阵 升级成为全尺寸变换矩阵
	void upgradeDeflatedQQTMatrix();

	//将单迭代变换矩阵 合并成为总体变换矩阵
	void updateQTMatrixTotal();
	void updateQMatrixTotal();

	//生成原始hessenberg操作矩阵
	void generateHessenbergOpMatrix();

	//更新原始hessenberg操作矩阵
	void updateHessenbergOpMatrix_By_QT();
	void updateHessenbergOpMatrix_By_Q();
	//计算特征值
	void calcEigenValue();

	//初始化特征值计算相关矩阵
	void initEigenCalcMatrix();

	//解算结束判断逻辑
	bool hasFinishedIteration();

	//2x2对角块是否为复数特征值判断
	bool isDiagonalBlockComplexEigen(BasicMatrix* p_Input_OpMatrix);

	//初步化为对角块以后，最后一步迭代，将对角线上的2x2对角块进行上三角化
	void lastStepIteration(int startIndex);

	//将矩阵重新缩小
	void resizeMatrixForDeflation();
	//将已降阶的中间矩阵 升级成为全尺寸变换矩阵
	void upgradeDeflatedTempMatrix();
	//BasicMatrix* getEigenValueMatrix();
	BasicMatrix* getOpHessenbergMatrix();

	//BasicMatrix* getQTMatrix_Iteration();
	//BasicMatrix* getQMatrix_Iteration();

	BasicMatrix* getQTMatrix_Total();
	BasicMatrix* getQMatrix_Total();

	BasicMatrix* getOpHessenbergMatrix_deflated();

	//BasicMatrix* getQTMatrix_Deflated_Iteration();
	//BasicMatrix* getQMatrix_Deflated_Iteration();

	//测试打印，QT_Total * OP * Q
	//void showQTOPxQ();
protected:

	//降阶起点索引指示，初始化为0
	int deflationStart;
	//降阶终点索引指示，初始化为n-1
	int deflationEnd;

	//原始操作矩阵
	BasicMatrix* p_OpMatrix;

	BasicVector* p_TransVector;

	//原始矩阵的总体变换矩阵 QT为左乘矩阵 Q为右乘矩阵
	BasicMatrix* p_QTMatrix_Total;
	BasicMatrix* p_QMatrix_Total;

	//原始矩阵的单迭代变换矩阵 QT为左乘矩阵 Q为右乘矩阵
	BasicMatrix* p_QQTMatrix_Iteration;

	//已降阶的 操作Hessenberg矩阵
	BasicMatrix* p_OpHessenbergMatrix_deflated;

	//中间过程矩阵
	BasicMatrix* p_TempMatrix_Trans;
	//中间过程矩阵
	BasicMatrix* p_TempMatrix;


	//原始操作矩阵对应的特征值矩阵
	//BasicMatrix* p_EigenValueMatrix;

	//原始操作Hessenberg矩阵
	//BasicMatrix* p_OpHessenbergMatrix;


	//原始矩阵的单迭代变换矩阵 QT为左乘矩阵 Q为右乘矩阵
	//BasicMatrix* p_QTMatrix_Iteration;
	//BasicMatrix* p_QMatrix_Iteration;





	//已降阶矩阵的单迭代变换矩阵Q\QT
	//BasicMatrix* p_QTMatrix_Deflated_Iteration;
	//BasicMatrix* p_QMatrix_Deflated_Iteration;

	//最后一步针对对角线2x2矩阵块的操作
	//BasicMatrix* p_LastStepMatrix_2x2;


	//双重步QR迭代器
	DoubleShiftQRIteration m_DoubleShifeQR;

	//单步QR迭代器
	SingleShiftQRIteration m_SingleShifeQR;

	//hessen格式化
	HessenbergFormular m_HessenbergForm;

	//乘法器
	MatrixMultiplier m_Multiplier;

	//转置器
	MatrixTransposer m_Transposer;


	//hessenberg降阶点查找器
	HessenbergDeflation m_HessenbergDeflation;

	//测试总体转置矩阵的临时测试对象
	//BasicMatrix* p_testForTemp_nxn;
	//MatrixMultiplier* p_testMulti;
};



#endif /* EIGEN_BASIC_NORMALEIGENSOLVERFORREAL_H_ */
