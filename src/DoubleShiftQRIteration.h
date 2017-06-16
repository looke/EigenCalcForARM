/*
 * DoubleShiftQRIteration.h
 *
 *  Created on: 2017年5月17日
 *      Author: looke
 */

#ifndef EIGEN_BASIC_DOUBLESHIFTQRITERATION_H_
#define EIGEN_BASIC_DOUBLESHIFTQRITERATION_H_

#include "BasicMatrix.h"
#include "MatrixTransposer.h"
#include "BasicVector.h"
#include "HessenbergFormular.h"
#include "GivensTransformation.h"
#include "HouseholderTransformation.h"

class DoubleShiftQRIteration
{
public:
	//DoubleShiftQRIteration();
	DoubleShiftQRIteration(BasicMatrix* p_input_OpMatrix,BasicVector* p_input_TransVector,BasicMatrix* p_input_QTMatrix_Total,BasicMatrix* p_input_QQTMatrix_Step,BasicMatrix* p_input_TempMatrix);

	void init(BasicMatrix* p_input_OpMatrix,BasicVector* p_input_TransVector,BasicMatrix* p_input_QTMatrix_Total,BasicMatrix* p_input_QQTMatrix_Step,BasicMatrix* p_input_TempMatrix);
	void reload(BasicMatrix* p_input_OpMatrix,BasicVector* p_input_TransVector,BasicMatrix* p_input_QTMatrix_Total,BasicMatrix* p_input_QQTMatrix_Step,BasicMatrix* p_input_TempMatrix);
	//virtual ~DoubleShiftQRIteration(){};

	//Wilkinson双位移QR迭代  显式
	//void wilkinson_EX_QRIteration();

	//Wilkinson双位移QR迭代 隐式
	void wilkinson_IM_QRIteration();
	//作为迭代接口单独使用
	void wilkinson_IM_QRIteration_Single();

	BasicMatrix* getOpHessenbergMatrix();

	//获取总体转换矩阵
	BasicMatrix* getQTMatrix_Total();
	//BasicMatrix* getQMatrix_Total();

	//Wilkinson双位移QR迭代 隐式 -单步
	void wilkinson_IM_QRIteration_Step();

	//Wilkinson位移QR迭代 隐式 初始化
	void initForWilkinsonImplicitQR();

	//Wilkinson双位移QR迭代 隐式 -收尾
	void endForWilkinsonImplicitQR();

	//生成hessenberg操作矩阵
	void generateHessenbergOpMatrix();

	//根据hessenberg矩阵最后2x2子矩阵 计算生成
	void generateWilkinsonShift();

	//B = (A-pI)*(A-tI); 根据hessenberg矩阵 生成双步移位后的矩阵第一列
	void generateBMatrixFirstColumn();

	//将Q子矩阵根据当前迭代进度,升级为全尺寸Q矩阵
	void upgradeQQTSubMatrix(int iterateNum);

	//将Qn-1子矩阵,升级为全尺寸Q矩阵
	void upgradeQQTLastSubMatrix();

	//隐式QR迭代 更新hessenberg操作矩阵
	//void updateHessenbergOpMatrix_IM_QRIteration();
	void updateHessenbergOpMatrix_By_QT_IM_QRIteration();
	void updateHessenbergOpMatrix_By_Q_IM_QRIteration();
	//隐式QR迭代 更新总体转换矩阵Q QT
	//void updateTotalQQT();
	void updateQT_Total();
	//void updateQ_Total();
protected:

	//p+t 两个wilkinson位移值的和
	double trace;

	//p*t 两个wilkinson位移值的乘积
	double determinant;

	// B = (A-pI)*(A-tI); 双步移位后的矩阵第一列
	double b11;
	double b21;
	double b31;

	//操作矩阵
	BasicMatrix* p_OpMatrix;

	//用于计算隐式双步QR跌代单步转换矩阵Qi/Qn-1 的向量---3维/2维
	BasicVector* p_TransVectorForQStep;

	//Q 矩阵 隐式迭代  总体Q用于右乘OP矩阵
	//BasicMatrix* p_QMatrix_Implicit_Total;
	//QT 矩阵 隐式迭代 总体QT用于左乘OP矩阵
	BasicMatrix* p_QTMatrix_Implicit_Total;

	//Q/QT 矩阵 隐式迭代 分步 QT用于右乘/左乘OP矩阵
	BasicMatrix* p_QQTMatrix_Implicit_Step;

	//中间过程矩阵
	BasicMatrix* p_TempMatrix;


	//Hessenberg矩阵
	//BasicMatrix* p_OpHessenbergMatrix;

	//用于计算隐式双步QR跌代单步转换矩阵Qi的向量---3维
	//BasicVector* p_TransVectorForQStep;
	//用于计算隐式双步QR跌代单步转换矩阵Qn-1的向量---2维
	//BasicVector* p_TransVectorForQ_LastStep;



	//Q 矩阵 隐式迭代 分步 Q用于右乘OP矩阵
	//BasicMatrix* p_QMatrix_Implicit_Step;
	//QT 矩阵 隐式迭代 分步 QT用于左乘OP矩阵
	//BasicMatrix* p_QTMatrix_Implicit_Step;


	//Qi 子矩阵 显式迭代 分步 Qi用于右乘OP矩阵 3x3
	//BasicMatrix* p_QSubMatrix_Implicit_Step;

	//Qn-1 子矩阵 显式迭代 分步 Qn-1用于右乘OP矩阵 2x2
	//BasicMatrix* p_QSubMatrix_Implicit_LastStep;

	//QTi 子矩阵 显式迭代 分步 QTi用于左乘OP矩阵 3x3
	//BasicMatrix* p_QTSubMatrix_Implicit_Step;

	//QTn-1 子矩阵 显式迭代 分步 QT用于左乘OP矩阵 2x2
	//BasicMatrix* p_QTSubMatrix_Implicit_LastStep;

	//QR分解
	//QRDecomposition* p_QRDecomp;

	//hessen格式化
	HessenbergFormular m_HessenbergForm;

	//Givens变换器
	GivensTransformation m_GivensTrans;

	//householder变换
	HouseholderTransformation m_HouseholderTrans;

	//乘法器
	MatrixMultiplier m_Multiplier;

	//转置器
	MatrixTransposer m_Transposer;

};



#endif /* EIGEN_BASIC_DOUBLESHIFTQRITERATION_H_ */
