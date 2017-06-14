/*
 * HessenbergTriangleFormular.h
 *
 *  Created on: 2017年4月28日
 *      Author: looke
 */

#ifndef TRANSFORMATION_BASIC_HESSENBERGTRIANGLEFORMULAR_H_
#define TRANSFORMATION_BASIC_HESSENBERGTRIANGLEFORMULAR_H_

#include "BasicMatrix.h"
#include "MatrixMultiplier.h"
#include "QRDecomposition.h"
#include "GivensTransformation.h"


class HessenbergTriangleFormular
{
public:
	//HessenbergTriangleFormular();
	HessenbergTriangleFormular(BasicMatrix* input_Matrix_A, BasicMatrix* input_Matrix_B,BasicMatrix* input_QMatrix_Total,BasicMatrix* input_ZMatrix_Total,BasicMatrix* input_QZMatrix_Step,BasicMatrix* input_TempMatrix_Trans,BasicMatrix* input_TempMatrix);

	void initABMatrix();
	void formularABMatrix();
	void formularColumnVector(int columnIndex);

	void updateOpMatrix_A_ByQ();
	void updateOpMatrix_B_ByQ();
	void updateOpMatrixByZ();

	void updateQMatrix_Total();
	void updateZMatrix_Total();

	BasicMatrix* getMatrixA();
	BasicMatrix* getMatrixB();
	BasicMatrix* getHessenbergMatrixA();
	BasicMatrix* getTriangleMatrixB();
	//BasicMatrix* getMatrixQ_Step();
	//BasicMatrix* getMatrixZ_Step();
	BasicMatrix* getMatrixQZ_Step();
	BasicMatrix* getMatrixQ_Total();
	BasicMatrix* getMatrixZ_Total();


	void init(BasicMatrix* input_Matrix_A, BasicMatrix* input_Matrix_B,BasicMatrix* input_QMatrix_Total,BasicMatrix* input_ZMatrix_Total,BasicMatrix* input_QZMatrix_Step,BasicMatrix* input_TempMatrix_Trans,BasicMatrix* input_TempMatrix);
	void reload(BasicMatrix* input_Matrix_A, BasicMatrix* input_Matrix_B,BasicMatrix* input_QMatrix_Total,BasicMatrix* input_ZMatrix_Total,BasicMatrix* input_QZMatrix_Step,BasicMatrix* input_TempMatrix_Trans,BasicMatrix* input_TempMatrix);

	//virtual ~HessenbergTriangleFormular(){};
protected:

	//原始操作方阵A
	BasicMatrix* p_OpMatrix_A;
	//原始操作方阵B
	BasicMatrix* p_OpMatrix_B;

	//Hessenberg方阵A
	//BasicMatrix* p_OpHessenbergMatrix_A;

	//上三角方阵B
	//BasicMatrix* p_OpTriangleMatrix_B;

	//总体Q变换矩阵,用于汇总左乘操作矩阵
	BasicMatrix* p_QMatrix_Total;

	//总体Z变换矩阵,用于汇总右乘操作矩阵
	BasicMatrix* p_ZMatrix_Total;

	//单步Q变换矩阵,用于左乘操作矩阵
	//BasicMatrix* p_QMatrix_Step;

	//单步Z变换矩阵,用于右乘操作矩阵
	//BasicMatrix* p_ZMatrix_Step;

	//Q* A Z--->Hessenberg矩阵
	//Q* B Z--->上三角矩阵

	//单步QZ变换矩阵,用于左乘操作矩阵
	BasicMatrix* p_QZMatrix_Step;

	//中间转换矩阵
	BasicMatrix* p_TempMatrix_Trans;

	//中间过程矩阵
	BasicMatrix* p_TempMatrix;

	//QR分解
	QRDecomposition m_QRDecomp;

	//Givens变换
	GivensTransformation m_GivensTrans;

	//乘法器
	MatrixMultiplier m_Multiplier;

	//转置器
	MatrixTransposer m_Transposer;

private:

};

#endif /* TRANSFORMATION_BASIC_HESSENBERGTRIANGLEFORMULAR_H_ */
