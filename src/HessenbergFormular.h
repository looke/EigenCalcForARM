/*
 * HessenbergFormular.h
 *
 *  Created on: 2017年5月8日
 *      Author: looke
 */

#ifndef TRANSFORMATION_BASIC_HESSENBERGFORMULAR_H_
#define TRANSFORMATION_BASIC_HESSENBERGFORMULAR_H_

#include "BasicMatrix.h"
#include "BasicVector.h"
#include "HouseholderTransformation.h"
#include "MatrixMultiplier.h"
#include "MatrixTransposer.h"

class HessenbergFormular
{
public:
	//HessenbergFormular();
	HessenbergFormular(BasicMatrix* p_Input_OpMatrix, BasicMatrix* p_Input_preTransMatrix,BasicMatrix* p_Input_transMatrix, BasicMatrix* p_Input_tempMatrix);

	void init(BasicMatrix* p_Input_OpMatrix, BasicMatrix* p_Input_preTransMatrix,BasicMatrix* p_Input_transMatrix, BasicMatrix* p_Input_tempMatrix);
	void reload(BasicMatrix* p_Input_OpMatrix, BasicMatrix* p_Input_preTransMatrix,BasicMatrix* p_Input_transMatrix, BasicMatrix* p_Input_tempMatrix);
	//void resizeSubMatrix(int rowAndColumnNumber);

	//根据当前迭代  初始化子变换矩阵
	//void initSubHouseholderTrans(int iterateNum);

	//Hessenberg格式化
	void formularUpperHessnbergMatrix();

	//根据当前迭代 计算生成子变换矩阵
	void generateSubHouseholderTrans(int iterateNum);

	//将子变换矩阵升级为全维度变换阵
	void upgradeSubHouseholderTrans(int iterateNum);

	void updatePreTransMatrix();
	//void updateAfterTransMatrix();

	void updateOpMatrix();

	BasicMatrix* getOpMatrix();
	BasicMatrix* getPreTransMatrix();
	//BasicMatrix* getAfterTransMatrix();
	BasicMatrix* getTransMatrix();
	//BasicMatrix* getSubTransMatrix();

	//virtual ~HessenbergFormular(){};
protected:
	BasicMatrix* p_OpMatrix;

	//总体变换阵--左乘
	BasicMatrix* p_preTransMatrix;

	//右乘变换阵 为左乘变换阵的转置
	//BasicMatrix* p_afterTransMatrix;

	//Householder变换阵
	BasicMatrix* p_transMatrix;

	//Householder子变换阵
	//BasicMatrix* p_subTransMatrix;

	//中间过程矩阵
	BasicMatrix* p_tempMatrix;

	//Householder变换
	HouseholderTransformation m_HouseholderTrans;

	//乘法器
	MatrixMultiplier m_Multiplier;

	//转置器
	MatrixTransposer m_Transposer;
};



#endif /* TRANSFORMATION_BASIC_HESSENBERGFORMULAR_H_ */
