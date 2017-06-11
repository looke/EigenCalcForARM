/*
 * HouseholderTransformation.h
 *
 *  Created on: 2017年4月24日
 *      Author: looke
 */

#ifndef TRANSFORMATION_BASIC_HOUSEHOLDERTRANSFORMATION_H_
#define TRANSFORMATION_BASIC_HOUSEHOLDERTRANSFORMATION_H_


#include "BasicMatrix.h"
#include "BasicVector.h"

class HouseholderTransformation
{
public:
	HouseholderTransformation();
	HouseholderTransformation(BasicVector* p_input_Vector);

	//生成Householder变换矩阵,将输入向量,转换为第一自然基向量e1
	//isReverseElement:指定是否对元素符号进行反转以便降低精度损失
	bool getHouseholderMatrixToE1_ReverseElement(BasicMatrix* p_HouseholderMatrix);

	//生成Householder变换矩阵,将输入向量,转换为第n自然基向量en
	//isReverseElement:指定是否对元素符号进行反转以便降低精度损失
	bool getHouseholderMatrixToEn_ReverseElement(BasicMatrix* p_HouseholderMatrix);


	void init(BasicVector* p_input_Vector);
	void reload(BasicVector* p_input_Vector);

	virtual ~HouseholderTransformation(){};

private:


protected:
	//BasicMatrix* p_HouseholderMatrix;
	BasicVector* p_HouseholderVector;

	//根据Householder Vector 计算生成Householder矩阵 W = I-U*UT
	void generateHouseholderMatrixByVector(BasicMatrix* p_HouseholderMatrix);

};


#endif /* TRANSFORMATION_BASIC_HOUSEHOLDERTRANSFORMATION_H_ */
