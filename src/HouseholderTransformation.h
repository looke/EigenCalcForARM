/*
 * HouseholderTransformation.h
 *
 *  Created on: 2017��4��24��
 *      Author: looke
 */

#ifndef TRANSFORMATION_BASIC_HOUSEHOLDERTRANSFORMATION_H_
#define TRANSFORMATION_BASIC_HOUSEHOLDERTRANSFORMATION_H_


#include "BasicMatrix.h"
#include "BasicVector.h"

class HouseholderTransformation
{
public:
	//HouseholderTransformation();
	HouseholderTransformation(BasicVector* p_input_Vector);

	//����Householder�任����,����������,ת��Ϊ��һ��Ȼ������e1
	//isReverseElement:ָ���Ƿ��Ԫ�ط��Ž��з�ת�Ա㽵�;�����ʧ
	bool getHouseholderMatrixToE1_ReverseElement(BasicMatrix* p_HouseholderMatrix);

	//����Householder�任����,����������,ת��Ϊ��n��Ȼ������en
	//isReverseElement:ָ���Ƿ��Ԫ�ط��Ž��з�ת�Ա㽵�;�����ʧ
	bool getHouseholderMatrixToEn_ReverseElement(BasicMatrix* p_HouseholderMatrix);


	void init(BasicVector* p_input_Vector);
	void reload(BasicVector* p_input_Vector);

	virtual ~HouseholderTransformation(){};

private:


protected:
	//BasicMatrix* p_HouseholderMatrix;
	BasicVector* p_HouseholderVector;

	//����Householder Vector ��������Householder���� W = I-U*UT
	void generateHouseholderMatrixByVector(BasicMatrix* p_HouseholderMatrix);

};


#endif /* TRANSFORMATION_BASIC_HOUSEHOLDERTRANSFORMATION_H_ */
