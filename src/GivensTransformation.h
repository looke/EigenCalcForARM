/*
 * GivensTransformation.h
 *
 *  Created on: 2017��4��22��
 *      Author: looke
 */

#ifndef TRANSFORMATION_GIVENSTRANSFORMATION_H_
#define TRANSFORMATION_GIVENSTRANSFORMATION_H_

#include "BasicMatrix.h"
#include "BasicVector.h"

class GivensTransformation
{
public:
	GivensTransformation();
	GivensTransformation(BasicVector* p_input_Vector);

	//����Givens�任����,������������,�Լ�Ҫ������Ԫ������.���ɾ��������������ָ��Ԫ����ԪΪ0
	//���Givens����ʹ��ǰ��Ԫ��
	bool getGivensMatrixPreMultiple(int elementIndexToZero, BasicMatrix* p_GivensMatrix);

	//����Givens�任����,������������,�Լ�Ҫ������Ԫ������.���ɾ������ҳ�������ָ��Ԫ����ԪΪ0
	//�ҳ�Givens����ʹ�ú���Ԫ��
	bool getGivensMatrixAfterMultiple(int elementIndexToZero, BasicMatrix* p_GivensMatrix);

	void init(BasicVector* p_input_Vector);
	void reload(BasicVector* p_input_Vector);

	virtual ~GivensTransformation(){};

	//void setIsUsingPreElement(bool input);

private:


protected:
	//�Ƿ�ʹ��ָ������Ԫ�ص�ǰһ��Ԫ�ع���Givens����,Ĭ�����ΪTrue;
	//�ڽ���Hessenberg-Tirangle��ʽ���Ĺ�����,��B�������Ԫ��Ҫ���˴�����ΪFalse;
	//bool isUsingPreElement;

	//BasicMatrix* p_GivensMatrix;
	BasicVector* p_GivensVector;
	void GenerateGivensMatrix(int i, double i_Value, int j, double j_Value, BasicMatrix* p_GivensMatrix);

	//string getMatrixSizeErrorMessage(int vectorSize, int rowNumber, int columnNumber);
};



#endif /* TRANSFORMATION_GIVENSTRANSFORMATION_H_ */
