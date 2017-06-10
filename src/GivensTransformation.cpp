/*
 * GivensTransformation.cpp
 *
 *  Created on: 2017��4��23��
 *      Author: looke
 */

#include "GivensTransformation.h"
#include "math.h"

GivensTransformation::GivensTransformation()
{
	//this->isUsingPreElement = true;
};

GivensTransformation::GivensTransformation(BasicVector* p_input_Vector)
{
	this->init(p_input_Vector);
};

/*
 * ������������������Givens���������ҳ˸�����������ָ��λ�õ�Ԫ��Ϊ0
 * �ҳ�Givens�����Ӧ������Ϊ����������
 */
bool GivensTransformation::getGivensMatrixAfterMultiple(int elementIndexToZero, BasicMatrix* p_GivensMatrix)
{
	int size = this->p_GivensVector->getDimension();
	if(p_GivensMatrix->rowNum != size || p_GivensMatrix->columnNum != size)
	{
		return false;
	}

	//getGivensMatrixPreMultiple(elementIndexToZero);
	int i,j;

	//�ҳ�Givens���� ʹ�ú���Ԫ�ؽ�����Ԫ
	if(this->p_GivensVector->getDimension()-1 == elementIndexToZero)
	{
		j = elementIndexToZero;
		i = 0;
	}
	else if(this->p_GivensVector->getDimension()-1 != elementIndexToZero)
	{
		j = elementIndexToZero;
		i = j+1;
	}

	double i_Value = this->p_GivensVector->getElement(i);
	double j_Value = this->p_GivensVector->getElement(j);

	//����Givens����
	GenerateGivensMatrix(i, i_Value, j, j_Value,p_GivensMatrix);

	//��p_GivensMatrix����ת��
	double temp;
	for(int i=0; i<p_GivensMatrix->rowNum;i++)
	{
		for(int j=i+1; j<p_GivensMatrix->columnNum;j++)
		{
			temp = p_GivensMatrix->getMatrixElement(i,j);

			p_GivensMatrix->setMatrixElement(i,j,p_GivensMatrix->getMatrixElement(j,i));
			p_GivensMatrix->setMatrixElement(j,i,temp);
		}
	}

	return true;
};

/*
 * ������������������Givens����������˸�����������ָ��λ�õ�Ԫ��Ϊ0
 * ���Givens�����Ӧ������Ϊ����������
 */
bool GivensTransformation::getGivensMatrixPreMultiple(int elementIndexToZero, BasicMatrix* p_GivensMatrix)
{
	int size = this->p_GivensVector->getDimension();
	if(p_GivensMatrix->rowNum != size || p_GivensMatrix->columnNum != size)
	{
		return false;
	}

	int i,j;
	//���Givens���� ʹ��ǰ��Ԫ�ؽ�����Ԫ
	if(0 == elementIndexToZero)
	{
		j = 0;
		i = (this->p_GivensVector->getDimension())-1;
	}
	else if(0 != elementIndexToZero)
	{
		j = elementIndexToZero;
		i = j-1;
	}

	double i_Value = this->p_GivensVector->getElement(i);
	double j_Value = this->p_GivensVector->getElement(j);

	//����Givens����
	GenerateGivensMatrix(i, i_Value, j, j_Value,p_GivensMatrix);

	return true;
};

void GivensTransformation::GenerateGivensMatrix(int i, double i_Value, int j, double j_Value, BasicMatrix* p_GivensMatrix)
{
	double cosine,sine;
	double denominator = sqrt(i_Value*i_Value + j_Value*j_Value);

	if(0 != denominator)
	{
		cosine = i_Value / denominator;
		sine = j_Value / denominator;
	}
	else
	{
		cosine = 1;
		sine = 0;
	}

	p_GivensMatrix->resetMatrixToI();

	p_GivensMatrix->setMatrixElement(i,i,cosine);
	p_GivensMatrix->setMatrixElement(j,j,cosine);

	p_GivensMatrix->setMatrixElement(j,i,0-sine);
	p_GivensMatrix->setMatrixElement(i,j,sine);
};

/*
void GivensTransformation::setIsUsingPreElement(bool input)
{
	this->isUsingPreElement = input;
};
*/

void GivensTransformation::init(BasicVector* p_input_Vector)
{
	this->p_GivensVector = p_input_Vector;
	//this->p_GivensMatrix = p_GivensMatrix;
};

void GivensTransformation::reload(BasicVector* p_input_Vector)
{
	this->init(p_input_Vector);
};
