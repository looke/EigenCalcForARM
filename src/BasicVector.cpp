/*
 * BasicVector.cpp
 *
 *  Created on: 2017��4��22��
 *      Author: looke
 */

#include "BasicVector.h"
#include "math.h"
#include "iostream"
using namespace std;

BasicVector::BasicVector()
{
	this->space = 20;
	this->dimension = 0;
	this->init(0);
};

BasicVector::BasicVector(int input_dimension)
{
	this->init(input_dimension);
};

//ȡ������Ԫ��
double BasicVector::getElement(int index)
{};

//����������Ԫ��
bool BasicVector::setElement(int index, double value)
{};

//��������ģ
double BasicVector::getNormOfVector()
{
	double temp = 0.0;
	double result = 0.0;
	for(int i=0; i<this->dimension; i++)
	{
		temp = this->getElement(i);
		result += temp*temp;
	}

	return sqrt(result);
};

//��������ģ^2
double BasicVector::getNormPowerOfVector()
{
	double temp = 0.0;
	double result = 0.0;
	for(int i=0; i<this->dimension; i++)
	{
		temp = this->getElement(i);
		result += temp*temp;
	}

	return result;
};

//����ά��
int BasicVector::getDimension()
{
	return this->dimension;
};

bool BasicVector::resetDimension(int input_newDimension)
{
	if(input_newDimension < 0 || input_newDimension > this->space)
	{
		return false;
	}
	this->dimension = input_newDimension;
	return true;
};

//����������Ԫ��
void BasicVector::copyVector(BasicVector* input_Vector)
{
	if(this->dimension == input_Vector->getDimension())
	{
		for(int i=0; i<this->dimension; i++)
		{
			this->setElement(i,input_Vector->getElement(i));
		}
	}

};
//����������һ������
void BasicVector::normalizationVector()
{
	double norm = this->getNormOfVector();
	double temp;

	for(int i=0; i<this->dimension;i++)
	{
		temp = this->getElement(i);
		temp = temp/norm;
		this->setElement(i,temp);
	}
};

//��ӡ����
void BasicVector::printVector()
{
	cout << "Vecotr: ";
	for(int i=0; i<this->dimension; i++)
	{
		cout << this->getElement(i) << "; ";
	}
	cout << endl;
};

/*
//�쳣ʹ�õĴ�����Ϣ
string BasicVector::getOutOfRangeErrorMessage(int lengthOrIndex)
{
	stringstream stream;
	stream<<lengthOrIndex;
	string detail_input_dimension = stream.str();
	stream.str("");
	stream<<this->space;
	string detail_space = stream.str();
	string exInfo = "Out of range. input:" +detail_input_dimension + ". Range: 0 to " + detail_space;
	return exInfo;
};
*/
void BasicVector::init(int input_dimension)
{};
