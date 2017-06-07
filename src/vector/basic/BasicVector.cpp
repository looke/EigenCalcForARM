/*
 * BasicVector.cpp
 *
 *  Created on: 2017��4��22��
 *      Author: looke
 */

#include "..\include\vector\basic\BasicVector.h"
#include "math.h"
#include "iostream"
using namespace std;

BasicVector::BasicVector()
{
	this->dimension = 0;
	this->init();
};

BasicVector::BasicVector(int input_dimension)
{
	if(input_dimension > 0)
	{
		this->dimension = input_dimension;
	}
	this->init();
};

//ȡ������Ԫ��
double BasicVector::getElement(int index)
{};

//����������Ԫ��
void BasicVector::setElement(int index, double value)
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

void BasicVector::resetDimension(int input_newDimension)
{
	this->dimension = input_newDimension;
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


void BasicVector::init()
{};
