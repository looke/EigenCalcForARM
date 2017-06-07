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
	this->space = 20;
	this->dimension = 0;
	this->init();
};

BasicVector::BasicVector(int input_dimension) throw(length_error)
{
	this->space = 20;
	if(input_dimension < 0 || input_dimension > this->space)
	{
		throw length_error("BasicVector Creation Exception: length out of range.");
	}
	if(input_dimension > 0)
	{
		this->dimension = input_dimension;
	}
	this->init();
};

//ȡ������Ԫ��
double BasicVector::getElement(int index) throw(length_error)
{};

//����������Ԫ��
void BasicVector::setElement(int index, double value) throw(length_error)
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

void BasicVector::resetDimension(int input_newDimension) throw(length_error)
{
	if(input_newDimension < 0 || input_newDimension > this->space)
	{
		throw length_error("BasicVector ResetDimension Exception: length out for range.");
	}
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
