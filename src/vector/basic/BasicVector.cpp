/*
 * BasicVector.cpp
 *
 *  Created on: 2017年4月22日
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

//取向量的元素
double BasicVector::getElement(int index) throw(length_error)
{};

//设置向量的元素
void BasicVector::setElement(int index, double value) throw(length_error)
{};

//求向量的模
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

//求向量的模^2
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

//向量维度
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

//拷贝向量的元素
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

//打印向量
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
