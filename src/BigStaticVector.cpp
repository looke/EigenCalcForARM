/*
 * BigStaticVector.cpp
 *
 *  Created on: 2017年6月24日
 *      Author: looke
 */

#include "BigStaticVector.h"

BigStaticVector::BigStaticVector()
{
	this->space = 20;
	this->init(0);
};

BigStaticVector::BigStaticVector(int input_dimension)
{
	this->init(input_dimension);
};


void BigStaticVector::init(int input_dimension)
{
	this->space = 20;
	for(int i=0; i<space; i++)
	{
		vector_static[i] = 1;
	}

	if(input_dimension < 0 || input_dimension > this->space)
	{
		this->dimension = this->space;
	}
	else
	{
		this->dimension = input_dimension;
	}

};

//取向量的元素
double BigStaticVector::getElement(int index)
{
	//if(index < 0 || index >= this->dimension)
	//{
		//throw out_of_range(getOutOfRangeErrorMessage(index));
	//}
	return this->vector_static[index];
};

//设置向量的元素
bool BigStaticVector::setElement(int index, double value)
{
	if(index < 0 || index >= this->dimension)
	{
		return false;
	}
	this->vector_static[index] = value;
	return true;
};


