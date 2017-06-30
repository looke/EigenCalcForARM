/*
 * StaticVector.cpp
 *
 *  Created on: 2017年4月25日
 *      Author: looke
 */

#include "StaticVector.h"

StaticVector::StaticVector()
{
	this->space = 10;
	this->init(10);
};

StaticVector::StaticVector(int input_dimension)
{
	this->init(input_dimension);
};


void StaticVector::init(int input_dimension)
{
	this->space = 10;
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
double StaticVector::getElement(int index)
{
	//if(index < 0 || index >= this->dimension)
	//{
		//throw out_of_range(getOutOfRangeErrorMessage(index));
	//}
	return this->vector_static[index];
};

//设置向量的元素
bool StaticVector::setElement(int index, double value)
{
	if(index < 0 || index >= this->dimension)
	{
		return false;
	}
	this->vector_static[index] = value;
	return true;
};
