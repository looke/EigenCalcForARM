/*
 * StaticVector.cpp
 *
 *  Created on: 2017年4月25日
 *      Author: looke
 */

#include "..\include\vector\static\StaticVector.h"

StaticVector::StaticVector()
{
	this->init(0);
};

StaticVector::StaticVector(int input_dimension) throw(out_of_range)
{
	if(input_dimension < 0 || input_dimension > this->space)
	{
		throw out_of_range(getOutOfRangeErrorMessage(input_dimension));
	}
	this->init(input_dimension);
};


void StaticVector::init(int input_dimension)
{
	this->space = 20;
	for(int i=0; i<space; i++)
	{
		vector_static[i] = 1;
	}

	this->dimension = input_dimension;
	//this->p_Vector = vector_static;
};

//取向量的元素
double StaticVector::getElement(int index) throw(out_of_range)
{
	if(index < 0 || index >= this->dimension)
	{
		throw out_of_range(getOutOfRangeErrorMessage(index));
	}
	return this->vector_static[index];
};

//设置向量的元素
void StaticVector::setElement(int index, double value) throw(out_of_range)
{
	if(index < 0 || index >= this->dimension)
	{
		throw out_of_range(getOutOfRangeErrorMessage(index));
	}
	this->vector_static[index] = value;
};
