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

StaticVector::StaticVector(int input_dimension) throw(length_error)
{
	if(input_dimension < 0 || input_dimension > this->dimension)
	{
		throw length_error("StaticVector Creation Exception: input_dimension out of range.");
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
double StaticVector::getElement(int index) throw(length_error)
{
	if(index < 0 || index > this->dimension)
	{
		throw length_error("StaticVector GetElement Exception: index out of range.");
	}
	return this->vector_static[index];
};

//设置向量的元素
void StaticVector::setElement(int index, double value) throw(length_error)
{
	if(index < 0 || index > this->dimension)
	{
		throw length_error("StaticVector SetElement Exception: index out of range.");
	}
	this->vector_static[index] = value;
};
