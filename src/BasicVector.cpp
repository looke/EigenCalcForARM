/*
 * BasicVector.cpp
 *
 *  Created on: 2017年4月22日
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

//取向量的元素
double BasicVector::getElement(int index)
{};

//设置向量的元素
bool BasicVector::setElement(int index, double value)
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

bool BasicVector::resetDimension(int input_newDimension)
{
	if(input_newDimension < 0 || input_newDimension > this->space)
	{
		return false;
	}
	this->dimension = input_newDimension;
	return true;
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
//对向量做归一化处理
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

//对向量所有元素置0
void BasicVector::resetVectorElementToZero()
{
	for(int i=0; i<this->dimension;i++)
	{
		this->setElement(i,0);
	}
};

/*
//异常使用的错误信息
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
