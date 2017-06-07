/*
 * StaticVector.h
 *
 *  Created on: 2017年4月25日
 *      Author: looke
 */

#ifndef VECTOR_STATIC_STATICVECTOR_H_
#define VECTOR_STATIC_STATICVECTOR_H_
#include "..\include\vector\basic\BasicVector.h"

class StaticVector : public BasicVector
{
public:
	StaticVector();
	StaticVector(int input_dimension);

	//取向量的元素
	double getElement(int index);

	//设置向量的元素
	void setElement(int index, double value);

protected:
	void init(int input_dimension);

	//int vector_length;
	double vector_static[20];

private:

};


#endif /* VECTOR_STATIC_STATICVECTOR_H_ */
