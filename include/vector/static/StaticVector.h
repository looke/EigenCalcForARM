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
	StaticVector(int input_dimension) throw(out_of_range);

	//取向量的元素
	double getElement(int index) throw(out_of_range);

	//设置向量的元素
	void setElement(int index, double value) throw(out_of_range);

protected:
	void init(int input_dimension);

	//int vector_length;
	double vector_static[20];

private:

};


#endif /* VECTOR_STATIC_STATICVECTOR_H_ */
