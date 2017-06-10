/*
 * StaticVector.h
 *
 *  Created on: 2017��4��25��
 *      Author: looke
 */

#ifndef VECTOR_STATIC_STATICVECTOR_H_
#define VECTOR_STATIC_STATICVECTOR_H_
#include "BasicVector.h"

class StaticVector : public BasicVector
{
public:
	StaticVector();
	StaticVector(int input_dimension);

	//ȡ������Ԫ��
	double getElement(int index);

	//����������Ԫ��
	bool setElement(int index, double value);

protected:
	void init(int input_dimension);

	//int vector_length;
	double vector_static[20];

private:

};


#endif /* VECTOR_STATIC_STATICVECTOR_H_ */
