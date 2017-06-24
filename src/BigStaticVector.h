/*
 * BigStaticVector.h
 *
 *  Created on: 2017年6月24日
 *      Author: looke
 */

#ifndef BIGSTATICVECTOR_H_
#define BIGSTATICVECTOR_H_

#include "BasicVector.h"

class BigStaticVector : public BasicVector
{
public:
	BigStaticVector();
	BigStaticVector(int input_dimension);

	//取向量的元素
	double getElement(int index);

	//设置向量的元素
	bool setElement(int index, double value);

protected:
	void init(int input_dimension);

	//int vector_length;
	double vector_static[20];

private:

};



#endif /* BIGSTATICVECTOR_H_ */
