/*
 * BigStaticVector.h
 *
 *  Created on: 2017��6��24��
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



#endif /* BIGSTATICVECTOR_H_ */
