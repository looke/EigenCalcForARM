/*
 * BasicVector.h
 *
 *  Created on: 2017��4��22��
 *      Author: looke
 */

#ifndef VECTOR_BASIC_BASICVECTOR_H_
#define VECTOR_BASIC_BASICVECTOR_H_

#include <stdexcept>
#include <string>
#include <sstream>
using namespace std;

class BasicVector
{

public:
	BasicVector();
	BasicVector(int input_dimension) throw(out_of_range);

	//��������ģ
	double getNormOfVector();
	//��������ģ^2
	double getNormPowerOfVector();

	//����ά��
	int getDimension();

	void resetDimension(int input_newDimension) throw(out_of_range);

	//ȡ������Ԫ��
	virtual double getElement(int index) throw(out_of_range);

	//����������Ԫ��
	virtual void setElement(int index, double value) throw(out_of_range);

	//����쳣��Ϣ
	string getOutOfRangeErrorMessage(int lengthOrIndex);
	//����������Ԫ��
	void copyVector(BasicVector* input_Vector);

	void printVector();

	virtual ~BasicVector(){};

private:



protected:
	int dimension;
	int space;
	//double* p_Vector;
	virtual void init();

};



#endif /* VECTOR_BASIC_BASICVECTOR_H_ */
