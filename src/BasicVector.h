/*
 * BasicVector.h
 *
 *  Created on: 2017��4��22��
 *      Author: looke
 */

#ifndef VECTOR_BASIC_BASICVECTOR_H_
#define VECTOR_BASIC_BASICVECTOR_H_

//#include <stdexcept>
//#include <string>
//#include <sstream>
//using namespace std;

class BasicVector
{

public:
	BasicVector();
	BasicVector(int input_dimension);// throw(out_of_range);

	//��������ģ
	double getNormOfVector();
	//��������ģ^2
	double getNormPowerOfVector();

	//����ά��
	int getDimension();

	bool resetDimension(int input_newDimension);

	//ȡ������Ԫ��
	virtual double getElement(int index);

	//����������Ԫ��
	virtual bool setElement(int index, double value);

	//����������һ������
	void normalizationVector();

	//����������Ԫ����0
	void resetVectorElementToZero();

	//����쳣��Ϣ
	//string getOutOfRangeErrorMessage(int lengthOrIndex);
	//����������Ԫ��
	void copyVector(BasicVector* input_Vector);

	void printVector();

	virtual ~BasicVector(){};

private:



protected:
	int dimension;
	int space;
	//double* p_Vector;
	virtual void init(int input_dimension);

};



#endif /* VECTOR_BASIC_BASICVECTOR_H_ */
