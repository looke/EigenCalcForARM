/*
 * BasicVector.h
 *
 *  Created on: 2017年4月22日
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

	//求向量的模
	double getNormOfVector();
	//求向量的模^2
	double getNormPowerOfVector();

	//向量维度
	int getDimension();

	bool resetDimension(int input_newDimension);

	//取向量的元素
	virtual double getElement(int index);

	//设置向量的元素
	virtual bool setElement(int index, double value);

	//输出异常信息
	//string getOutOfRangeErrorMessage(int lengthOrIndex);
	//拷贝向量的元素
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
