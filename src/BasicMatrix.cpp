/*
 * BasicMatrix.cpp
 *
 *  Created on: 2017��4��1��
 *      Author: looke
 */


#include "BasicMatrix.h"
#include <iostream>

using namespace std;


BasicMatrix::BasicMatrix()
{
	//this->precision = 0.0000000001;
};


BasicMatrix::BasicMatrix(int inputRowNum, int inputColumnNum)
{
	//this->precision = 0.0000000001;
	if(inputRowNum > 0 && inputColumnNum >0)
	{
		this->rowNum = inputRowNum;
		this->columnNum = inputColumnNum;
		//call the default constructor
		//this->initMatrix();
	}
	else
	{
		cout<<"Illegal in put row/column numbers!" << endl;
		cout<<"Using default row/column number construct Matrix!" << endl;
		//BasicMatrix();
	}
};

void BasicMatrix::initMatrix()
{};

void BasicMatrix::printMatrix()
{

};


void BasicMatrix::setMatrixElement(int rNum, int cNum, double val)
{

};

double BasicMatrix::getMatrixElement(int rNum, int cNum)
{
};

//��ȡ����ָ��Ԫ�ص�ֵ(��Ԫ��ֵ�������Σ�С�ھ��ȵ�ֱֵ�ӷ���0)
double BasicMatrix::getMatrixElementRegulared(int rNum, int cNum, double lowEdge)
{};

/*
 * ����ָ������
 */
void BasicMatrix::swapRow(int from, int to)
{

};

/*
 *����ָ������
 */
void BasicMatrix::swapColumn(int from, int to)
{

};

//�����Խ�����Ԫ
void BasicMatrix::swapDiagElement(int from, int to)
{

};


/*
 * ����������Ϊ��λ��
 */
void BasicMatrix::resetMatrixToI()
{

};

/*
 * ����������Ϊ0��
 */
void BasicMatrix::resetMatrixToZero()
{

};

//��ȡָ��������
void BasicMatrix::getColumnVector(int columnNo, BasicVector* p_Vector)
{};

//��ȡָ��������
void BasicMatrix::getRowVector(int rowNo, BasicVector* p_Vector)
{};

//��ȡָ���Խ��Ӿ���������
void BasicMatrix::getSubMatrixColumnVector(int subMatrixIndex, int columnNo, BasicVector* p_Vector)
{};

//��ȡָ���Խ��Ӿ���������
void BasicMatrix::getSubMatrixRowVector(int subMatrixIndex, int rowNo, BasicVector* p_Vector)
{};
//��ȡָ���Խ��Ӿ���hessenberg������
void BasicMatrix::getSubMatrixHessenColumnVector(int subMatrixIndex, BasicVector* p_Vector)
{};
/*
 * ��������ӦԪ��������(���������)
 */
double BasicMatrix::calcMaxDifferentialWithCheck(BasicMatrix* targetMatrix)
{
	if(targetMatrix->rowNum == this->rowNum && targetMatrix->columnNum == this->columnNum)
	{
		double aMatrixElement = 0;
		double bMatrixElement = 0;
		double tempDiff = 0;
		double maxDiff = 0;

		for(int i=0; i<this->rowNum; i++)
		{
			for(int j=0; j<this->columnNum; j++)
			{
				//�������������ӦԪ�صĲ�ֵ�ľ���ֵ
				aMatrixElement = this->getMatrixElement(i,j);
				bMatrixElement = targetMatrix->getMatrixElement(i,j);
				tempDiff = aMatrixElement-bMatrixElement;
				if(tempDiff < 0)
				{
					tempDiff = 0 - tempDiff;
				}

				//�Ա�����ֵ�����ݽ����������ֵ
				if(maxDiff < tempDiff)
				{
					maxDiff = tempDiff;
				}
			}
		}
		return maxDiff;
	}
	else
	{
		return -1;
	}
};

/*
 * ��������ӦԪ��������(�����������)
 */
double BasicMatrix::calcMaxDifferentialNoCheck(BasicMatrix* targetMatrix)
{
	double aMatrixElement = 0;
	double bMatrixElement = 0;
	double tempDiff = 0;
	double maxDiff = 0;

	for(int i=0; i<this->rowNum; i++)
	{
		for(int j=0; j<this->columnNum; j++)
		{
			//�������������ӦԪ�صĲ�ֵ�ľ���ֵ
			aMatrixElement = this->getMatrixElement(i,j);
			bMatrixElement = targetMatrix->getMatrixElement(i,j);
			tempDiff = aMatrixElement-bMatrixElement;
			if(tempDiff < 0)
			{
				tempDiff = 0 - tempDiff;
			}

			//�Ա�����ֵ�����ݽ����������ֵ
			if(maxDiff < tempDiff)
			{
				maxDiff = tempDiff;
			}
		}
	}
	return maxDiff;
};

/*
 * ��������Ԫ�� ���������
 *
 * ���������������һ�£���ִ�п���������1
 *
 * �����������������һ�£���ֱ�ӷ���0
 */
bool BasicMatrix::copyMatrixElementWithCheck(BasicMatrix* input_Matrix)
{
	if(this->rowNum == input_Matrix->rowNum &&
		this->columnNum== input_Matrix->columnNum)
	{
		for(int i=0; i<this->rowNum;i++)
		{
			for(int j=0; j<this->columnNum; j++)
			{
				this->setMatrixElement(i,j, input_Matrix->getMatrixElement(i,j));
			}
		}
		return 1;
	}
	return 0;
};

/*
 * ��������Ԫ�� �����������
 *
 * ���������������һ�£���ִ�п���
 *
 * �����������������һ�£��������쳣
 */
void BasicMatrix::copyMatrixElementNoCheck(BasicMatrix* input_Matrix)
{
	for(int i=0; i<this->rowNum;i++)
	{
		for(int j=0; j<this->columnNum; j++)
		{
			this->setMatrixElement(i,j, input_Matrix->getMatrixElement(i,j));
		}
	}
};

/*
 * �ж��Ƿ�Ϊ�ԳƷ���
 */
bool BasicMatrix::isInputSymmetryMatrix()
{
	if(this->rowNum == this->columnNum)
	{
		for(int i=0; i<this->rowNum; i++)
		{
			for(int j=i; j<this->columnNum; j++)
			{
				if(i==j)
				{
					continue;
				}
				if(this->getMatrixElement(i,j) != this->getMatrixElement(j,i))
				{
					return false;
				}
			}
		}
		return true;
	}
	else
	{
		return false;
	}
};

//��������Ԫ�س������
void BasicMatrix::rMatrix(double r)
{
	for(int i=0; i<this->rowNum;i++)
	{
		for(int j=0; j<this->columnNum; j++)
		{
			this->setMatrixElement(i,j, this->getMatrixElement(i,j)*r);
		}
	}
};

//����ָ���е�����Ԫ�ط��ŷ�ת
void BasicMatrix::reverseSignOfColumn(int columnIndex)
{
	for(int i=0; i<this->rowNum;i++)
	{
		this->setMatrixElement(i,columnIndex, 0 - this->getMatrixElement(i,columnIndex));
	}
};

//�Ծ�������Ԫ�ؽ�������,������ֵС�ھ��ȵ�Ԫ������Ϊ0
void BasicMatrix::regularZeroElement()
{
	double normF = this->FrobeniousNorm();
	double lowEdge = normF * this->precision;
	double temp;
	for(int i=0; i<this->rowNum;i++)
	{
		for(int j=0; j<this->columnNum;j++)
		{
			temp = fabs(this->getMatrixElement(i,j));
			if(temp < lowEdge)
			{
				this->setMatrixElement(i,j,0.0);
			}
		}
	}
};

//Ѱ�����Խ�����Ϊ0��Ԫ��,��������λ��
//����ж��0Ԫ,��������ֵ�����Ǹ�(������½ǵ�)
int BasicMatrix::indexOfZeroOnDiagonal()
{
	double temp;
	for(int i=this->rowNum-1; i>=0;i--)
	{
		temp = fabs(this->getMatrixElement(i,i));
		if(temp < this->precision)
		{
			return i;
		}
	}
	return -1;
};

/*
 * �Խ���ȫ��Ԫ�ؼ�ȥָ��ֵ
 */
void BasicMatrix::diagonalSubtraction(double subValue)
{
	double temp;
	for(int i=0; i<this->rowNum; i++)
	{
		temp = this->getMatrixElement(i,i);
		temp = temp - subValue;
		this->setMatrixElement(i,i,temp);
	}
};

/*
 * �Խ���ȫ��Ԫ�ؼ���ָ��ֵ
 */
void BasicMatrix::diagonalAddition(double addValue)
{
	double temp;
	for(int i=0; i<this->rowNum; i++)
	{
		temp = this->getMatrixElement(i,i);
		temp = temp + addValue;
		this->setMatrixElement(i,i,temp);
	}
};

/*
 * ��ʼ��������Ϣ
 */
void BasicMatrix::initPrecision()
{
	int size = sizeof(double);

	switch(size)
	{
		case 4: //4�ֽ�32λ 23λ����λ log2^23 С�����7λ����
		{
			this->precision = 0.0000001;
			break;
		}
		case 8: //8�ֽ�64λ 52λ����λ log2^52 С�����15λ����
		{
			this->precision = 0.000000000000001;
			break;
		}
		default:
			this->precision = 0.0000001;
	}
};

//get������Ϣ
double BasicMatrix::getPrecision()
{
	return this->precision;
};

/*
 * ��������Frobenious���� ||A||f ����Ԫ�ؾ���ֵ��ƽ�����ٿ�ƽ��
 */
double BasicMatrix::FrobeniousNorm()
{
	double normValue = 0.0;
	double temp;
	for(int i=0; i<this->rowNum; i++)
	{
		for(int j=0; j<this->columnNum; j++)
		{
			temp = fabs(this->getMatrixElement(i,j));
			normValue = normValue + temp*temp;
		}
	}

	normValue = sqrt(normValue);
	return normValue;
};

/*
 * �����趨����ά��
 */
void BasicMatrix::resizeMatrix(int row, int column)
{
	if(row < this->space && column < this->space)
	{
		this->rowNum = row;
		this->columnNum = column;
		this->resetMatrixToI();
	}
};

/*
 * //�жϾ����Ƿ�����upper hessenberg���� ���ڼ�Сֵ��0����
 */
bool BasicMatrix::isUpperHessenbergMatrix()
{
	if(this->rowNum != this->columnNum)
	{
		return false;
	}

	double normF = this->FrobeniousNorm();
	double lowEdge = normF * this->precision;

	double temp;
	for(int i=rowNum-1; i>=0; i--)
	{
		for(int j=i-2; j>=0; j--)
		{
			temp = this->getMatrixElementRegulared(i,j,lowEdge);
			if(0 != temp)
			{
				return false;
			}
		}
	}
	return true;

};

/*
 * �жϾ����Ƿ�����upper Triangle���� ���ڼ�Сֵ��0����
 */
bool BasicMatrix::isUpperTriangleMatrix()
{
	if(this->rowNum != this->columnNum)
	{
		return false;
	}
	double normF = this->FrobeniousNorm();
	double lowEdge = normF * this->precision;

	double temp;
	for(int i=rowNum-1; i>=0; i--)
	{
		for(int j=i-1; j>=0; j--)
		{
			temp = this->getMatrixElementRegulared(i,j,lowEdge);
			if(0 != temp)
			{
				return false;
			}
		}
	}
	return true;
};
