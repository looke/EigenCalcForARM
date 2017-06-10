/*
 * StaticMatrix.cpp
 *
 *  Created on: 2017��4��25��
 *      Author: looke
 */

#include "StaticMatrix.h"
#include <iostream>
#include "math.h"

using namespace std;

StaticMatrix::StaticMatrix()
{
	//this->precision = 0.0000000001;
	this->space = 20;
	this->rowNum = 1;
	this->columnNum = 1;
};

StaticMatrix::StaticMatrix(int inputRowNum, int inputColumnNum)
{
	initPrecision();
	this->space = 20;
	this->rowNum = inputRowNum;
	this->columnNum = inputColumnNum;
	for(int i=0; i<this->space; i++)
	{
		for(int j=0; j<this->space; j++)
		{
			this->matrixNxN[i][j] = 0;
		}
	}

	this->initMatrix();
};


void StaticMatrix::initMatrix()
{
	double tempValue = 1;
	for(int i=0; i<this->rowNum; i++)
	{
		for(int j=0; j<this->columnNum; j++)
		{
			this->matrixNxN[i][j] = tempValue;
			tempValue++;
		}
	}

	//��ʼ��������������
	//this->columnVector = StaticVector(this->rowNum);

	//this->rowVector = StaticVector(this->columnNum);

};

void StaticMatrix::setMatrixElement(int rNum, int cNum, double val)
{
	this->matrixNxN[rNum][cNum] = val;
};


double StaticMatrix::getMatrixElement(int rNum, int cNum)
{
	return this->matrixNxN[rNum][cNum];
};

/*
 * ��ȡ����ָ��Ԫ�ص�ֵ(��Ԫ��ֵ�������Σ�С�ھ��ȵ�ֱֵ�ӷ���0)
 */
double StaticMatrix::getMatrixElementRegulared(int rNum, int cNum, double lowEdge)
{
	if(fabs(this->matrixNxN[rNum][cNum]) < lowEdge)
	{
		return 0;
	}
	return this->matrixNxN[rNum][cNum];
};

void StaticMatrix::printMatrix()
{
	for(int i=0;i<this->rowNum;i++)
	{
		for(int j=0;j<this->columnNum;j++)
		{
			cout.width(10);
			cout<<this->matrixNxN[i][j]<<"\t";
		}
		cout<<endl;
	}
};

//������
void StaticMatrix::swapRow(int from, int to)
{
	if(from < this->rowNum && to < this->rowNum)
	{
		double temp[this->columnNum];

		for(int i=0; i < this->columnNum; i++)
		{
			temp[i] = matrixNxN[from][i];
		}

		for(int i=0; i < this->columnNum; i++)
		{
			matrixNxN[from][i] = matrixNxN[to][i];
		}

		for(int i=0; i < this->columnNum; i++)
		{
			matrixNxN[to][i] = temp[i];
		}
	}
	else
	{
		cout<<"Illegal row number." << endl;
	}
};

//������
void StaticMatrix::swapColumn(int from, int to)
{
	if(from < this->columnNum && to < this->columnNum)
	{
		double temp[this->rowNum];

		for(int i=0; i < this->rowNum; i++)
		{
			temp[i] = matrixNxN[i][from];
		}

		for(int i=0; i < this->rowNum; i++)
		{
			matrixNxN[i][from] = matrixNxN[i][to];
		}

		for(int i=0; i < this->rowNum; i++)
		{
			matrixNxN[i][to] = temp[i];
		}
	}
	else
	{
		cout<<"Illegal column number." << endl;
	}

};

//�����Խ�����Ԫ
void StaticMatrix::swapDiagElement(int from, int to)
{
	if(from < this->columnNum && to < this->columnNum && from < this->rowNum && to < this->rowNum)
	{
		double temp;
		temp = matrixNxN[from][from];

		matrixNxN[from][from] = matrixNxN[to][to];
		matrixNxN[to][to] = temp;
	}
};

//����������Ϊ��λ����
void StaticMatrix::resetMatrixToI()
{
	for(int i=0; i<this->rowNum; i++)
	{
		for(int j=0; j<this->columnNum; j++)
		{
			if(i==j)
			{
				matrixNxN[i][j] = 1;
			}
			else
			{
				matrixNxN[i][j] = 0;
			}
		}
	}
};

//����������Ϊ��λ����
void StaticMatrix::resetMatrixToZero()
{
	for(int i=0; i<this->rowNum; i++)
	{
		for(int j=0; j<this->columnNum; j++)
		{

			matrixNxN[i][j] = 0;

		}
	}
};

//��ȡָ��������
void StaticMatrix::getColumnVector(int columnNo, BasicVector* p_Vector)
{
	int newColumnVectorDimension = this->rowNum ;
	p_Vector->resetDimension(newColumnVectorDimension);

	for(int i=0; i<this->rowNum; i++)
	{
		p_Vector->setElement(i,matrixNxN[i][columnNo]);

	}
};

//��ȡָ��������
void StaticMatrix::getRowVector(int rowNo, BasicVector* p_Vector)
{
	int newRowVectorDimension = this->columnNum ;
	p_Vector->resetDimension(newRowVectorDimension);

	for(int i=0; i<this->columnNum; i++)
	{
		p_Vector->setElement(i,matrixNxN[rowNo][i]);
	}
};


//��ȡ�Խ����Ӿ���ָ��������
void StaticMatrix::getSubMatrixColumnVector(int subMatrixIndex, int columnNo, BasicVector* p_Vector)
{
	int newColumnVectorDimension = this->rowNum - subMatrixIndex;
	p_Vector->resetDimension(newColumnVectorDimension);

	int index = 0;
	for(int i=subMatrixIndex; i<this->rowNum; i++)
	{
		p_Vector->setElement(index,matrixNxN[i][subMatrixIndex+columnNo]);
		index++;
	}

};

//��ȡ�Խ����Ӿ���ָ��������
void StaticMatrix::getSubMatrixRowVector(int subMatrixIndex, int rowNo, BasicVector* p_Vector)
{
	int newRowVectorDimension = this->columnNum - subMatrixIndex;
	p_Vector->resetDimension(newRowVectorDimension);
	int index = 0;
	for(int i=subMatrixIndex; i<this->columnNum; i++)
	{
		p_Vector->setElement(index,matrixNxN[subMatrixIndex+rowNo][i]);
		index++;
	}

};

//��ȡָ���Խ��Ӿ���hessenberg������
void StaticMatrix::getSubMatrixHessenColumnVector(int subMatrixIndex, BasicVector* p_Vector)
{
	int newColumnVectorDimension = this->rowNum - subMatrixIndex - 1;
	p_Vector->resetDimension(newColumnVectorDimension);

	int index = 0;
	for(int i=subMatrixIndex+1; i<this->rowNum; i++)
	{
		p_Vector->setElement(index,matrixNxN[i][subMatrixIndex]);
		index++;
	}

};
