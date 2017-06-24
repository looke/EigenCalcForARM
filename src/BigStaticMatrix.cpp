/*
 * BigStaticMatrix.cpp
 *
 *  Created on: 2017年6月24日
 *      Author: looke
 */

#include "BigStaticMatrix.h"
#include <iostream>
#include "math.h"

using namespace std;

BigStaticMatrix::BigStaticMatrix()
{
	initPrecision();
	this->space = 20;

	this->rowNum = 0;
	this->columnNum = 0;

	this->initMatrix();
	//初始化行列向量对象
	this->matrixVector = BigStaticVector(this->rowNum);
};

BigStaticMatrix::BigStaticMatrix(int inputRowNum, int inputColumnNum)
{
	initPrecision();
	this->space = 20;

	if(inputRowNum >= 0 && inputColumnNum >=0 && inputRowNum<=space && inputColumnNum<=space)
	{
		this->rowNum = inputRowNum;
		this->columnNum = inputColumnNum;
	}
	else
	{
		this->rowNum = 0;
		this->columnNum = 0;
	}

	this->initMatrix();
	//初始化行列向量对象
	this->matrixVector = BigStaticVector(this->rowNum);
};

/*
 * 按1,2,3,4...递增顺序初始化矩阵各元素
 */
void BigStaticMatrix::initMatrix()
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
};

bool BigStaticMatrix::setMatrixElement(int rNum, int cNum, double val)
{
	if(rNum >= 0 && rNum < rowNum && cNum >= 0 && cNum < columnNum)
	{
		this->matrixNxN[rNum][cNum] = val;
		return true;
	}
	return false;
};


double BigStaticMatrix::getMatrixElement(int rNum, int cNum)
{
	return this->matrixNxN[rNum][cNum];
};

/*
 * 获取矩阵指定元素的值(对元素值进行整形，小于精度的值直接返回0)
 */
double BigStaticMatrix::getMatrixElementRegulared(int rNum, int cNum, double lowEdge)
{
	if(fabs(this->matrixNxN[rNum][cNum]) < lowEdge)
	{
		return 0;
	}
	return this->matrixNxN[rNum][cNum];
};

void BigStaticMatrix::printMatrix()
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

//交换行
bool BigStaticMatrix::swapRow(int from, int to)
{
	if(from >=0 && from < this->rowNum && to >=0 && to < this->rowNum)
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
		return true;
	}
	else
	{
		return false;
	}
};

//交换列
bool BigStaticMatrix::swapColumn(int from, int to)
{
	if(from >=0 && from < this->columnNum && to >=0 && to < this->columnNum)
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
		return true;
	}
	else
	{
		return false;
	}

};

/*
 * 交换对角线主元
 * 仅方阵专用，非方阵返回false，并不会交换
 */
bool BigStaticMatrix::swapDiagElement(int from, int to)
{
	//非方阵
	if(this->rowNum != this->columnNum)
	{
		return false;
	}

	if(from>=0 && from < this->rowNum && from < this->columnNum && to >=0 && to < this->rowNum && to < this->columnNum )
	{
		double temp;
		temp = matrixNxN[from][from];

		matrixNxN[from][from] = matrixNxN[to][to];
		matrixNxN[to][to] = temp;
		return true;
	}
	return false;
};

/*
 * 将矩阵重置为单位矩阵
 * 方阵专用，非方阵返回false 并保持原状
 */
bool BigStaticMatrix::resetMatrixToI()
{
	//非方阵
	if(this->rowNum != this->columnNum)
	{
		return false;
	}

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
	return true;
};

//将矩阵重置为0矩阵
void BigStaticMatrix::resetMatrixToZero()
{
	for(int i=0; i<this->rowNum; i++)
	{
		for(int j=0; j<this->columnNum; j++)
		{

			matrixNxN[i][j] = 0;

		}
	}
};

//获取指定列向量
BasicVector* BigStaticMatrix::getColumnVector(int columnNo)
{
	if(columnNo>=0 && columnNo<this->columnNum)
	{
		int newColumnVectorDimension = this->rowNum ;
		matrixVector.resetDimension(newColumnVectorDimension);

		for(int i=0; i<this->rowNum; i++)
		{
			matrixVector.setElement(i,matrixNxN[i][columnNo]);

		}
		return &matrixVector;
	}
	else
	{
		matrixVector.resetDimension(0);
		return &matrixVector;
	}

};

//获取指定行向量
BasicVector* BigStaticMatrix::getRowVector(int rowNo)
{
	if(rowNo>=0 && rowNo<this->rowNum)
	{
		int newRowVectorDimension = this->columnNum ;
		matrixVector.resetDimension(newRowVectorDimension);

		for(int i=0; i<this->columnNum; i++)
		{
			matrixVector.setElement(i,matrixNxN[rowNo][i]);
		}
		return &matrixVector;
	}
	else
	{
		matrixVector.resetDimension(0);
		return &matrixVector;
	}
};


/*
 * 获取对角线子矩阵指定列向量
 * 仅对方阵生效，非方阵返回0元素向量
 */
BasicVector* BigStaticMatrix::getSubMatrixColumnVector(int subMatrixIndex, int columnNo)
{
	//非方阵
	if(this->rowNum != this->columnNum)
	{
		matrixVector.resetDimension(0);
		return &matrixVector;
	}

	//入参越界判断
	if(subMatrixIndex >=0 && columnNo>=0 && subMatrixIndex+columnNo < this->columnNum)
	{
		int newColumnVectorDimension = this->rowNum - subMatrixIndex;
		matrixVector.resetDimension(newColumnVectorDimension);

		int index = 0;
		for(int i=subMatrixIndex; i<this->rowNum; i++)
		{
			matrixVector.setElement(index,matrixNxN[i][subMatrixIndex+columnNo]);
			index++;
		}
		return &matrixVector;
	}
	else
	{
		matrixVector.resetDimension(0);
		return &matrixVector;
	}

};

/*
 * 获取对角线子矩阵指定行向量
 * 仅对方阵生效，非方阵返回0元素向量
 */
BasicVector* BigStaticMatrix::getSubMatrixRowVector(int subMatrixIndex, int rowNo)
{
	//非方阵
	if(this->rowNum != this->columnNum)
	{
		matrixVector.resetDimension(0);
		return &matrixVector;
	}

	//入参越界判断
	if(subMatrixIndex >=0 && rowNo>=0 && subMatrixIndex+rowNo < this->rowNum)
	{
		int newRowVectorDimension = this->columnNum - subMatrixIndex;
		matrixVector.resetDimension(newRowVectorDimension);
		int index = 0;
		for(int i=subMatrixIndex; i<this->columnNum; i++)
		{
			matrixVector.setElement(index,matrixNxN[subMatrixIndex+rowNo][i]);
			index++;
		}
		return &matrixVector;
	}
	else
	{
		matrixVector.resetDimension(0);
		return &matrixVector;
	}
};

/*
 * 获取指定对角子矩阵hessenberg列向量
 * 仅对方阵生效，非方阵返回0元素向量
 */
BasicVector* BigStaticMatrix::getSubMatrixHessenColumnVector(int subMatrixIndex)
{
	//非方阵
	if(this->rowNum != this->columnNum)
	{
		matrixVector.resetDimension(0);
		return &matrixVector;
	}
	//入参越界判断
	if(subMatrixIndex >=0 && subMatrixIndex < this->columnNum)
	{
		int newColumnVectorDimension = this->rowNum - subMatrixIndex - 1;
		matrixVector.resetDimension(newColumnVectorDimension);

		int index = 0;
		for(int i=subMatrixIndex+1; i<this->rowNum; i++)
		{
			matrixVector.setElement(index,matrixNxN[i][subMatrixIndex]);
			index++;
		}
		return &matrixVector;
	}
	else
	{
		matrixVector.resetDimension(0);
		return &matrixVector;
	}
};


