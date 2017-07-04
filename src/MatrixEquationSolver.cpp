/*
 * MatrixEquationSolver.cpp
 *
 *  Created on: 2017��2��19��
 *      Author: looke
 */

#include "MatrixEquationSolver.h"
//#include <iostream>
//#include <malloc.h>
#include "math.h"

using namespace std;

MatrixEquationSolver::MatrixEquationSolver(BasicMatrix* input_opMatrix)
{
	this->init(input_opMatrix);
}

void MatrixEquationSolver::init(BasicMatrix* input_opMatrix)
{
	this->p_opMatrix = input_opMatrix;
	this->rowNumber = input_opMatrix->rowNum;
	this->columnNumber = input_opMatrix->columnNum;
};

void MatrixEquationSolver::reload(BasicMatrix* input_opMatrix)
{
	this->init(input_opMatrix);
};

/**
 *  ��˹����Ԫ��Ԫ�� ���й��̺���
 *  ��������ָ������
 *  from��0�п�ʼ
 *
 */
void MatrixEquationSolver::swapRow(int from, int to)
{
	p_opMatrix->swapRow(from,to);

}

/*
 * ��˹����Ԫ��Ԫ�� ������Ԫ�����й��̺���
 *
 * ��������Ԫ������
 *
 * ��startRowָ���п�ʼ�������кŴ��ڵ���ָ���е�����Ԫ��startRow��0��ʼ
 * ��������Ԫ�����к�
 *
 */
int MatrixEquationSolver::findPivotRow(int startRow)
{
	if(startRow >= this->rowNumber)
	{
		//cout<<"Illegal row number." << endl;
		return -1;
	}
	if(startRow >= this->columnNumber)
	{
		//cout<<"Pivot not exist in this sub matrix" << endl;
		return -1;
	}

	int pivotRow = startRow;
	double maxElement = 0.0;

	for(int i=startRow; i<this->rowNumber; i++)
	{
		if(i >= this->columnNumber)
		{
			break;
		}
		if(fabs(p_opMatrix->getMatrixElement(i,startRow)) > fabs(maxElement))
		{
			maxElement = p_opMatrix->getMatrixElement(i,startRow);
			pivotRow = i;
		}
	}

	return pivotRow;

}

/*
 * ��˹����Ԫ��Ԫ�� ��һ�����̺���
 *
 * ��һ��ָ���У�����Ӧ��Ϊ����Ԫ������
 * �˺���Ӧ���ڵ���findPivotRow�ҵ�����Ԫ�����в����任���ϲ��Ժ����
 *
 * pivotRow��0��ʼ
 */
void MatrixEquationSolver::normalizePivotRow(int pivotRow)
{
	if(pivotRow >= this->rowNumber)
	{
		//cout<<"Illegal row number." << endl;
		return ;
	}
	if(pivotRow >= this->columnNumber)
	{
		//cout<<"Can not normalize sub matrix." << endl;
		//cout<<"PivotRow:" << pivotRow <<endl;
		//cout<<"Column number:" << this->columnNum <<endl;
		return ;
	}

	//ȡ����Ԫ��ֵ��Ϊ��ĸ
	double pivotValue = p_opMatrix->getMatrixElement(pivotRow,pivotRow);

	//�жϷ�ĸΪ0�����
	if(pivotValue == 0)
	{
		//������һ��
		return;
	}

	for(int i=pivotRow;i<this->columnNumber;i++)
	{
		p_opMatrix->setMatrixElement(pivotRow, i, p_opMatrix->getMatrixElement(pivotRow,i)/pivotValue);
	}
};

/*
 * ��˹����Ԫ��Ԫ�� ��Ԫ���̺���
 *
 * ��ָ������Ԫ�����п�ʼ�������·����н�����Ԫ
 * �˺���Ӧ���ڵ����й�һ���Ժ����
 *
 * pivotRow��0��ʼ
 */
void MatrixEquationSolver::eliminateSubMatrix(int pivotRow)

{
	if(pivotRow >= this->rowNumber)
	{
		//cout<<"Illegal row number." << endl;
		return ;
	}
	if(pivotRow >= this->columnNumber)
	{
		//cout<<"Can not eliminate sub matrix." << endl;
		//cout<<"PivotRow:" << pivotRow <<endl;
		//cout<<"Column number:" << this->columnNum <<endl;
		return ;
	}
	for(int i=pivotRow+1;i<this->rowNumber;i++)
	{
		double currentPivotValue = p_opMatrix->getMatrixElement(i,pivotRow);
		for(int j=pivotRow; j<this->columnNumber; j++)
		{
			//matrixPointer[i][j] = matrixPointer[i][j] - currentPivotValue*matrixPointer[pivotRow][j];
			p_opMatrix->setMatrixElement(i, j,
			p_opMatrix->getMatrixElement(i,j) - currentPivotValue*p_opMatrix->getMatrixElement(pivotRow,j));
		}
	}
};

/*
 * ��˹����Ԫ��Ԫ�� ת���Խ�����̺���
 *
 * ������Ϊ�Խ�����Ԫ��Ϊ1���϶Խ���
 *
 */
void MatrixEquationSolver::gaussElim_ColmnPrin()
{
	int pivotRow;
	for(int i=0; i<this->rowNumber; i++)
	{
		//cout << "LOOP:" << i << endl;

		pivotRow = this->findPivotRow(i);
		//cout << "PivotRow:" << pivotRow << endl;

		this->swapRow(i, pivotRow);
		//cout << "After swap"<< endl;
		//this->p_opMatrix->printMatrix();

		this->normalizePivotRow(i);
		//cout << "After Normalization" << endl;
		//this->p_opMatrix->printMatrix();

		this->eliminateSubMatrix(i);
		//cout << "After elimination" << endl;
		//this->p_opMatrix->printMatrix();
	}
};

/*
 * ��˹����Ԫ��Ԫ�� �ж�������ָ���е�Ԫ���Ƿ�ȫΪ0���̺���
 *
 * �ж�������ָ���е�Ԫ���Ƿ�ȫΪ0
 * arguRowNo��0��ʼ
 */
bool MatrixEquationSolver::isAllZeroForArgumentRow(int arguRowNo)
{
	if(arguRowNo >= this->rowNumber)
	{
		//cout << "Illegal input row No." << endl;
		return false;
	}
	int length = this->columnNumber;

	int zeroNumber = 0;

	for(int i=0; i<length; i++)
	{
		if(p_opMatrix->getMatrixElement(arguRowNo,i) == 0)
		{
			zeroNumber++;
		}
	}

	if(zeroNumber == length)
	{
		return true;
	}
	else
	{
		return false;
	}
};

/*
 * ��˹����Ԫ��Ԫ�� �ж�ϵ����ָ���е�Ԫ���Ƿ�ȫΪ0���̺���
 *
 * �ж�ϵ����ָ���е�Ԫ���Ƿ�ȫΪ0
 * coRowNo��0��ʼ
 */
bool MatrixEquationSolver::isAllZeroForCoRow(int coRowNo)
{
	if(coRowNo >= this->rowNumber)
	{
			//cout << "Illegal input row No." << endl;
			return false;
	}
	int length = this->columnNumber-1;

	int zeroNumber = 0;

	for(int i=0; i<length; i++)
	{
		if(p_opMatrix->getMatrixElement(coRowNo,i) == 0)
		{
			zeroNumber++;
		}
	}

	if(zeroNumber == length)
	{
		return true;
	}
	else
	{
		return false;
	}

};

/*
 * ��˹����Ԫ��Ԫ�� ϵ�������ȹ��̺���
 *
 * ����ϵ���������
 *
 * ��������Ԫ���̽��������
 */
int MatrixEquationSolver::coefficientMatrixRank()
{
	int rank = 0;
	for(int i=0; i<this->rowNumber; i++)
	{
		if(!isAllZeroForCoRow(i))
		{
			rank++;
		}
	}
	return rank;
};

/*
 * ��˹����Ԫ��Ԫ�� ���������ȹ��̺���
 *
 * ��������������
 *
 * ��������Ԫ���̽��������
 */
int MatrixEquationSolver::augmentedMatrixRank()
{
	int rank = 0;
	for(int i=0; i<this->rowNumber; i++)
	{
		if(!isAllZeroForArgumentRow(i))
		{
			rank++;
		}
	}
	return rank;
};


/*
 * ��˹����Ԫ��Ԫ�� �������
 *
 * �������ĸ�
 *
 * ��������Ԫ���̽��������
 */
bool MatrixEquationSolver::calcRoots()
{
	if(this->rowNumber < this->columnNumber-1)
	{
		//cout << "Do not have enough row to solve." << endl;
		return false;
	}

	int augRank = this->augmentedMatrixRank();
	int coRank = this->coefficientMatrixRank();
	//cout << "Argument Rank of Matrix:" << arguRank << " . Co Rank of Matrix:" << coRank <<endl;
	if(augRank > coRank)
	{
		//cout << "No root exist." <<endl;
		return false;
	}

	if(augRank == coRank)
	{
		//cout << "Only one root exist." <<endl;
	}

	//����ĩ�п�ʼ���������Ǿ��������Ԫ�����������Խ�Ԫ��
	for(int i=this->rowNumber-1; i > 0; i--)
	{
		double temp = p_opMatrix->getMatrixElement(i,columnNumber-1);

		for(int j=i-1; j>=0; j--)
		{
			//matrixPointer[j][this->columnNumber-1] = matrixPointer[j][this->columnNumber-1] - matrixPointer[j][i]*temp;
			p_opMatrix->setMatrixElement(j, this->columnNumber-1,
			p_opMatrix->getMatrixElement(j,columnNumber-1) - p_opMatrix->getMatrixElement(j,i)*temp);
			p_opMatrix->setMatrixElement(j,i,0);
		}
	}
	return true;
	//for(int i=0; i < rowNumber; i++)
	//{
	//	this->roots[i] = opMatrix.getMatrixElement(i,columnNumber-1);
	//}
};

/*
 * ��˹����Ԫ��Ԫ�� ��ȡ������
 *
 * ��ȡ����ĸ�ֵ
 *
 * ������������̽��������
 */
//double* MatrixEquationSolver::getRoots()
//{
//	return this->roots;
//};


/*
 * ��˹����Ԫ��Ԫ�� �ⷽ����
 *
 * �����󷽳���
 */
bool MatrixEquationSolver::solveMatrixEquation()
{
	this->gaussElim_ColmnPrin();
	return this->calcRoots();
};
