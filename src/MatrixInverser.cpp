/*
 * MatrixInverser.cpp
 *
 *  Created on: 2017��2��27��
 *      Author: looke
 */

#include "MatrixInverser.h"
//#include <iostream>
using namespace std;

MatrixInverser::MatrixInverser(BasicMatrix* p_input_opMatrix,BasicMatrix* p_input_inverseMatrix)
{
	this->init(p_input_opMatrix, p_input_inverseMatrix);
};

/*
 * 	��ʼ���������������
 */
void MatrixInverser::init(BasicMatrix* p_input_opMatrix,BasicMatrix* p_input_inverseMatrix)
{
	this->isFullRank = false;
	this->isSquareMatrix = false;

	this->rowNumber = p_input_opMatrix->rowNum;
	this->columnNumber = p_input_opMatrix->columnNum;


	//ֻ�з����������
	if(rowNumber == columnNumber)
	{
		isSquareMatrix = true;
	}

	this->p_operateMatrix = p_input_opMatrix;
	this->p_inverseMatrix = p_input_inverseMatrix;

};
/*
 * ����װ��������󣬳�ʼ���������
 * ���������������
 */
void MatrixInverser::reload(BasicMatrix* p_input_opMatrix,BasicMatrix* p_input_inverseMatrix)
{
	this->init(p_input_opMatrix, p_input_inverseMatrix);
};

/*
 * ��ӡ�����
 */
//void MatrixInverser::printInverseMatrix()
//{
//	this->p_inverseMatrix->printMatrix();
//};

//void MatrixInverser::printOpMatrix()
//{
//	this->p_operateMatrix->printMatrix();
//};



/*
 * ���������
 * ��������ʽ�任
 */
void MatrixInverser::generateInverseMatrix()
{
	//�Ƿ��󲻽�������������
	if(!this->isSquareMatrix)
	{
		return;
	}
	//cout << "Start Invers" << endl;


	//�Բ���������й�һ���Լ���Ԫ �������ͬ������
	for(int i=0; i<this->rowNumber; i++)
	{
		this->unZeroPivotRow(i);
		this->normalizePivotRow(i);
		this->eliminateSubMatrix(i);
		//cout << "After normalize and sub eliminate: " << i << endl;
		//this->p_operateMatrix->printMatrix();
		//cout << "-----------------------" << endl;
		//this->p_inverseMatrix->printMatrix();
	}



	for(int i=this->rowNumber-1; i>0; i--)
	{
		this->eliminateUpperMatrix(i);
		//cout << "After up eliminate: " << i << endl;
		//this->p_operateMatrix->printMatrix();
		//cout << "-----------------------" << endl;
		//this->p_inverseMatrix->printMatrix();
	}

	//�ж��Ƿ�����
	int rank = 0;
	for(int i=0; i<this->rowNumber; i++)
	{
		if(0 != this->p_operateMatrix->getMatrixElement(i,i))
		{
			rank++;
		}
	}
	if(rank == this->rowNumber)
	{
		this->isFullRank = true;
	}

};

/*
 * �ж�ָ���еĶԽ�����Ԫ�Ƿ�Ϊ0,��Ϊ0�򽻻���������Ԫ��0��
 *
 */
void MatrixInverser::unZeroPivotRow(int pivotRow)
{
	//���ȶԾ����������,�����׶Խ�����ԪΪ0
	if(0 == this->p_operateMatrix->getMatrixElement(pivotRow,pivotRow))
	{
		int i=pivotRow;
		while(i<this->rowNumber)
		{
			if(0 != this->p_operateMatrix->getMatrixElement(i,pivotRow))
			{
				break;
			}
			i++;
		}
		this->p_operateMatrix->swapRow(0,i);
		this->p_inverseMatrix->swapRow(0,i);
	}
	//cout << "After swap" << endl;
	//this->p_operateMatrix->printMatrix();
	//cout << "-----------------------" << endl;
	//this->p_inverseMatrix->printMatrix();
};


/*
 * ��һ��ָ���У�����Ӧ��Ϊ����Ԫ������
 *
 * ͬʱҲ�ᴦ��Э����
 * pivotRow��0��ʼ
 */
void MatrixInverser::normalizePivotRow(int pivotRow)
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
	//double pivotValue = this->opMatrixPointer[pivotRow][pivotRow];
	double pivotValue = this->p_operateMatrix->getMatrixElement(pivotRow, pivotRow);
	//�жϷ�ĸΪ0�����
	if(pivotValue == 0)
	{
		//������һ��
		return;
	}

	for(int i=pivotRow;i<this->columnNumber;i++)
	{
		//�������������
		//this->opMatrixPointer[pivotRow][i] = this->opMatrixPointer[pivotRow][i]/pivotValue;
		p_operateMatrix->setMatrixElement(pivotRow, i, p_operateMatrix->getMatrixElement(pivotRow, i)/pivotValue);

	}
	//����Э����
	for(int j=0;j<this->columnNumber;j++)
	{
		//this->invMatrixPointer[pivotRow][j] = this->invMatrixPointer[pivotRow][j]/pivotValue;
		p_inverseMatrix->setMatrixElement(pivotRow, j, p_inverseMatrix->getMatrixElement(pivotRow, j)/pivotValue);

	}

};


/*
 * ��ָ������Ԫ�����п�ʼ�������·����н�����Ԫ
 * �˺���Ӧ���ڵ����й�һ���Ժ����
 *
 * pivotRow��0��ʼ
 */
void MatrixInverser::eliminateSubMatrix(int pivotRow)

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
	double currentPivotValue;
	for(int i=pivotRow+1; i<this->rowNumber; i++)
	{
		//currentPivotValue = opMatrixPointer[i][pivotRow];
		currentPivotValue = p_operateMatrix->getMatrixElement(i, pivotRow);
		for(int j=0; j<this->columnNumber; j++)
		{
			//opMatrixPointer[i][j] = opMatrixPointer[i][j] - currentPivotValue*opMatrixPointer[pivotRow][j];
			p_operateMatrix->setMatrixElement(i,j, p_operateMatrix->getMatrixElement(i,j) - currentPivotValue*p_operateMatrix->getMatrixElement(pivotRow,j));

			//invMatrixPointer[i][j] = invMatrixPointer[i][j]- currentPivotValue*invMatrixPointer[pivotRow][j];
			p_inverseMatrix->setMatrixElement(i,j,p_inverseMatrix->getMatrixElement(i,j) - currentPivotValue*p_inverseMatrix->getMatrixElement(pivotRow,j));
		}
	}
};


/*
 * ��ָ���п�ʼ�������Ϸ����н�����Ԫ
 * �˺���Ӧ���ڵ����Ӿ�����Ԫ�Ժ����
 *
 * startRow�Ӿ�����ĩ�п�ʼ
 */
void MatrixInverser::eliminateUpperMatrix(int startRow)
{
	double temp;
	//��ָ���п�ʼ���������Ǿ��������Ԫ�����������Խ�Ԫ��
	for(int i=startRow-1; i >= 0; i--)
	{
		//temp = this->opMatrixPointer[i][startRow];
		temp = p_operateMatrix->getMatrixElement(i,startRow);
		//this->opMatrixPointer[i][startRow] = 0;
		p_operateMatrix->setMatrixElement(i,startRow,0);
		for(int j=0; j<this->columnNumber; j++)
		{
			//this->invMatrixPointer[i][j] = invMatrixPointer[i][j] - this->invMatrixPointer[startRow][j]*temp;
			p_inverseMatrix->setMatrixElement(i,j, p_inverseMatrix->getMatrixElement(i,j) - p_inverseMatrix->getMatrixElement(startRow,j)*temp);
		}
	}
};

//BasicMatrix* MatrixInverser::getInverseMatrix()
//{
//	return this->p_inverseMatrix;
//};
bool MatrixInverser::isOpMatrixFullRank()
{
	return this->isFullRank;
};
