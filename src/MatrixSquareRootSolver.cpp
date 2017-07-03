/*
 * MatrixSquareRootSolver.cpp
 *
 *  Created on: 2017��3��11��
 *      Author: looke
 */
#include "MatrixSquareRootSolver.h"
#include <iostream>
using namespace std;
/*
 * ��ƽ����������ʼ����������ξ��������
 * �ж��Ƿ���
 * ��ʼ��������
 * ��ʼ������
 * ���ɰ��浥λ����Z
 */
MatrixSquareRootSolver::MatrixSquareRootSolver(BasicMatrix* input_opYMatrix, BasicMatrix* input_opZMatrix, BasicMatrix* input_YMatrix_Temp, BasicMatrix* input_ZMatrix_Temp, BasicMatrix* input_tempMatrix)
:m_inverser(input_YMatrix_Temp,input_ZMatrix_Temp)
{
	this->init(input_opYMatrix,input_opZMatrix,input_YMatrix_Temp,input_ZMatrix_Temp, input_tempMatrix);
};

void MatrixSquareRootSolver::init(BasicMatrix* input_opYMatrix, BasicMatrix* input_opZMatrix, BasicMatrix* input_YMatrix_Temp, BasicMatrix* input_ZMatrix_Temp, BasicMatrix* input_tempMatrix)
{
	//��ʼ��������
	this->rowNumber = input_opYMatrix->rowNum;
	this->columnNumber = input_opYMatrix->columnNum;
	this->maxDiff = 10000;
	this->threshold = 0.01;	//����Ҫ��
	this->maxIterationTime = 30; //������30��
	this->iterationTime = 0; //��ǰ������

	//�ж��Ƿ���
	if(this->rowNumber == this->columnNumber)
	{
		this->isSquareMatrix = true;
	}
	else
	{
		this->isSquareMatrix = false;
	}

	//��ʼ����������
	this->p_YMatrix = input_opYMatrix;

	//���ɰ��浥λ����Z
	this->p_ZMatrix = input_opZMatrix;
	//this->p_ZMatrix->resetMatrixToI();

	//��ʼ���м����ֵ����
	this->p_tempYMatrix = input_YMatrix_Temp;
	//this->p_tempYMatrix->resetMatrixToI();

	this->p_tempZMatrix = input_ZMatrix_Temp;
	//this->p_tempZMatrix->resetMatrixToI();

	this->p_tempMatrix = input_tempMatrix;
};

void MatrixSquareRootSolver::reload(BasicMatrix* input_opYMatrix, BasicMatrix* input_opZMatrix, BasicMatrix* input_YMatrix_Temp, BasicMatrix* input_ZMatrix_Temp, BasicMatrix* input_tempMatrix)
{
	this->init(input_opYMatrix,input_opZMatrix,input_YMatrix_Temp,input_ZMatrix_Temp, input_tempMatrix);
};
/*
 * ���ݲ�������Y�Ͱ������Z�����һ��D-B����
 * �����������������ʱ����tempYMatrix��tempZMatrix
 */
void MatrixSquareRootSolver::DenmanBeaversIteration()
{
	double tempValue;

	//����Y�������
	p_tempMatrix->copyMatrixElementNoCheck(p_YMatrix);
	p_tempZMatrix->resetMatrixToI();
	this->m_inverser.reload(p_tempMatrix,p_tempZMatrix);
	this->m_inverser.generateInverseMatrix();

	//����Y���������ʱ����tempZMatrix
	//p_tempZMatrix->copyMatrixElementNoCheck(m_inverser->getInverseMatrix());
	/*
	for(int i=0; i<this->rowNumber;i++)
	{
		for(int j=0; j<this->columnNumber;j++)
		{
			this->tempZMatrix.setMatrixElement(i,j,tempInvers[i][j]);
		}
	}
	*/
	/////���Դ�ӡ���
	cout << "Y Invers:" << endl;
	this->p_tempZMatrix->printMatrix();

	//����Z�������
	p_tempMatrix->copyMatrixElementNoCheck(p_ZMatrix);
	p_tempYMatrix->resetMatrixToI();
	this->m_inverser.reload(p_tempMatrix,p_tempYMatrix);
	this->m_inverser.generateInverseMatrix();

	//tempInvers = this->inver.getInverseMatrix();
	//����Z���������ʱ����tempYMatrix
	//p_tempYMatrix->copyMatrixElementNoCheck(this->p_inverser->getInverseMatrix());

	/*
	for(int i=0; i<this->rowNumber;i++)
	{
		for(int j=0; j<this->columnNumber;j++)
		{
			this->tempYMatrix.setMatrixElement(i,j,tempInvers[i][j]);
		}
	}
	*/
	cout << "Z Invers:" << endl;
	this->p_tempYMatrix->printMatrix();

	//�������֮�Ͳ��������Ԫ�س���0.5�����DB����
	for(int i=0; i<this->rowNumber;i++)
	{
		for(int j=0; j<this->columnNumber;j++)
		{
			//���� Yn+1����ĸ���Ԫ�� ��Zn�����Yn֮���ٳ���0.5
			//tempValue = this->tempYMatrixPointer[i][j] + this->YMatrixPointer[i][j];
			tempValue = this->p_tempYMatrix->getMatrixElement(i,j) + this->p_YMatrix->getMatrixElement(i,j);
			tempValue = 0.5*tempValue;
			this->p_tempYMatrix->setMatrixElement(i,j,tempValue);

			//���� Zn+1����ĸ���Ԫ�� ��Yn�����Zn֮���ٳ���0.5
			//tempValue = this->tempZMatrixPointer[i][j] + this->ZMatrixPointer[i][j];
			tempValue = this->p_tempZMatrix->getMatrixElement(i,j) + this->p_ZMatrix->getMatrixElement(i,j);
			tempValue = 0.5*tempValue;
			this->p_tempZMatrix->setMatrixElement(i,j,tempValue);
		}
	}
	cout << "Yn+1:" << endl;
	this->p_tempYMatrix->printMatrix();

	cout << "Zn+1:" << endl;
	this->p_tempZMatrix->printMatrix();

	//����������Yn+1��Yn��Ԫ������ֵ
	//this->maxDiff = calcMaxDifferential(this->YMatrix, this->tempYMatrix);
	this->maxDiff = this->p_YMatrix->calcMaxDifferentialNoCheck(this->p_tempYMatrix);
};


/*
 * ��ӡY���󣬾���������������ƽ����
 */
void MatrixSquareRootSolver::printYMatrix()
{
	for(int i=0;i<this->rowNumber;i++)
	{
		for(int j=0;j<this->columnNumber;j++)
		{
			//cout<<this->YMatrixPointer[i][j]<<"\t";
			cout << this->p_YMatrix->getMatrixElement(i,j) << "\t";
		}
		cout<<endl;
	}
};

/*
 * ��ӡZ���󣬾���������������ƽ�����������
 */
void MatrixSquareRootSolver::printZMatrix()
{
	for(int i=0;i<this->rowNumber;i++)
	{
		for(int j=0;j<this->columnNumber;j++)
		{
			//cout<<this->ZMatrixPointer[i][j]<<"\t";
			cout << this->p_ZMatrix->getMatrixElement(i,j) << "\t";
		}
		cout<<endl;
	}
};

bool MatrixSquareRootSolver::generateSquareRootMatrix()
{
	bool isSuccess = true;
	this->iterationTime=0;

	//���ɰ��浥λ����Z
	this->p_ZMatrix->resetMatrixToI();

	//��ʼ���м����ֵ����
	this->p_tempYMatrix->resetMatrixToI();
	this->p_tempZMatrix->resetMatrixToI();
	this->p_tempMatrix->resetMatrixToI();

	while(this->maxDiff > this->threshold)
	{
		this->iterationTime++;
		//���������������˳� ��ָʾ���ʧ��
		if(this->iterationTime > this->maxIterationTime)
		{
			isSuccess = false;
			break;
		}

		cout << "D-B Iteration:" << this->iterationTime << endl;
		this->DenmanBeaversIteration();
		//����Yn+1��Yn  Zn+1��Zn
		p_YMatrix->copyMatrixElementNoCheck(p_tempYMatrix);
		p_ZMatrix->copyMatrixElementNoCheck(p_tempZMatrix);
		/*
		for(int i=0;i<this->rowNumber;i++)
		{
			for(int j=0;j<this->columnNumber;j++)
			{
				this->YMatrix.setMatrixElement(i,j, this->tempYMatrixPointer[i][j]);
				this->YMatrix.setMatrixElement(i,j, this->tempYMatrix.getMatrixElement(i,j));


				this->ZMatrix.setMatrixElement(i,j, this->tempZMatrixPointer[i][j]);
				this->ZMatrix.setMatrixElement(i,j, this->tempYMatrix.getMatrixElement(i,j));
			}
		}
		*/
	}

	return isSuccess;
};

int MatrixSquareRootSolver::getIterationTime()
{
	return this->iterationTime;
};

