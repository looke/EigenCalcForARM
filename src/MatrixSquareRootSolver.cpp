/*
 * MatrixSquareRootSolver.cpp
 *
 *  Created on: 2017年3月11日
 *      Author: looke
 */
#include "MatrixSquareRootSolver.h"
//#include <iostream>
using namespace std;
/*
 * 求平方根函数初始化，复制入参矩阵的内容
 * 判断是否方阵
 * 初始化行列数
 * 初始化矩阵
 * 生成伴随单位矩阵Z
 */
MatrixSquareRootSolver::MatrixSquareRootSolver(BasicMatrix* input_opYMatrix, BasicMatrix* input_opZMatrix, BasicMatrix* input_YMatrix_Temp, BasicMatrix* input_ZMatrix_Temp, BasicMatrix* input_tempMatrix)
:m_inverser(input_YMatrix_Temp,input_ZMatrix_Temp)
{
	this->init(input_opYMatrix,input_opZMatrix,input_YMatrix_Temp,input_ZMatrix_Temp, input_tempMatrix);
};

void MatrixSquareRootSolver::init(BasicMatrix* input_opYMatrix, BasicMatrix* input_opZMatrix, BasicMatrix* input_YMatrix_Temp, BasicMatrix* input_ZMatrix_Temp, BasicMatrix* input_tempMatrix)
{
	//初始化行列数
	this->rowNumber = input_opYMatrix->rowNum;
	this->columnNumber = input_opYMatrix->columnNum;
	this->maxDiff = 10000;
	this->threshold = 0.01;	//精度要求
	this->maxIterationTime = 30; //最多迭代30次
	this->iterationTime = 0; //当前迭代数

	//判断是否方阵
	if(this->rowNumber == this->columnNumber)
	{
		this->isSquareMatrix = true;
	}
	else
	{
		this->isSquareMatrix = false;
	}

	//初始化操作矩阵
	this->p_YMatrix = input_opYMatrix;

	//生成伴随单位矩阵Z
	this->p_ZMatrix = input_opZMatrix;
	//this->p_ZMatrix->resetMatrixToI();

	//初始化中间计算值矩阵
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
 * 根据操作矩阵Y和伴随矩阵Z，完成一次D-B迭代
 * 将迭代结果保存在临时矩阵tempYMatrix和tempZMatrix
 */
void MatrixSquareRootSolver::DenmanBeaversIteration()
{
	double tempValue;

	//生成Y的逆矩阵
	p_tempMatrix->copyMatrixElementNoCheck(p_YMatrix);
	p_tempZMatrix->resetMatrixToI();
	this->m_inverser.reload(p_tempMatrix,p_tempZMatrix);
	this->m_inverser.generateInverseMatrix();

	//拷贝Y的逆矩阵到临时矩阵tempZMatrix
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
	/////测试打印输出
	//cout << "Y Invers:" << endl;
	//this->p_tempZMatrix->printMatrix();

	//生成Z的逆矩阵
	p_tempMatrix->copyMatrixElementNoCheck(p_ZMatrix);
	p_tempYMatrix->resetMatrixToI();
	this->m_inverser.reload(p_tempMatrix,p_tempYMatrix);
	this->m_inverser.generateInverseMatrix();

	//tempInvers = this->inver.getInverseMatrix();
	//拷贝Z的逆矩阵到临时矩阵tempYMatrix
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
	//cout << "Z Invers:" << endl;
	//this->p_tempYMatrix->printMatrix();

	//计算矩阵之和并将矩阵各元素乘以0.5，完成DB迭代
	for(int i=0; i<this->rowNumber;i++)
	{
		for(int j=0; j<this->columnNumber;j++)
		{
			//计算 Yn+1矩阵的各个元素 由Zn的逆和Yn之和再乘以0.5
			//tempValue = this->tempYMatrixPointer[i][j] + this->YMatrixPointer[i][j];
			tempValue = this->p_tempYMatrix->getMatrixElement(i,j) + this->p_YMatrix->getMatrixElement(i,j);
			tempValue = 0.5*tempValue;
			this->p_tempYMatrix->setMatrixElement(i,j,tempValue);

			//计算 Zn+1矩阵的各个元素 由Yn的逆和Zn之和再乘以0.5
			//tempValue = this->tempZMatrixPointer[i][j] + this->ZMatrixPointer[i][j];
			tempValue = this->p_tempZMatrix->getMatrixElement(i,j) + this->p_ZMatrix->getMatrixElement(i,j);
			tempValue = 0.5*tempValue;
			this->p_tempZMatrix->setMatrixElement(i,j,tempValue);
		}
	}
	//cout << "Yn+1:" << endl;
	//this->p_tempYMatrix->printMatrix();

	//cout << "Zn+1:" << endl;
	//this->p_tempZMatrix->printMatrix();

	//计算迭代后的Yn+1与Yn的元素最大差值
	//this->maxDiff = calcMaxDifferential(this->YMatrix, this->tempYMatrix);
	this->maxDiff = this->p_YMatrix->calcMaxDifferentialNoCheck(this->p_tempYMatrix);
};


/*
 * 打印Y矩阵，经过迭代会趋近于平方根
 */
/*
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
*/
/*
 * 打印Z矩阵，经过迭代会趋近于平方根的逆矩阵
 */
/*
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
*/
bool MatrixSquareRootSolver::generateSquareRootMatrix()
{
	bool isSuccess = true;
	this->iterationTime=0;

	//生成伴随单位矩阵Z
	this->p_ZMatrix->resetMatrixToI();

	//初始化中间计算值矩阵
	this->p_tempYMatrix->resetMatrixToI();
	this->p_tempZMatrix->resetMatrixToI();
	this->p_tempMatrix->resetMatrixToI();

	while(this->maxDiff > this->threshold)
	{
		this->iterationTime++;
		//超过迭代次数后退出 并指示求解失败
		if(this->iterationTime > this->maxIterationTime)
		{
			isSuccess = false;
			break;
		}

		//cout << "D-B Iteration:" << this->iterationTime << endl;
		this->DenmanBeaversIteration();
		//拷贝Yn+1到Yn  Zn+1到Zn
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

