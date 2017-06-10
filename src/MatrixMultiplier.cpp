/*
 * MatrixMultiplier.cpp
 *
 *  Created on: 2017年2月25日
 *      Author: looke
 */

#include "MatrixMultiplier.h"
#include <iostream>
using namespace std;

MatrixMultiplier::MatrixMultiplier(BasicMatrix* leftOp, BasicMatrix* rightOp, BasicMatrix* resultOp)
{
	this->init(leftOp, rightOp, resultOp);
};

void MatrixMultiplier::init(BasicMatrix* leftOp, BasicMatrix* rightOp, BasicMatrix* resultOp)
{
	this->p_leftOpMatrix = leftOp;
	this->p_rightOpMatrix = rightOp;
	this->p_MultiResult = resultOp;
};

void MatrixMultiplier::reload(BasicMatrix* leftOp, BasicMatrix* rightOp, BasicMatrix* resultOp)
{
	this->init(leftOp, rightOp, resultOp);
};
/*
 * 重新装填
 * 使用新的左右操作矩阵，重新初始化乘法器
 * 如果新的左右操作矩阵行列数不变，则不会重新申请操作结果矩阵，以节省内存开销
 */
/*
void MatrixMultiplier::reload(BasicMatrix* leftOp, BasicMatrix* rightOp)
{
	//如果重新装填的矩阵与当前初始化的矩阵行列数一致，则不需要重新生成乘积结果矩阵，将原乘积结果矩阵内容重置即可，降低内存开销
	if(	this->leftRow == leftOp->rowNum &&
		this->leftColumn == leftOp->columnNum &&
		this->rightRow == rightOp->rowNum &&
		this->rightColumn == rightOp->columnNum
	  )
	{
		this->leftOpMatrix = leftOp;
		this->rightOpMatrix = rightOp;
		this->MultiResult.resetMatrixToI();
	}
	else//如果重新装填的矩阵与当前初始化的矩阵行列数不同，则需要重新生成乘积结果矩阵
	{
		this->init(leftOp, rightOp);
	}
};
*/


//BasicMatrix* MatrixMultiplier::getMultiplyResult()
//{
	//if(isLegalMultiply)
	//{
//		return this->p_MultiResult;
	//}
	//else
	//{
	//	//当乘法不成立时返回空指针
	//	return 0;
	//}
//};

bool MatrixMultiplier::multiplyCalc()
{
	if(p_leftOpMatrix->columnNum == p_rightOpMatrix->rowNum && p_leftOpMatrix->rowNum == p_MultiResult->rowNum && p_rightOpMatrix->columnNum == p_MultiResult->columnNum)
	{
		double temp = 0.0;
		double leftOp = 0.0;
		double rightOp = 0.0;
		double multiTemp = 0.0;
		for(int i=0; i<p_leftOpMatrix->rowNum; i++)
		{
			for(int j=0; j<p_rightOpMatrix->columnNum; j++)
			{
				temp=0;

				for(int k=0; k<p_leftOpMatrix->columnNum; k++ )
				{
					leftOp = p_leftOpMatrix->getMatrixElement(i,k);
					rightOp = p_rightOpMatrix->getMatrixElement(k,j);
					multiTemp = leftOp*rightOp;
					temp = temp + multiTemp;
				}
				p_MultiResult->setMatrixElement(i,j, temp);
			}
		}
		return true;
	}
	else
		return false;
};

//打印乘积结果矩阵
void MatrixMultiplier::printMultiplyResult()
{
	p_MultiResult->printMatrix();
};
/*
string MatrixMultiplier::getMatrixLengthErrorMessage(BasicMatrix* leftOp, BasicMatrix* rightOp, BasicMatrix* resultOp)
{
	stringstream stream;

	stream<<leftOp->rowNum;
	string left_Row = stream.str();
	stream.str("");
	stream<<leftOp->columnNum;
	string left_Column = stream.str();

	stream.str("");
	stream<<rightOp->rowNum;
	string right_Row = stream.str();
	stream.str("");
	stream<<rightOp->columnNum;
	string right_Column = stream.str();

	stream.str("");
	stream<<resultOp->rowNum;
	string result_Row = stream.str();
	stream.str("");
	stream<<resultOp->columnNum;
	string result_Column = stream.str();

	string exInfo = "Matrix length Error. LeftOp:" +left_Row + "*" + left_Column + " RightOp:" + right_Row + "*" + right_Column + " ResultOp:" + result_Row + "*" + result_Column + ".";
	return exInfo;
};
*/
