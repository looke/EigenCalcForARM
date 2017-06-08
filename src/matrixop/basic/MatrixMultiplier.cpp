/*
 * MatrixMultiplier.cpp
 *
 *  Created on: 2017��2��25��
 *      Author: looke
 */

#include "..\include\matrixop\basic\MatrixMultiplier.h"
#include <iostream>
using namespace std;

MatrixMultiplier::MatrixMultiplier(BasicMatrix* leftOp, BasicMatrix* rightOp, BasicMatrix* resultOp)
{
	this->init(leftOp, rightOp, resultOp);
};

void MatrixMultiplier::init(BasicMatrix* leftOp, BasicMatrix* rightOp, BasicMatrix* resultOp)  throw(length_error)
{

};
/*
void MatrixMultiplier::init(BasicMatrix* leftOp, BasicMatrix* rightOp)
{
	this->leftOpMatrix = leftOp;
	this->leftRow = leftOp->rowNum;
	this->leftColumn = leftOp->columnNum;

	this->rightOpMatrix = rightOp;
	this->rightRow = rightOp->rowNum;
	this->rightColumn = rightOp->columnNum;

	if(this->leftColumn == this->rightRow)
	{
		this->isLegalMultiply = true;
		this->MultiResult = DynamicMatrix(this->leftRow, this->rightColumn);
	}
	else
	{
		this->isLegalMultiply = false;
	}
};
*/
void MatrixMultiplier::reload(BasicMatrix* leftOp, BasicMatrix* rightOp, BasicMatrix* resultOp)
{

};
/*
 * ����װ��
 * ʹ���µ����Ҳ����������³�ʼ���˷���
 * ����µ����Ҳ����������������䣬�򲻻��������������������Խ�ʡ�ڴ濪��
 */
/*
void MatrixMultiplier::reload(BasicMatrix* leftOp, BasicMatrix* rightOp)
{
	//�������װ��ľ����뵱ǰ��ʼ���ľ���������һ�£�����Ҫ�������ɳ˻�������󣬽�ԭ�˻���������������ü��ɣ������ڴ濪��
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
	else//�������װ��ľ����뵱ǰ��ʼ���ľ�����������ͬ������Ҫ�������ɳ˻��������
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
	//	//���˷�������ʱ���ؿ�ָ��
	//	return 0;
	//}
//};

void MatrixMultiplier::multiplyCalc()
{
	//if(isLegalMultiply)
	//{
		double temp = 0.0;
		double leftOp = 0.0;
		double rightOp = 0.0;
		double multiTemp = 0.0;
		for(int i=0; i<this->leftRow; i++)
		{
			for(int j=0; j<this->rightColumn; j++)
			{
				temp=0;

				for(int k=0; k<this->leftColumn; k++ )
				{
					leftOp = p_leftOpMatrix->getMatrixElement(i,k);
					rightOp = p_rightOpMatrix->getMatrixElement(k,j);
					multiTemp = leftOp*rightOp;
					temp = temp + multiTemp;
				}
				p_MultiResult->setMatrixElement(i,j, temp);
			}
		}
	//}
};


void MatrixMultiplier::printMultiplyResult()
{
	//if(isLegalMultiply)
	//{
		p_MultiResult->printMatrix();

	//}
};

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
