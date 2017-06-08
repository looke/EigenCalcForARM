/*
 * StaticMatrixMultiplier.cpp
 *
 *  Created on: 2017��2��25��
 *      Author: looke
 */

#include "..\include\matrixop\static\StaticMatrixMultiplier.h"

using namespace std;

StaticMatrixMultiplier::StaticMatrixMultiplier(BasicMatrix* leftOp, BasicMatrix* rightOp, BasicMatrix* resultOp):MatrixMultiplier(leftOp,rightOp,resultOp)
{
	this->init(leftOp, rightOp, resultOp);
};

void StaticMatrixMultiplier::init(BasicMatrix* leftOp, BasicMatrix* rightOp, BasicMatrix* resultOp) throw(length_error)
{
	if(leftOp->columnNum == rightOp->rowNum && leftOp->rowNum == resultOp->rowNum && rightOp->columnNum == resultOp->columnNum)
	{
		this->p_leftOpMatrix = leftOp;
		this->leftRow = leftOp->rowNum;
		this->leftColumn = leftOp->columnNum;

		this->p_rightOpMatrix = rightOp;
		this->rightRow = rightOp->rowNum;
		this->rightColumn = rightOp->columnNum;

		this->p_MultiResult = resultOp;
	}
	else
	{
		throw length_error(getMatrixLengthErrorMessage(leftOp, rightOp, resultOp));
	}
};

/*
 * ����װ��
 * ʹ���µ����Ҳ����������³�ʼ���˷���
 * ����µ����Ҳ����������������䣬�򲻻��������������������Խ�ʡ�ڴ濪��
 */
void StaticMatrixMultiplier::reload(BasicMatrix* leftOp, BasicMatrix* rightOp, BasicMatrix* resultOp)
{
	//�������װ��ľ����뵱ǰ��ʼ���ľ���������һ�£�����Ҫ�������ɳ˻�������󣬽�ԭ�˻���������������ü��ɣ������ڴ濪��
	//if(	this->leftRow == leftOp->rowNum &&
	//	this->leftColumn == leftOp->columnNum &&
	//	this->rightRow == rightOp->rowNum &&
	//	this->rightColumn == rightOp->columnNum
	//  )
	//{
	//	this->p_leftOpMatrix = leftOp;
	//	this->p_rightOpMatrix = rightOp;
		//this->MultiResult.resetMatrixToI();
	//}
	//else//�������װ��ľ����뵱ǰ��ʼ���ľ�����������ͬ������Ҫ�������ɳ˻��������
	//{
		this->init(leftOp, rightOp, resultOp);
	//}
};


