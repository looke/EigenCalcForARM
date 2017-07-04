/*
 * HessenbergFormular.cpp
 *
 *  Created on: 2017��5��8��
 *      Author: looke
 */

#include "HessenbergFormular.h"
//#include "iostream"
using namespace std;

//HessenbergFormular::HessenbergFormular()
//{};

HessenbergFormular::HessenbergFormular(BasicMatrix* p_Input_OpMatrix, BasicMatrix* p_Input_preTransMatrix,BasicMatrix* p_Input_transMatrix, BasicMatrix* p_Input_tempMatrix):m_HouseholderTrans(p_Input_OpMatrix->getColumnVector(0)),m_Multiplier(p_Input_OpMatrix,p_Input_OpMatrix,p_Input_tempMatrix),m_Transposer()
{
	this->init(p_Input_OpMatrix, p_Input_preTransMatrix, p_Input_transMatrix, p_Input_tempMatrix);
};

/*
 * Hessenberg��ʽ��
 *
 */
void HessenbergFormular::formularUpperHessnbergMatrix()
{
	this->p_preTransMatrix->resetMatrixToI();
	//this->p_afterTransMatrix->resetMatrixToI();

	if(this->p_OpMatrix->isUpperHessenbergMatrix())
	{
		return;
	}
	for(int i=0; i<this->p_OpMatrix->columnNum-2; i++)
	{
		generateSubHouseholderTrans(i);
		upgradeSubHouseholderTrans(i);
		updateOpMatrix();
		updatePreTransMatrix();
		//updateAfterTransMatrix();
	}
};

//���ݵ�ǰ����  ��ʼ���ӱ任����
//void HessenbergFormular::initSubHouseholderTrans(int iterateNum)
//{
//	int rowAndColumnNumber = this->p_OpMatrix->rowNum - iterateNum;
//	this->resizeSubMatrix(rowAndColumnNumber);
//};

/*
 * ���ݵ�ǰ��������,�����Ӿ����ά��,����任�Ӿ���
 */
void HessenbergFormular::generateSubHouseholderTrans(int iterateNum)
{
	//��ʼ���任�Ӿ���
	this->p_transMatrix->resetMatrixToI();
	int rowAndColumnNumber = this->p_OpMatrix->rowNum - iterateNum -1;
	this->p_transMatrix->resizeMatrix(rowAndColumnNumber, rowAndColumnNumber);

	BasicVector* p_hesseColumnVector = this->p_OpMatrix->getSubMatrixHessenColumnVector(iterateNum);
	//cout << "generateSubHouseholderTrans-----Current Hessen Vector:" << endl;
	//p_hesseColumnVector->printVector();

	this->m_HouseholderTrans.reload(p_hesseColumnVector);
	this->m_HouseholderTrans.getHouseholderMatrixToE1_ReverseElement(p_transMatrix);

	//cout << "generateSubHouseholderTrans-----Sub trans matrix for iterate:" << iterateNum <<endl;
	//this->p_transMatrix->printMatrix();

};


/*
 * ���ӱ任��������Ϊȫά�ȱ任��
 * ���ݵ�ǰ��������,�����Ӿ����ά��,���任�Ӿ�������
 */
void HessenbergFormular::upgradeSubHouseholderTrans(int iterateNum)
{
	//p_transMatrix->resetMatrixToI();

	//int startIndex = iterateNum+1;
	//int subTransMatrixRowColumnNumber = this->p_OpMatrix->rowNum - startIndex;

	//double temp;
	//for(int i=0; i<subTransMatrixRowColumnNumber; i++)
	//{
	//	for(int j=0; j<subTransMatrixRowColumnNumber;j++)
	//	{
	//		temp = this->p_subTransMatrix->getMatrixElement(i,j);
	//		this->p_transMatrix->setMatrixElement(startIndex+i,startIndex+j,temp);
	//	}
	//}
	p_transMatrix->resizeMatrix(this->p_OpMatrix->rowNum, this->p_OpMatrix->columnNum);
	p_transMatrix->moveDiagonalSubMatrixDown(0,this->p_OpMatrix->columnNum-iterateNum-2,iterateNum+1);
	//cout << "upgradeSubHouseholderTrans-----Complete trans matrix for iterate:" << iterateNum <<endl;
	//this->p_transMatrix->printMatrix();
};

/*
 * ������˱任����
 */
void HessenbergFormular::updatePreTransMatrix()
{
	this->m_Multiplier.reload(this->p_transMatrix, this->p_preTransMatrix, this->p_tempMatrix);
	this->m_Multiplier.multiplyCalc();
	this->p_preTransMatrix->copyMatrixElementNoCheck(p_tempMatrix);

	//cout << "updatePreTransMatrix-----Pre trans matrix for total" << endl;
	//this->p_preTransMatrix->printMatrix();

};

/*
 * �����ҳ˱任����
 */
//void HessenbergFormular::updateAfterTransMatrix()
//{
//	this->p_Multiplier->reload(this->p_afterTransMatrix, this->p_transMatrix);
//	this->p_Multiplier->multiplyCalc();
//	this->p_afterTransMatrix->copyMatrixElementNoCheck(p_Multiplier->getMultiplyResult());

//	cout << "updateAfterTransMatrix-----After trans matrix for total" << endl;
//	this->p_afterTransMatrix->printMatrix();
//};

/*
 * ���²�������  pre*OP*after
 */
void HessenbergFormular::updateOpMatrix()
{
	this->m_Multiplier.reload(this->p_transMatrix, this->p_OpMatrix, this->p_tempMatrix);
	this->m_Multiplier.multiplyCalc();
	this->p_OpMatrix->copyMatrixElementNoCheck(p_tempMatrix);
	//cout << "updateOpMatrix-----OP matrix pre tans" << endl;
	//this->p_OpMatrix->printMatrix();

	this->m_Multiplier.reload(this->p_OpMatrix, this->p_transMatrix, this->p_tempMatrix);
	this->m_Multiplier.multiplyCalc();
	this->p_OpMatrix->copyMatrixElementNoCheck(p_tempMatrix);

	//cout << "updateOpMatrix-----OP matrix after tans" << endl;
	//this->p_OpMatrix->printMatrix();
};

BasicMatrix* HessenbergFormular::getOpMatrix()
{
	return this->p_OpMatrix;
};

BasicMatrix* HessenbergFormular::getPreTransMatrix()
{
	return this->p_preTransMatrix;
};

//BasicMatrix* HessenbergFormular::getAfterTransMatrix()
//{
//	return this->p_afterTransMatrix;
//};

BasicMatrix* HessenbergFormular::getTransMatrix()
{
	return this->p_transMatrix;
};

//BasicMatrix* HessenbergFormular::getSubTransMatrix()
//{
//	return this->p_subTransMatrix;
//};

void HessenbergFormular::init(BasicMatrix* p_Input_OpMatrix, BasicMatrix* p_Input_preTransMatrix,BasicMatrix* p_Input_transMatrix, BasicMatrix* p_Input_tempMatrix)
{
	this->p_OpMatrix = p_Input_OpMatrix;
	this->p_preTransMatrix = p_Input_preTransMatrix;
	this->p_transMatrix = p_Input_transMatrix;
	this->p_tempMatrix = p_Input_tempMatrix;
};
void HessenbergFormular::reload(BasicMatrix* p_Input_OpMatrix, BasicMatrix* p_Input_preTransMatrix,BasicMatrix* p_Input_transMatrix, BasicMatrix* p_Input_tempMatrix)
{
	this->init(p_Input_OpMatrix, p_Input_preTransMatrix, p_Input_transMatrix, p_Input_tempMatrix);
};

//void HessenbergFormular::resizeSubMatrix(int rowAndColumnNumber)
//{
//	p_transMatrix->resizeMatrix(rowAndColumnNumber, rowAndColumnNumber);
//};
