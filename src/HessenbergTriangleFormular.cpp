/*
 * HessenbergTriangleFormular.cpp
 *
 *  Created on: 2017年4月28日
 *      Author: looke
 */

#include "HessenbergTriangleFormular.h"
#include <iostream>
using namespace std;

//HessenbergTriangleFormular::HessenbergTriangleFormular()
//{};

HessenbergTriangleFormular::HessenbergTriangleFormular(BasicMatrix* input_Matrix_A, BasicMatrix* input_Matrix_B,BasicMatrix* input_QMatrix_Total,BasicMatrix* input_ZMatrix_Total,BasicMatrix* input_QZMatrix_Step,BasicMatrix* input_TempMatrix_Trans,BasicMatrix* input_TempMatrix):m_Transposer(),m_Multiplier(input_Matrix_A, input_Matrix_B,input_TempMatrix),m_GivensTrans(input_Matrix_A->getColumnVector(0)),m_QRDecomp(input_Matrix_B, input_QZMatrix_Step, input_TempMatrix_Trans, input_TempMatrix)
{
	this->init(input_Matrix_A, input_Matrix_B, input_QMatrix_Total, input_ZMatrix_Total, input_QZMatrix_Step, input_TempMatrix_Trans, input_TempMatrix);
};

void HessenbergTriangleFormular::init(BasicMatrix* input_Matrix_A, BasicMatrix* input_Matrix_B,BasicMatrix* input_QMatrix_Total,BasicMatrix* input_ZMatrix_Total,BasicMatrix* input_QZMatrix_Step,BasicMatrix* input_TempMatrix_Trans,BasicMatrix* input_TempMatrix)
{
	//原始操作方阵A
	this->p_OpMatrix_A = input_Matrix_A;
	//原始操作方阵B
	this->p_OpMatrix_B = input_Matrix_B;

	//总体Q变换矩阵,用于汇总左乘操作矩阵
	this->p_QMatrix_Total = input_QMatrix_Total;

	//总体Z变换矩阵,用于汇总右乘操作矩阵
	this->p_ZMatrix_Total = input_ZMatrix_Total;

	//单步QZ变换矩阵,用于左乘操作矩阵
	this->p_QZMatrix_Step = input_QZMatrix_Step;

	//中间转换矩阵
	this->p_TempMatrix_Trans = input_TempMatrix_Trans;

	//中间过程矩阵
	this->p_TempMatrix = input_TempMatrix;
};

void HessenbergTriangleFormular::reload(BasicMatrix* input_Matrix_A, BasicMatrix* input_Matrix_B,BasicMatrix* input_QMatrix_Total,BasicMatrix* input_ZMatrix_Total,BasicMatrix* input_QZMatrix_Step,BasicMatrix* input_TempMatrix_Trans,BasicMatrix* input_TempMatrix)
{
	this->init(input_Matrix_A, input_Matrix_B, input_QMatrix_Total, input_ZMatrix_Total, input_QZMatrix_Step, input_TempMatrix_Trans, input_TempMatrix);
};

BasicMatrix* HessenbergTriangleFormular::getMatrixA()
{
	return this->p_OpMatrix_A;
};

BasicMatrix* HessenbergTriangleFormular::getMatrixB()
{
	return this->p_OpMatrix_B;
};
BasicMatrix* HessenbergTriangleFormular::getHessenbergMatrixA()
{
	return this->p_OpMatrix_A;
};
BasicMatrix* HessenbergTriangleFormular::getTriangleMatrixB()
{
	return this->p_OpMatrix_B;
};
BasicMatrix* HessenbergTriangleFormular::getMatrixQ_Total()
{
	return this->p_QMatrix_Total;
};

BasicMatrix* HessenbergTriangleFormular::getMatrixZ_Total()
{
	return this->p_ZMatrix_Total;
};

//BasicMatrix* HessenbergTriangleFormular::getMatrixQ_Step()
//{
//	return this->p_QMatrix_Step;
//};

//BasicMatrix* HessenbergTriangleFormular::getMatrixZ_Step()
//{
//	return this->p_ZMatrix_Step;
//};

BasicMatrix* HessenbergTriangleFormular::getMatrixQZ_Step()
{
	return this->p_QZMatrix_Step;
};

/*
 * 形成A-Hessenberg B-Triangle 矩阵对
 */
void HessenbergTriangleFormular::formularABMatrix()
{
	p_QMatrix_Total->resetMatrixToI();
	p_ZMatrix_Total->resetMatrixToI();
	//p_QMatrix_Step->resetMatrixToI();
	p_QZMatrix_Step->resetMatrixToI();
	if(this->p_OpMatrix_A->isUpperHessenbergMatrix() && this->p_OpMatrix_B->isUpperTriangleMatrix())
	{
		//this->p_OpHessenbergMatrix_A->copyMatrixElementNoCheck(p_OpMatrix_A);
		//this->p_OpTriangleMatrix_B->copyMatrixElementNoCheck(p_OpMatrix_B);
		return;
	}

	this->initABMatrix();
	cout << "After HessenBergTriangle init--- A:" << endl;
	p_OpMatrix_A->printMatrix();
	cout << "After HessenBergTriangle init--- B:" << endl;
	p_OpMatrix_B->printMatrix();
	cout << "After HessenBergTriangle init--- Q Total:" << endl;
	p_QMatrix_Total->printMatrix();
	cout << "After HessenBergTriangle init--- Z Total:" << endl;
	p_ZMatrix_Total->printMatrix();

	for(int i=0; i<this->p_OpMatrix_A->columnNum; i++)
	{
		this->formularColumnVector(i);
		//p_OpHessenbergMatrix_A->regularZeroElement();
		//p_OpTriangleMatrix_B->regularZeroElement();
	}
};


/*
 * 对A矩阵的指定列进行hessenberg格式化
 */
void HessenbergTriangleFormular::formularColumnVector(int columnIndex)
{
	double element,lowEdge;
	BasicVector* p_ColumnVector;
	//BasicMatrix* p_GivensMatrix;
	//for(int i=this->p_OpMatrix_A->rowNum-1; i>columnIndex+1; i--)

	int i=this->p_OpMatrix_A->rowNum-1;

	while(i>columnIndex+1)
	{
		p_QZMatrix_Step->resetMatrixToI();

		lowEdge = p_OpMatrix_A->getLowEdge();
		//先计算更新A矩阵
		element = p_OpMatrix_A->getMatrixElementRegulared(i,columnIndex,lowEdge);
		if(0 != element)
		{
			p_ColumnVector = p_OpMatrix_A->getColumnVector(columnIndex);

			this->m_GivensTrans.reload(p_ColumnVector);
			this->m_GivensTrans.getGivensMatrixPreMultiple(i, p_QZMatrix_Step);
			cout << "After HessenBergTriangle--- Left Givens Matrix:" << i << endl;
			p_QZMatrix_Step->printMatrix();
		}


		updateOpMatrix_A_ByQ();
		cout << "formularColumnVector: columnIndex:" <<columnIndex <<"OP Matrix A By Q" << endl;
		p_OpMatrix_A->printMatrix();
		updateOpMatrix_B_ByQ();
		cout << "formularColumnVector: columnIndex:" <<columnIndex <<"OP Matrix B By Q" << endl;
		p_OpMatrix_B->printMatrix();
		updateQMatrix_Total();

		p_QZMatrix_Step->resetMatrixToI();

		//接着计算更新B矩阵
		lowEdge = p_OpMatrix_B->getLowEdge();
		element = p_OpMatrix_B->getMatrixElementRegulared(i,i-1,lowEdge);
		if(0 != element)
		{
			p_ColumnVector = p_OpMatrix_B->getRowVector(i);
			cout << "p_OpTriangleMatrix_B:getRowVector:" << i << endl;
			p_ColumnVector->printVector();
			this->m_GivensTrans.reload(p_ColumnVector);
			//this->p_GivensTrans->setIsUsingPreElement(false);
			this->m_GivensTrans.getGivensMatrixAfterMultiple(i-1, p_QZMatrix_Step);
			cout << "After HessenBergTriangle--- Right Givens Matrix:" << i << endl;
			p_QZMatrix_Step->printMatrix();
		}

		updateOpMatrixByZ();
		updateZMatrix_Total();

		cout << "After HessenBergTriangle--- A:" << i << endl;
		p_OpMatrix_A->printMatrix();
		cout << "After HessenBergTriangle--- B:" << i<< endl;
		p_OpMatrix_B->printMatrix();
		cout << "After HessenBergTriangle--- Q:" << i << endl;
		p_QMatrix_Total->printMatrix();
		cout << "After HessenBergTriangle--- Z:" << i << endl;
		p_ZMatrix_Total->printMatrix();

		//确保当前要消去的元素值可以被忽略
		if(0 == p_OpMatrix_A->getMatrixElementRegulared(i,columnIndex,lowEdge))
		{
			i--;
		}
	}
};


/*
 * 使用QT矩阵 更新操作矩阵
 *  QT * A
 *  QT * B
 */
void HessenbergTriangleFormular::updateOpMatrix_A_ByQ()
{
	//更新操作矩阵A
	m_Multiplier.reload(p_QZMatrix_Step, p_OpMatrix_A, p_TempMatrix);
	m_Multiplier.multiplyCalc();
	this->p_OpMatrix_A->copyMatrixElementNoCheck(p_TempMatrix);
	//this->p_OpMatrix_A->regularZeroElement();
};

void HessenbergTriangleFormular::updateOpMatrix_B_ByQ()
{
	//更新操作矩阵B
	m_Multiplier.reload(p_QZMatrix_Step, p_OpMatrix_B, p_TempMatrix);
	m_Multiplier.multiplyCalc();
	this->p_OpMatrix_B->copyMatrixElementNoCheck(p_TempMatrix);
	//this->p_OpTriangleMatrix_B->regularZeroElement();
};
/*
 * 使用Z矩阵 更新操作矩阵
 *  A * Z
 *  B * Z
 */
void HessenbergTriangleFormular::updateOpMatrixByZ()
{
	//更新操作矩阵A
	m_Multiplier.reload(p_OpMatrix_A, p_QZMatrix_Step, p_TempMatrix);
	m_Multiplier.multiplyCalc();
	this->p_OpMatrix_A->copyMatrixElementNoCheck(p_TempMatrix);
	//this->p_OpHessenbergMatrix_A->regularZeroElement();

	//更新操作矩阵B
	m_Multiplier.reload(p_OpMatrix_B, p_QZMatrix_Step, p_TempMatrix);
	m_Multiplier.multiplyCalc();
	this->p_OpMatrix_B->copyMatrixElementNoCheck(p_TempMatrix);
	//this->p_OpTriangleMatrix_B->regularZeroElement();

};

void HessenbergTriangleFormular::updateQMatrix_Total()
{
	//更新Q矩阵
	m_Multiplier.reload(p_QZMatrix_Step, p_QMatrix_Total, p_TempMatrix);
	m_Multiplier.multiplyCalc();
	this->p_QMatrix_Total->copyMatrixElementNoCheck(p_TempMatrix);
};

void HessenbergTriangleFormular::updateZMatrix_Total()
{
	//更新Z矩阵
	m_Multiplier.reload(p_ZMatrix_Total, p_QZMatrix_Step, p_TempMatrix);
	m_Multiplier.multiplyCalc();
	this->p_ZMatrix_Total->copyMatrixElementNoCheck(p_TempMatrix);
};

/*
 * 原始操作矩阵B 进行QR分解 化为上三角矩阵，并将QT左乘A
 */
void HessenbergTriangleFormular::initABMatrix()
{
	//this->p_OpHessenbergMatrix_A->copyMatrixElementNoCheck(this->p_OpMatrix_A);
	//this->p_OpTriangleMatrix_B->copyMatrixElementNoCheck(this->p_OpMatrix_B);

	//对B矩阵进行QR分解，
	this->m_QRDecomp.reload(p_OpMatrix_B, p_QZMatrix_Step, p_TempMatrix_Trans, p_TempMatrix);
	this->m_QRDecomp.calcQRMatrix();
	//Q变换矩阵为QR分解出的Q矩阵的转置
	m_Transposer.transposeSquareMatrix(p_QZMatrix_Step);

	//this->p_QMatrix_Step->copyMatrixElementNoCheck(this->p_QRDecomp->getQTMatrix());
	updateOpMatrix_A_ByQ();
	updateQMatrix_Total();

	//Z变换矩阵为单位阵
	//this->p_QZMatrix_Step->resetMatrixToI();
	//updateOpMatrixByZ();
	//updateZMatrix_Total();
};
