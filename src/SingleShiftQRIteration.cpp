/*
 * SingelShiftQRIteration.cpp
 *
 *  Created on: 2017年5月13日
 *      Author: looke
 */

#include "SingleShiftQRIteration.h"
//#include <iostream>
//using namespace std;

//SingleShiftQRIteration::SingleShiftQRIteration()
//{};

SingleShiftQRIteration::SingleShiftQRIteration(BasicMatrix* p_input_OpMatrix, BasicMatrix* p_input_QTMatrix_Implicit_Total,BasicMatrix* p_input_Q_QT_Matrix_Implicit_Step,BasicMatrix* p_input_TempMatrix):m_Transposer(),m_Multiplier(p_input_OpMatrix,p_input_OpMatrix,p_input_TempMatrix),m_GivensTrans(p_input_OpMatrix->getColumnVector(0)),m_HessenbergForm(p_input_OpMatrix, p_input_QTMatrix_Implicit_Total, p_input_Q_QT_Matrix_Implicit_Step, p_input_TempMatrix),m_QRDecomp(p_input_OpMatrix, p_input_QTMatrix_Implicit_Total, p_input_Q_QT_Matrix_Implicit_Step, p_input_TempMatrix)
{
	this->init(p_input_OpMatrix, p_input_QTMatrix_Implicit_Total, p_input_Q_QT_Matrix_Implicit_Step, p_input_TempMatrix);
};

void SingleShiftQRIteration::init(BasicMatrix* p_input_OpMatrix,BasicMatrix* p_input_QTMatrix_Implicit_Total,BasicMatrix* p_input_Q_QT_Matrix_Implicit_Step,BasicMatrix* p_input_TempMatrix)
{
	//操作矩阵
	this->p_OpMatrix = p_input_OpMatrix;

	//Q 矩阵 隐式迭代 综合 Q用于右乘OP矩阵
	//this->p_QMatrix_Implicit_Total = p_input_QMatrix_Implicit_Total;
	//QT 矩阵 隐式迭代 综合 QT用于左乘OP矩阵
	this->p_QTMatrix_Implicit_Total = p_input_QTMatrix_Implicit_Total;

	//Q/QT 矩阵 隐式迭代 分步 Q用于右乘OP矩阵 QT用于左乘OP矩阵
	this->p_Q_QT_Matrix_Implicit_Step = p_input_Q_QT_Matrix_Implicit_Step;

	//中间过程矩阵
	this->p_TempMatrix = p_input_TempMatrix;

	generateHessenbergOpMatrix();
};

void SingleShiftQRIteration::reload(BasicMatrix* p_input_OpMatrix, BasicMatrix* p_input_QTMatrix_Implicit_Total,BasicMatrix* p_input_Q_QT_Matrix_Implicit_Step,BasicMatrix* p_input_TempMatrix)
{
	this->init(p_input_OpMatrix, p_input_QTMatrix_Implicit_Total, p_input_Q_QT_Matrix_Implicit_Step, p_input_TempMatrix);
};

/*
 * 生成hessenberg操作矩阵
 */
void SingleShiftQRIteration::generateHessenbergOpMatrix()
{
	m_HessenbergForm.reload(p_OpMatrix, p_QTMatrix_Implicit_Total, p_Q_QT_Matrix_Implicit_Step, p_TempMatrix);
	m_HessenbergForm.formularUpperHessnbergMatrix();
	//this->p_OpHessenbergMatrix->copyMatrixElementNoCheck(p_HessenbergForm->getOpMatrix());
	//p_QMatrix_Implicit_Total->copyMatrixElementNoCheck(p_QTMatrix_Implicit_Total);
	//m_Transposer.transposeSquareMatrix(p_QMatrix_Implicit_Total);
};


/*
 * 单值位移QR迭代 显式
 */
/*
void SingleShiftQRIteration::explicit_QRIteration(double input_shiftValue)
{
	int iterationNumber = 10;

	//this->generateHessenbergOpMatrix();
	cout << "SingelShiftQRIteration--explicit_QRIteration----OP Hessenberg Matrix" << endl;
	this->p_OpHessenbergMatrix->printMatrix();

	int i=0;
	while (i<iterationNumber)
	{
		explicit_QRIteration_Step(input_shiftValue);
		cout << "SingelShiftQRIteration--explicit_QRIteration----OP Hessenberg Matrix after:" << i << "iteration" <<endl;
		p_OpHessenbergMatrix->printMatrix();
		i++;
	}
};
*/
/*
//单值位移QR迭代 显式单步
void SingleShiftQRIteration::explicit_QRIteration_Step(double input_shiftValue)
{
	//对角元减去移动值
	p_OpHessenbergMatrix->diagonalSubtraction(input_shiftValue);
	cout << "SingelShiftQRIteration--explicit_QRIteration_Step----OP Hessenberg Matrix after Subtraction" <<endl;
	p_OpHessenbergMatrix->printMatrix();

	this->p_QRDecomp->reload(p_OpHessenbergMatrix);
	this->p_QRDecomp->calcQRMatrix();
	this->p_QMatrix_Explicit->copyMatrixElementNoCheck(this->p_QRDecomp->getQMatrix());
	p_OpHessenbergMatrix->copyMatrixElementNoCheck(this->p_QRDecomp->getRMatrix());

	cout << "SingelShiftQRIteration--explicit_QRIteration_Step----OP Hessenberg Matrix after QRDecomp" <<endl;
	p_OpHessenbergMatrix->printMatrix();

	this->p_Multiplier->reload(p_OpHessenbergMatrix, p_QMatrix_Explicit);
	this->p_Multiplier->multiplyCalc();
	p_OpHessenbergMatrix->copyMatrixElementNoCheck(p_Multiplier->getMultiplyResult());

	cout << "SingelShiftQRIteration--explicit_QRIteration_Step----OP Hessenberg Matrix after R*Q" <<endl;
	p_OpHessenbergMatrix->printMatrix();
	//对角元加上移动值
	p_OpHessenbergMatrix->diagonalAddition(input_shiftValue);
	cout << "SingelShiftQRIteration--explicit_QRIteration_Step----OP Hessenberg Matrix after Addition" <<endl;
	p_OpHessenbergMatrix->printMatrix();
};
*/

//单值位移QR迭代 隐式 初始化
void SingleShiftQRIteration::initForImplicitQR(double input_shiftValue)
{
	this->p_Q_QT_Matrix_Implicit_Step->resetMatrixToI();
	//单步移位
	double temp_old = this->p_OpMatrix->getMatrixElement(0,0);
	double temp_new = temp_old - input_shiftValue;
	this->p_OpMatrix->setMatrixElement(0,0,temp_new);


	BasicVector* p_firstColumnVector = this->p_OpMatrix->getColumnVector(0);
	//还原OpMatrix
	this->p_OpMatrix->setMatrixElement(0,0,temp_old);

	//生成隐式Q及QT
	//生成隐式QT 左乘矩阵
	m_GivensTrans.reload(p_firstColumnVector);
	m_GivensTrans.getGivensMatrixPreMultiple(1,p_Q_QT_Matrix_Implicit_Step);

	//更新操作矩阵 A-Hess
	this->updateOpMatrix_By_QT_IM_QRIteration();
	//更新总体转换矩阵QT
	this->updateQTMatrix_Total_IM_QRIteration();

	//cout << "SingelShiftQRIteration--initForImplicitQR----QTMatrix_Implicit_Step" <<endl;
	//p_Q_QT_Matrix_Implicit_Step->printMatrix();

	//this->p_QTMatrix_Implicit_Step->copyMatrixElementNoCheck(p_GivensTrans->getGivensMatrixPreMultiple(1));

	//生成隐式Q 右乘矩阵
	this->m_Transposer.transposeSquareMatrix(p_Q_QT_Matrix_Implicit_Step);

	//更新操作矩阵 A-Hess
	this->updateOpMatrix_By_Q_IM_QRIteration();
	//更新总体转换矩阵Q
	//this->updateQMatrix_Total_IM_QRIteration();

	//cout << "SingelShiftQRIteration--initForImplicitQR----QMatrix_Implicit_Step" <<endl;
	//p_Q_QT_Matrix_Implicit_Step->printMatrix();

	//cout << "SingelShiftQRIteration--initForImplicitQR----QTMatrix_Implicit_Total" <<endl;
	//p_QTMatrix_Implicit_Total->printMatrix();
	//cout << "SingelShiftQRIteration--initForImplicitQR----QMatrix_Implicit_Total" <<endl;
	//p_QMatrix_Implicit_Total->printMatrix();
	//此处已经形成隐式单步QR迭代的初始化操作矩阵
	//cout << "SingelShiftQRIteration--initForImplicitQR----OP Hessenberg Matrix after init"<<endl;
	//p_OpMatrix->printMatrix();
};

//隐式QR迭代 更新hessenberg操作矩阵 QT*A
void SingleShiftQRIteration::updateOpMatrix_By_QT_IM_QRIteration()
{
	//计算QT * OP-Hessenberg
	this->m_Multiplier.reload(p_Q_QT_Matrix_Implicit_Step, p_OpMatrix,p_TempMatrix);
	this->m_Multiplier.multiplyCalc();
	p_OpMatrix->copyMatrixElementNoCheck(p_TempMatrix);

};

//隐式QR迭代 更新hessenberg操作矩阵 A*Q
void SingleShiftQRIteration::updateOpMatrix_By_Q_IM_QRIteration()
{
	//计算 OP-Hessenberg * Q
	this->m_Multiplier.reload(p_OpMatrix, p_Q_QT_Matrix_Implicit_Step, p_TempMatrix);
	this->m_Multiplier.multiplyCalc();
	p_OpMatrix->copyMatrixElementNoCheck(p_TempMatrix);
};

//单值位移QR迭代 隐式
void SingleShiftQRIteration::implicit_QRIteration(double input_shiftValue)
{
	//操作矩阵转换为Hessenberg矩阵
	//this->generateHessenbergOpMatrix();
	//cout << "SingelShiftQRIteration--initForImplicitQR----OP Hessenberg Matrix" << endl;
	//this->p_OpMatrix->printMatrix();

	int iteratNumber = 10;
	int i = 0;
	while(i < iteratNumber)
	{
		implicit_QRIteration_Step(input_shiftValue);
		//cout << "SingelShiftQRIteration--implicit_QRIteration----OP Hessenberg Matrix after: " << i << " iteration" <<endl;
		//p_OpMatrix->printMatrix();
		i++;
	}
};

//单值位移QR迭代 隐式
void SingleShiftQRIteration::implicit_QRIteration_Step(double input_shiftValue)
{
	//p_QTMatrix_Implicit_Total->resetMatrixToI();
	//p_QMatrix_Implicit_Total->resetMatrixToI();
	initForImplicitQR(input_shiftValue);

	int bulgeNumber = this->p_OpMatrix->columnNum-2;
	for(int j=0; j<bulgeNumber; j++)
	{
		//生成隐式Q及QT
		BasicVector* p_columnVector = this->p_OpMatrix->getColumnVector(j);
		m_GivensTrans.reload(p_columnVector);
		//计算QT
		m_GivensTrans.getGivensMatrixPreMultiple(j+2,p_Q_QT_Matrix_Implicit_Step);

		//更新A-Hessenberg矩阵
		this->updateOpMatrix_By_QT_IM_QRIteration();
		//更新总体转换矩阵 QT
		this->updateQTMatrix_Total_IM_QRIteration();

		//计算Q
		m_Transposer.transposeSquareMatrix(p_Q_QT_Matrix_Implicit_Step);
		//this->p_Transposer->transposeMatrix();
		//this->->copyMatrixElementNoCheck(this->p_Transposer->getTransposeMatrix());
		//更新A-Hessenberg矩阵
		this->updateOpMatrix_By_Q_IM_QRIteration();

		//更新总体转换矩阵Q
		//this->updateQMatrix_Total_IM_QRIteration();
	}
};

/*
 * 瑞利商位移QR迭代 显式
 */
/*
void SingleShiftQRIteration::rayleigh_Quotient_EX_QRIteration()
{
	int iterationNumber = 10;

	//this->generateHessenbergOpMatrix();
	cout << "SingelShiftQRIteration--rayleigh_Quotient_EX_QRIteration----OP Hessenberg Matrix" << endl;
	this->p_OpHessenbergMatrix->printMatrix();

	int rayleighValueIndex = p_OpHessenbergMatrix->rowNum - 1;
	double rayleighValue;
	int i=0;
	while (i<iterationNumber)
	{
		rayleighValue = p_OpHessenbergMatrix->getMatrixElement(rayleighValueIndex,rayleighValueIndex);
		explicit_QRIteration_Step(rayleighValue);
		cout << "SingelShiftQRIteration--rayleigh_Quotient_EX_QRIteration----OP Hessenberg Matrix after:" << i << "iteration" <<endl;
		p_OpHessenbergMatrix->printMatrix();
		i++;
	}
};
*/

//瑞利商位移QR迭代 隐式 10次迭代
void SingleShiftQRIteration::rayleigh_Quotient_IM_QRIteration()
{
	this->rayleigh_Quotient_IM_QRIteration(10);
}

//瑞利商位移QR迭代 隐式 指定迭代次数
void SingleShiftQRIteration::rayleigh_Quotient_IM_QRIteration(int iterateNum)
{
	int iterationNumber = iterateNum;

	//this->generateHessenbergOpMatrix();
	//cout << "SingelShiftQRIteration--rayleigh_Quotient_IM_QRIteration----OP Hessenberg Matrix" << endl;
	//this->p_OpMatrix->printMatrix();

	int rayleighValueIndex = p_OpMatrix->rowNum - 1;
	double rayleighValue;
	int i=0;
	while (i<iterationNumber)
	{
		rayleighValue = p_OpMatrix->getMatrixElement(rayleighValueIndex,rayleighValueIndex);
		implicit_QRIteration_Step(rayleighValue);
		//cout << "SingelShiftQRIteration--rayleigh_Quotient_IM_QRIteration----OP Hessenberg Matrix after:" << i << "iteration" <<endl;
		//p_OpMatrix->printMatrix();
		i++;
	}
};

//瑞利商位移QR迭代 隐式 单步--每次计算都重置Total转置矩阵
void SingleShiftQRIteration::rayleigh_Quotient_IM_QRIteration_Step()
{
	//this->generateHessenbergOpMatrix();
	//cout << "SingelShiftQRIteration--rayleigh_Quotient_IM_QRIteration Step----OP Hessenberg Matrix" << endl;
	//this->p_OpMatrix->printMatrix();

	int rayleighValueIndex = p_OpMatrix->rowNum - 1;
	double rayleighValue= p_OpMatrix->getMatrixElement(rayleighValueIndex,rayleighValueIndex);
	p_QTMatrix_Implicit_Total->resetMatrixToI();
	implicit_QRIteration_Step(rayleighValue);
	//cout << "SingelShiftQRIteration--rayleigh_Quotient_IM_QRIteration Step----OP Hessenberg Matrix" <<endl;
	//p_OpMatrix->printMatrix();
};

/*
 * 更新综合转换矩阵QT Total
 */
void SingleShiftQRIteration::updateQTMatrix_Total_IM_QRIteration()
{
	//计算左乘综合矩阵 QT_Step * QT_Total
	this->m_Multiplier.reload(p_Q_QT_Matrix_Implicit_Step, p_QTMatrix_Implicit_Total,p_TempMatrix);
	this->m_Multiplier.multiplyCalc();
	p_QTMatrix_Implicit_Total->copyMatrixElementNoCheck(p_TempMatrix);
};

/*
 * 更新综合转换矩阵Q Total
 */
//void SingleShiftQRIteration::updateQMatrix_Total_IM_QRIteration()
//{
	//计算右乘综合矩阵 Q_Total * Q_Step
//	this->m_Multiplier.reload(p_QMatrix_Implicit_Total, p_Q_QT_Matrix_Implicit_Step, p_TempMatrix);
//	this->m_Multiplier.multiplyCalc();
//	p_QMatrix_Implicit_Total->copyMatrixElementNoCheck(p_TempMatrix);
//};

//BasicMatrix* SingleShiftQRIteration::getOpHessenbergMatrix()
//{
//	return this->p_OpHessenbergMatrix;
//};

/*
 * 获取总体转换矩阵 QT
 */
BasicMatrix* SingleShiftQRIteration::getQTMatrix_Total()
{
	return this->p_QTMatrix_Implicit_Total;
};

/*
 * 获取总体转换矩阵 Q
 */
//BasicMatrix* SingleShiftQRIteration::getQMatrix_Total()
//{
//	return this->p_QMatrix_Implicit_Total;
//};
