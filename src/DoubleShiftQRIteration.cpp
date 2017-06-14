/*
 * DoubleShiftQRIteration.cpp
 *
 *  Created on: 2017年5月17日
 *      Author: looke
 */

#include "DoubleShiftQRIteration.h"
#include <iostream>
using namespace std;

//DoubleShiftQRIteration::DoubleShiftQRIteration()
//{};

DoubleShiftQRIteration::DoubleShiftQRIteration(BasicMatrix* p_input_OpMatrix,BasicVector* p_input_TransVector,BasicMatrix* p_input_QMatrix_Total,BasicMatrix* p_input_QTMatrix_Total,BasicMatrix* p_input_QQTMatrix_Step,BasicMatrix* p_input_TempMatrix)
:m_Transposer(),m_Multiplier(p_input_OpMatrix,p_input_OpMatrix,p_input_TempMatrix),m_GivensTrans(p_input_TransVector),m_HouseholderTrans(p_input_TransVector),m_HessenbergForm(p_input_OpMatrix,p_input_QTMatrix_Total,p_input_QQTMatrix_Step,p_input_TempMatrix)
{
	this->init(p_input_OpMatrix,p_input_TransVector,p_input_QMatrix_Total,p_input_QTMatrix_Total,p_input_QQTMatrix_Step,p_input_TempMatrix);
};

void DoubleShiftQRIteration::init(BasicMatrix* p_input_OpMatrix,BasicVector* p_input_TransVector,BasicMatrix* p_input_QMatrix_Total,BasicMatrix* p_input_QTMatrix_Total,BasicMatrix* p_input_QQTMatrix_Step,BasicMatrix* p_input_TempMatrix)
{
	//操作矩阵
	this->p_OpMatrix = p_input_OpMatrix;

	//用于计算隐式双步QR跌代单步转换矩阵Qi/Qn-1 的向量---3维/2维
	this->p_TransVectorForQStep = p_input_TransVector;

	//Q 矩阵 隐式迭代  总体Q用于右乘OP矩阵
	this->p_QMatrix_Implicit_Total = p_input_QMatrix_Total;
	//QT 矩阵 隐式迭代 总体QT用于左乘OP矩阵
	this->p_QTMatrix_Implicit_Total = p_input_QTMatrix_Total;

	//Q/QT 矩阵 隐式迭代 分步 QT用于右乘/左乘OP矩阵
	this->p_QQTMatrix_Implicit_Step = p_input_QQTMatrix_Step;

	//中间过程矩阵
	this->p_TempMatrix = p_input_TempMatrix;

	this->generateHessenbergOpMatrix();
};

void DoubleShiftQRIteration::reload(BasicMatrix* p_input_OpMatrix,BasicVector* p_input_TransVector,BasicMatrix* p_input_QMatrix_Total,BasicMatrix* p_input_QTMatrix_Total,BasicMatrix* p_input_QQTMatrix_Step,BasicMatrix* p_input_TempMatrix)
{
	this->init(p_input_OpMatrix,p_input_TransVector,p_input_QMatrix_Total,p_input_QTMatrix_Total,p_input_QQTMatrix_Step,p_input_TempMatrix);
};


/*
 * 生成hessenberg操作矩阵
 */
void DoubleShiftQRIteration::generateHessenbergOpMatrix()
{
	m_HessenbergForm.reload(this->p_OpMatrix, this->p_QTMatrix_Implicit_Total,this->p_QQTMatrix_Implicit_Step,this->p_TempMatrix);
	m_HessenbergForm.formularUpperHessnbergMatrix();
	//this->p_OpHessenbergMatrix->copyMatrixElementNoCheck(m_HessenbergForm.getOpMatrix());
	//此时p_OpMatrix已变成上Hessenberg矩阵
	p_QMatrix_Implicit_Total->copyMatrixElementNoCheck(p_QTMatrix_Implicit_Total);
	this->m_Transposer.transposeSquareMatrix(p_QMatrix_Implicit_Total);
};

/*
 * 根据hessenberg矩阵最后2x2子矩阵 计算生成特征值的和以及乘积
 */
void DoubleShiftQRIteration::generateWilkinsonShift()
{
	int lastRowIndex = this->p_OpMatrix->rowNum - 1;

	double An_1n_1 = this->p_OpMatrix->getMatrixElement(lastRowIndex-1,lastRowIndex-1);
	double Ann = this->p_OpMatrix->getMatrixElement(lastRowIndex,lastRowIndex);

	double An_1n = this->p_OpMatrix->getMatrixElement(lastRowIndex-1,lastRowIndex);
	double Ann_1 = this->p_OpMatrix->getMatrixElement(lastRowIndex,lastRowIndex-1);

	this->trace = An_1n_1 + Ann;
	this->determinant = An_1n_1*Ann - An_1n*Ann_1;
};

/*
 * B = (A-pI)*(A-tI); 根据hessenberg矩阵 生成双步移位后的矩阵第一列
 */
void DoubleShiftQRIteration::generateBMatrixFirstColumn()
{
	generateWilkinsonShift();
	double a11, a12, a21, a22, a32;
	a11 = this->p_OpMatrix->getMatrixElement(0,0);
	a12 = this->p_OpMatrix->getMatrixElement(0,1);
	a21 = this->p_OpMatrix->getMatrixElement(1,0);
	a22 = this->p_OpMatrix->getMatrixElement(1,1);
	a32 = this->p_OpMatrix->getMatrixElement(2,1);

	this->b11 = a11*a11 - a11*this->trace + this->determinant + a12*a21;
	this->b21 = a21*a11 + a21*a22 -a21*this->trace;
	this->b31 = a21*a32;
};

/*
 * Wilkinson位移QR迭代 隐式 初始化
 */
void DoubleShiftQRIteration::initForWilkinsonImplicitQR()
{
	//this->p_QQTSubMatrix_Implicit_Step->resetMatrixToI();
	//this->p_QSubMatrix_Implicit_Step->resetMatrixToI();
	//重置转换矩阵
	this->p_QQTMatrix_Implicit_Step->resetMatrixToI();
	//设置转换矩阵，此处应为3x3转换子矩阵
	this->p_QQTMatrix_Implicit_Step->resizeMatrix(3,3);

	//更新wilkinson偏移 trace以及determinant
	//generateWilkinsonShift();

	//hessenberg矩阵 生成双步移位后的矩阵第一列
	generateBMatrixFirstColumn();

	//双步移位后的矩阵首列三个元素构成三元向量
	p_TransVectorForQStep->resetDimension(3);
	p_TransVectorForQStep->setElement(0,this->b11);
	p_TransVectorForQStep->setElement(1,this->b21);
	p_TransVectorForQStep->setElement(2,this->b31);

	//计算Q1T
	this->m_HouseholderTrans.reload(this->p_TransVectorForQStep);
	//p_QTSubMatrix_Implicit_Step->copyMatrixElementNoCheck(this->p_HouseholderTrans->getHouseholderMatrixToE1(true));
	this->m_HouseholderTrans.getHouseholderMatrixToE1_ReverseElement(p_QQTMatrix_Implicit_Step);

	cout << "initForWilkinsonImplicitQR---Q1T sub" << endl;
	p_QQTMatrix_Implicit_Step->printMatrix();
	//将子矩阵QT 升级为全尺寸QT矩阵
	upgradeQQTSubMatrix(0);
	cout << "initForWilkinsonImplicitQR---Q1T" << endl;
	p_QQTMatrix_Implicit_Step->printMatrix();
	//更新总体转换矩阵QT total
	updateQT_Total();
	//更新hessenberg矩阵 By QT
	updateHessenbergOpMatrix_By_QT_IM_QRIteration();

	//生成Q1
	this->m_Transposer.transposeSquareMatrix(p_QQTMatrix_Implicit_Step);
	//this->p_Transposer->transposeMatrix();
	//this->p_QSubMatrix_Implicit_Step->copyMatrixElementNoCheck(this->p_Transposer->getTransposeMatrix());

	//cout << "initForWilkinsonImplicitQR---Q1 sub" << endl;
	//p_QQTMatrix_Implicit_Step->printMatrix();
	//将子矩阵Q 升级为全尺寸Q矩阵
	//upgradeQQTSubMatrix(0);
	cout << "initForWilkinsonImplicitQR---Q1" << endl;
	p_QQTMatrix_Implicit_Step->printMatrix();
	//更新总体转换矩阵Q total
	updateQ_Total();
	//更新hessenberg矩阵 By Q
	updateHessenbergOpMatrix_By_Q_IM_QRIteration();

	cout << "initForWilkinsonImplicitQR---Q1T * Hessenberg * Q1" << endl;
	this->p_OpMatrix->printMatrix();
};

/*
 * 将转换子矩阵Q/QT 根据当前迭代进度,升级为全尺寸Q/QT矩阵
 */
void DoubleShiftQRIteration::upgradeQQTSubMatrix(int iterateNum)
{
	this->p_QQTMatrix_Implicit_Step->resizeMatrix(this->p_OpMatrix->rowNum, this->p_OpMatrix->columnNum);

	//左上对角线3阶子矩阵下移
	p_QQTMatrix_Implicit_Step->moveDiagonalSubMatrixDown(0,2,iterateNum);

	/*
	for(int i=0,m=iterateNum; i<p_QSubMatrix_Implicit_Step->rowNum; i++,m++)
	{
		for(int j=0,n=iterateNum; j<p_QSubMatrix_Implicit_Step->columnNum; j++,n++)
		{
			p_QTMatrix_Implicit_Step->setMatrixElement(m,n,p_QTSubMatrix_Implicit_Step->getMatrixElement(i,j));
			p_QMatrix_Implicit_Step->setMatrixElement(m,n,p_QSubMatrix_Implicit_Step->getMatrixElement(i,j));
		}
	}
	*/
};

/*
 * 将Qn-1子矩阵,升级为全尺寸Q矩阵
 */
void DoubleShiftQRIteration::upgradeQQTLastSubMatrix()
{
	int lastQIndex = p_OpMatrix->rowNum-1;

	this->p_QQTMatrix_Implicit_Step->resizeMatrix(this->p_OpMatrix->rowNum, this->p_OpMatrix->columnNum);
	this->p_QQTMatrix_Implicit_Step->moveDiagonalSubMatrixDown(0,1,lastQIndex-1);
	//左上对角线2阶子矩阵下移至右下


	//this->p_QTMatrix_Implicit_Step->resetMatrixToI();
	//this->p_QMatrix_Implicit_Step->resetMatrixToI();
	/*
	for(int i=0,m=lastQIndex-1; i<p_QSubMatrix_Implicit_LastStep->rowNum; i++,m++)
	{
		for(int j=0,n=lastQIndex-1; j<p_QSubMatrix_Implicit_LastStep->columnNum; j++,n++)
		{
			p_QTMatrix_Implicit_Step->setMatrixElement(m,n,p_QTSubMatrix_Implicit_LastStep->getMatrixElement(i,j));
			p_QMatrix_Implicit_Step->setMatrixElement(m,n,p_QSubMatrix_Implicit_LastStep->getMatrixElement(i,j));
		}
	}
	*/
};


/*
 * 隐式QR迭代 更新hessenberg操作矩阵 计算QT * OP-Hessenberg
 */
void DoubleShiftQRIteration::updateHessenbergOpMatrix_By_QT_IM_QRIteration()
{
	//计算QT * OP-Hessenberg
	this->m_Multiplier.reload(p_QQTMatrix_Implicit_Step, p_OpMatrix, p_TempMatrix);
	this->m_Multiplier.multiplyCalc();
	p_OpMatrix->copyMatrixElementNoCheck(p_TempMatrix);
};

/*
 * 隐式QR迭代 更新hessenberg操作矩阵 计算 OP-Hessenberg * Q
 */
void DoubleShiftQRIteration::updateHessenbergOpMatrix_By_Q_IM_QRIteration()
{
	//计算QT * OP-Hessenberg * Q
	this->m_Multiplier.reload(p_OpMatrix, p_QQTMatrix_Implicit_Step, p_TempMatrix);
	this->m_Multiplier.multiplyCalc();
	p_OpMatrix->copyMatrixElementNoCheck(p_TempMatrix);
};

/*
 * Wilkinson双位移QR迭代 隐式(不需要先进行矩阵Hessenberg格式化，init和reload时已经完成转换)
 */
void DoubleShiftQRIteration::wilkinson_IM_QRIteration()
{
	//操作矩阵转换为Hessenberg矩阵
	//this->generateHessenbergOpMatrix();
	cout << "DoubleShiftQRIteration--initForImplicitQR----OP Hessenberg Matrix" << endl;
	this->p_OpMatrix->printMatrix();

	int iteratNumber = 10;
	int i = 0;
	while(i < iteratNumber)
	{
		wilkinson_IM_QRIteration_Step();
		i++;
	}
};

/*
 * Wilkinson双位移QR迭代 隐式 作为迭代接口单独使用(不需要先进行矩阵Hessenberg格式化，init和reload时已经完成转换)
 */
void DoubleShiftQRIteration::wilkinson_IM_QRIteration_Single()
{

	cout << "DoubleShiftQRIteration--initForImplicitQR----OP Hessenberg Matrix" << endl;
	this->p_OpMatrix->printMatrix();

	this->wilkinson_IM_QRIteration_Step();
};
/*
 * Wilkinson双位移QR迭代 隐式 -单步
 */
void DoubleShiftQRIteration::wilkinson_IM_QRIteration_Step()
{
	//this->p_QTMatrix_Implicit_Total->resetMatrixToI();
	//this->p_QMatrix_Implicit_Total->resetMatrixToI();

	initForWilkinsonImplicitQR();

	//chasing bulge
	int bulgeNumber = this->p_OpMatrix->columnNum-3;
	for(int j=0; j<bulgeNumber; j++)
	{
		//生成隐式Q及QT
		double a1,a2,a3;
		a1 = this->p_OpMatrix->getMatrixElement(j+1,j);
		a2 = this->p_OpMatrix->getMatrixElement(j+2,j);
		a3 = this->p_OpMatrix->getMatrixElement(j+3,j);

		//双步移位 三个元素构成三元向量
		this->p_TransVectorForQStep->resetDimension(3);
		this->p_TransVectorForQStep->setElement(0,a1);
		this->p_TransVectorForQStep->setElement(1,a2);
		this->p_TransVectorForQStep->setElement(2,a3);

		//重置转换矩阵
		this->p_QQTMatrix_Implicit_Step->resetMatrixToI();
		//设置转换矩阵，此处应为3x3转换子矩阵
		this->p_QQTMatrix_Implicit_Step->resizeMatrix(3,3);

		//计算Q1T
		this->m_HouseholderTrans.reload(this->p_TransVectorForQStep);
		//p_QTSubMatrix_Implicit_Step->copyMatrixElementNoCheck(this->p_HouseholderTrans->getHouseholderMatrixToE1(true));
		this->m_HouseholderTrans.getHouseholderMatrixToE1_ReverseElement(p_QQTMatrix_Implicit_Step);

		cout << "wilkinson_IM_QRIteration_Step---QiT Sub:" << j <<endl;
		p_QQTMatrix_Implicit_Step->printMatrix();

		//将子矩阵QiT 升级为全尺寸QiT矩阵
		upgradeQQTSubMatrix(j+1);

		cout << "wilkinson_IM_QRIteration_Step---QiT Full:" << j <<endl;
		p_QQTMatrix_Implicit_Step->printMatrix();

		//更新总体转换矩阵QT total
		updateQT_Total();
		//更新hessenberg矩阵 By QT
		updateHessenbergOpMatrix_By_QT_IM_QRIteration();

		//生成全尺寸Qi
		this->m_Transposer.transposeSquareMatrix(p_QQTMatrix_Implicit_Step);
		//this->p_Transposer->transposeMatrix();
		//this->p_QSubMatrix_Implicit_Step->copyMatrixElementNoCheck(this->p_Transposer->getTransposeMatrix());
		cout << "wilkinson_IM_QRIteration_Step---Qi Full:" << j <<endl;
		p_QQTMatrix_Implicit_Step->printMatrix();

		//更新总体转换矩阵Q total
		updateQ_Total();
		//更新hessenberg矩阵
		updateHessenbergOpMatrix_By_Q_IM_QRIteration();
		cout << "wilkinson_IM_QRIteration_Step---QiT * Hessenberg * Qi" << endl;
		p_OpMatrix->printMatrix();
	}

	endForWilkinsonImplicitQR();
};

/*
 * Wilkinson位移QR迭代 隐式 收尾最后一次迭代
 */
void DoubleShiftQRIteration::endForWilkinsonImplicitQR()
{
	int lastQIndex = p_OpMatrix->rowNum-1;
	//生成隐式Q及QT
	double an1,an2;
	an1 = this->p_OpMatrix->getMatrixElement(lastQIndex-1,lastQIndex-2);
	an2 = this->p_OpMatrix->getMatrixElement(lastQIndex,lastQIndex-2);

	//双步移位最后一步 2个元素构成2元向量
	this->p_TransVectorForQStep->resetDimension(2);
	this->p_TransVectorForQStep->setElement(0,an1);
	this->p_TransVectorForQStep->setElement(1,an2);

	//重置转换矩阵
	this->p_QQTMatrix_Implicit_Step->resetMatrixToI();
	//设置转换矩阵，此处应为2x2转换子矩阵
	this->p_QQTMatrix_Implicit_Step->resizeMatrix(2,2);

	this->m_GivensTrans.reload(p_TransVectorForQStep);
	//this->p_QTSubMatrix_Implicit_LastStep->copyMatrixElementNoCheck(p_GivensTrans->getGivensMatrixPreMultiple(1));
	this->m_GivensTrans.getGivensMatrixPreMultiple(1,p_QQTMatrix_Implicit_Step);

	cout << "endForWilkinsonImplicitQR---Qn_1T Sub:" <<endl;
	p_QQTMatrix_Implicit_Step->printMatrix();
	//将子矩阵升级为全尺寸转换矩阵
	upgradeQQTLastSubMatrix();
	cout << "endForWilkinsonImplicitQR---Qn_1T Full:" <<endl;
	p_QQTMatrix_Implicit_Step->printMatrix();

	//更新总体转换矩阵QT total
	updateQT_Total();

	//更新hessenberg矩阵 By QT
	updateHessenbergOpMatrix_By_QT_IM_QRIteration();

	this->m_Transposer.transposeSquareMatrix(p_QQTMatrix_Implicit_Step);
	//this->p_Transposer->transposeMatrix();
	//this->p_QSubMatrix_Implicit_LastStep->copyMatrixElementNoCheck(this->p_Transposer->getTransposeMatrix());


	cout << "endForWilkinsonImplicitQR---Qn_1 Full:" << endl;
	p_QQTMatrix_Implicit_Step->printMatrix();

	//将子矩阵升级为全尺寸转换矩阵
	//upgradeQQTLastSubMatrix();
	//cout << "endForWilkinsonImplicitQR---Qn_1T" <<endl;
	//p_QTMatrix_Implicit_Step->printMatrix();
	//cout << "endForWilkinsonImplicitQR---Qn_1" << endl;
	//p_QMatrix_Implicit_Step->printMatrix();

	//更新总体转换矩阵Q total
	updateQ_Total();

	//更新hessenberg矩阵 By Q
	updateHessenbergOpMatrix_By_Q_IM_QRIteration();

	cout << "endForWilkinsonImplicitQR---Qn_1T * Hessenberg * Qn_1" << endl;
	p_OpMatrix->printMatrix();
};

/*
 * 更新总体转换矩阵Q QT
 */
void DoubleShiftQRIteration::updateQT_Total()
{
	//计算QT_Step * QT_Total
	this->m_Multiplier.reload(p_QQTMatrix_Implicit_Step, p_QTMatrix_Implicit_Total, p_TempMatrix);
	this->m_Multiplier.multiplyCalc();
	this->p_QTMatrix_Implicit_Total->copyMatrixElementNoCheck(p_TempMatrix);

};

/*
 * 更新总体转换矩阵Q QT
 */
void DoubleShiftQRIteration::updateQ_Total()
{
	//计算Q_Total * Q_Step
	this->m_Multiplier.reload(p_QMatrix_Implicit_Total, p_QQTMatrix_Implicit_Step, p_TempMatrix);
	this->m_Multiplier.multiplyCalc();
	this->p_QMatrix_Implicit_Total->copyMatrixElementNoCheck(p_TempMatrix);
};

BasicMatrix* DoubleShiftQRIteration::getOpHessenbergMatrix()
{
	return this->p_OpMatrix;
};

//获取总体转换矩阵
BasicMatrix* DoubleShiftQRIteration::getQTMatrix_Total()
{
	return this->p_QTMatrix_Implicit_Total;
};

BasicMatrix* DoubleShiftQRIteration::getQMatrix_Total()
{
	return this->p_QMatrix_Implicit_Total;
};
