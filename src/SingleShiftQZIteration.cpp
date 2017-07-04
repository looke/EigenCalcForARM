/*
 * SingleShiftQZIteration.cpp
 *
 *  Created on: 2017年5月27日
 *      Author: looke
 */

#include "SingleShiftQZIteration.h"
//#include <iostream>
using namespace std;

//SingleShiftQZIteration::SingleShiftQZIteration()
//{};

SingleShiftQZIteration::SingleShiftQZIteration(BasicMatrix* p_input_OpMatrix_A, BasicMatrix* p_input_OpMatrix_B,BasicMatrix* p_input_QMatrix_Total,BasicMatrix* p_input_ZMatrix_Total,BasicMatrix* p_input_QZMatrix_Step,BasicMatrix* p_input_TempMatrix_Trans,BasicMatrix* p_input_TempMatrix)
:m_ABInvCalc(),m_GivensTrans(p_input_OpMatrix_A->getColumnVector(0)),m_Multiplier(p_input_OpMatrix_A,p_input_OpMatrix_B,p_input_TempMatrix),m_HessenbergTriangleFormular(p_input_OpMatrix_A,p_input_OpMatrix_B,p_input_QMatrix_Total,p_input_ZMatrix_Total,p_input_QZMatrix_Step,p_input_TempMatrix_Trans,p_input_TempMatrix)
{
	this->init(p_input_OpMatrix_A,p_input_OpMatrix_B,p_input_QMatrix_Total,p_input_ZMatrix_Total,p_input_QZMatrix_Step,p_input_TempMatrix_Trans,p_input_TempMatrix);
};

void SingleShiftQZIteration::init(BasicMatrix* p_input_OpMatrix_A, BasicMatrix* p_input_OpMatrix_B,BasicMatrix* p_input_QMatrix_Total,BasicMatrix* p_input_ZMatrix_Total,BasicMatrix* p_input_QZMatrix_Step,BasicMatrix* p_input_TempMatrix_Trans,BasicMatrix* p_input_TempMatrix)
{
	//原始操作矩阵A Hessenberg
	this->p_OpMatrix_A = p_input_OpMatrix_A;
	//原始操作矩阵B Triangle
	this->p_OpMatrix_B = p_input_OpMatrix_B;

	//Q 矩阵 隐式迭代 综合 Z用于右乘OP矩阵
	this->p_QMatrix_Implicit_Total = p_input_QMatrix_Total;
	//Z 矩阵 隐式迭代 综合 Q用于左乘OP矩阵
	this->p_ZMatrix_Implicit_Total = p_input_ZMatrix_Total;

	//Q/Z 矩阵 隐式迭代 分步 Q用于左乘OP矩阵 Z用于右乘OP矩阵
	this->p_QZMatrix_Implicit_Step = p_input_QZMatrix_Step;

	//中间过程矩阵
	this->p_TempMatrix_Trans = p_input_TempMatrix_Trans;
	//中间过程矩阵
	this->p_TempMatrix = p_input_TempMatrix;

	//this->generateHessenTriangleOpMatrix();
};

void SingleShiftQZIteration::reload(BasicMatrix* p_intput_OpMatrix_A, BasicMatrix* p_intput_OpMatrix_B,BasicMatrix* p_input_QMatrix_Total,BasicMatrix* p_input_ZMatrix_Total,BasicMatrix* p_input_QZMatrix_Step,BasicMatrix* p_input_TempMatrix_Trans,BasicMatrix* p_input_TempMatrix)
{
	this->init(p_intput_OpMatrix_A,p_intput_OpMatrix_B,p_input_QMatrix_Total,p_input_ZMatrix_Total,p_input_QZMatrix_Step,p_input_TempMatrix_Trans,p_input_TempMatrix);
};


/*
 * 生成H-T操作矩阵对
 */
void SingleShiftQZIteration::generateHessenTriangleOpMatrix()
{
	this->m_HessenbergTriangleFormular.reload(p_OpMatrix_A, p_OpMatrix_B,p_QMatrix_Implicit_Total,p_ZMatrix_Implicit_Total,p_QZMatrix_Implicit_Step,p_TempMatrix_Trans,p_TempMatrix);
	this->m_HessenbergTriangleFormular.formularABMatrix();
	//this->p_OpMatrix_Hessenberg->copyMatrixElementNoCheck(p_HessenbergTriangleFormular->getHessenbergMatrixA());
	//this->p_OpMatrix_Triangle->copyMatrixElementNoCheck(p_HessenbergTriangleFormular->getTriangleMatrixB());

	//cout << "SingleShiftQZIteration--generateHessenTriangleOpMatrix----OP Hessenberg Matrix" << endl;
	//this->p_OpMatrix_A->printMatrix();
	//cout << "SingleShiftQZIteration--generateHessenTriangleOpMatrix----OP Triangle Matrix" << endl;
	//this->p_OpMatrix_B->printMatrix();
};

/*
 * 使用Q矩阵左乘H-T矩阵对
 */
void SingleShiftQZIteration::updateHTMatrixByQ()
{
	//cout << "SingleShiftQZIteration--updateHTMatrixByQ----OP Hessenberg Matrix Before" << endl;
	//this->p_OpMatrix_A->printMatrix();
	//cout << "SingleShiftQZIteration--updateHTMatrixByQ----OP Triangle Matrix Before" << endl;
	//this->p_OpMatrix_B->printMatrix();

	//更新矩阵Hessenberg A
	this->m_Multiplier.reload(p_QZMatrix_Implicit_Step, p_OpMatrix_A, p_TempMatrix);
	this->m_Multiplier.multiplyCalc();
	p_OpMatrix_A->copyMatrixElementNoCheck(p_TempMatrix);

	//更新矩阵Triangle B
	this->m_Multiplier.reload(p_QZMatrix_Implicit_Step, p_OpMatrix_B,p_TempMatrix);
	this->m_Multiplier.multiplyCalc();
	p_OpMatrix_B->copyMatrixElementNoCheck(p_TempMatrix);

	//cout << "SingleShiftQZIteration--updateHTMatrixByQ----OP Hessenberg Matrix After" << endl;
	//this->p_OpMatrix_A->printMatrix();
	//cout << "SingleShiftQZIteration--updateHTMatrixByQ----OP Triangle Matrix After" << endl;
	//this->p_OpMatrix_B->printMatrix();
};

/*
 * 使用Z矩阵右乘H-T矩阵对
 */
void SingleShiftQZIteration::updateHTMatrixByZ()
{
	//cout << "SingleShiftQZIteration--updateHTMatrixByZ----OP Hessenberg Matrix Before" << endl;
	//this->p_OpMatrix_A->printMatrix();
	//cout << "SingleShiftQZIteration--updateHTMatrixByZ----OP Triangle Matrix Before" << endl;
	//this->p_OpMatrix_B->printMatrix();

	//更新矩阵Hessenberg A
	this->m_Multiplier.reload(p_OpMatrix_A, p_QZMatrix_Implicit_Step, p_TempMatrix);
	this->m_Multiplier.multiplyCalc();
	p_OpMatrix_A->copyMatrixElementNoCheck(p_TempMatrix);

	//更新矩阵Triangle B
	this->m_Multiplier.reload(p_OpMatrix_B, p_QZMatrix_Implicit_Step, p_TempMatrix);
	this->m_Multiplier.multiplyCalc();
	p_OpMatrix_B->copyMatrixElementNoCheck(p_TempMatrix);

	//cout << "SingleShiftQZIteration--updateHTMatrixByZ----OP Hessenberg Matrix After" << endl;
	//this->p_OpMatrix_A->printMatrix();
	//cout << "SingleShiftQZIteration--updateHTMatrixByZ----OP Triangle Matrix After" << endl;
	//this->p_OpMatrix_B->printMatrix();
};

/*
 * 更新综合转换矩阵Q/Z Total(将Step合并入Total)
 */
void SingleShiftQZIteration::updateQMatrix_Total()
{
	//更新矩阵Q total
	this->m_Multiplier.reload(p_QZMatrix_Implicit_Step, p_QMatrix_Implicit_Total, p_TempMatrix);
	this->m_Multiplier.multiplyCalc();
	p_QMatrix_Implicit_Total->copyMatrixElementNoCheck(p_TempMatrix);

};

/*
 * 更新综合转换矩阵Q/Z Total(将Step合并入Total)
 */
void SingleShiftQZIteration::updateZMatrix_Total()
{
	//更新矩阵Z total
	this->m_Multiplier.reload(p_ZMatrix_Implicit_Total, p_QZMatrix_Implicit_Step, p_TempMatrix);
	this->m_Multiplier.multiplyCalc();
	p_ZMatrix_Implicit_Total->copyMatrixElementNoCheck(p_TempMatrix);
};

/*
 * 初始化隐式迭代
 * 执行隐式迭代第一步，在H-T矩阵对的次对角线上形成一个突出部位
 */
void SingleShiftQZIteration::initForImplicitQZ(double input_ShiftValue)
{
	//cout << "SingleShiftQZIteration--initForImplicitQZ" << endl;
	//p_QMatrix_Implicit_Step->resetMatrixToI();
	p_QZMatrix_Implicit_Step->resetMatrixToI();

	//QZ单步移位 A*B^-1 第一列
	double temp_old_a00 = this->p_OpMatrix_A->getMatrixElement(0,0);
	double temp_old_a10 = this->p_OpMatrix_A->getMatrixElement(1,0);

	double temp_old_b00 = this->p_OpMatrix_B->getMatrixElement(0,0);


	double temp_new_a00 = temp_old_a00/temp_old_b00 - input_ShiftValue;
	double temp_new_a10 = temp_old_a10/temp_old_b00;

	//更新第一列
	this->p_OpMatrix_A->setMatrixElement(0,0,temp_new_a00);
	this->p_OpMatrix_A->setMatrixElement(1,0,temp_new_a10);

	//生成隐式Q
	BasicVector* p_firstColumnVector = this->p_OpMatrix_A->getColumnVector(0);
	m_GivensTrans.reload(p_firstColumnVector);
	m_GivensTrans.getGivensMatrixPreMultiple(1,p_QZMatrix_Implicit_Step);

	//cout << "SingleShiftQZIteration--initForImplicitQZ-----Q-Step" << endl;
	//p_QZMatrix_Implicit_Step->printMatrix();
	//还原Hessenberg矩阵
	this->p_OpMatrix_A->setMatrixElement(0,0,temp_old_a00);
	this->p_OpMatrix_A->setMatrixElement(1,0,temp_old_a10);

	//更新H-T矩阵对
	this->updateHTMatrixByQ();
	//更新综合转换矩阵Q Total
	updateQMatrix_Total();

	//生成隐式Z
	BasicVector* p_secondRowVector = this->p_OpMatrix_B->getRowVector(1);
	m_GivensTrans.reload(p_secondRowVector);
	m_GivensTrans.getGivensMatrixAfterMultiple(0,p_QZMatrix_Implicit_Step);

	//cout << "SingleShiftQZIteration--initForImplicitQZ-----Z-Step" << endl;
	//p_QZMatrix_Implicit_Step->printMatrix();
	//更新H-T矩阵对
	this->updateHTMatrixByZ();
	//更新综合转换矩阵Z Total
	updateZMatrix_Total();

	//cout << "SingleShiftQZIteration--initForImplicitQZ-----END" << endl;

};

//单值位移QZ迭代 隐式单步
void SingleShiftQZIteration::implicit_QZIteration_Step(double input_shiftValue)
{
	//p_QMatrix_Implicit_Total->resetMatrixToI();
	//p_ZMatrix_Implicit_Total->resetMatrixToI();

	this->initForImplicitQZ(input_shiftValue);

	int bulgeNumber = this->p_OpMatrix_A->columnNum-2;

	for(int j=0; j<bulgeNumber; j++)
	{
		//生成隐式Q
		BasicVector* p_columnVector = this->p_OpMatrix_A->getColumnVector(j);
		m_GivensTrans.reload(p_columnVector);
		//this->p_QMatrix_Implicit_Step->copyMatrixElementNoCheck(p_GivensTrans->getGivensMatrixPreMultiple(j+2));
		m_GivensTrans.getGivensMatrixPreMultiple(j+2,p_QZMatrix_Implicit_Step);
		//更新H-T矩阵对
		this->updateHTMatrixByQ();
		//更新综合转换矩阵Q Total
		updateQMatrix_Total();

		//生成隐式Z
		BasicVector* p_secondRowVector = this->p_OpMatrix_B->getRowVector(j+2);
		m_GivensTrans.reload(p_secondRowVector);
		//this->p_ZMatrix_Implicit_Step->copyMatrixElementNoCheck(p_GivensTrans->getGivensMatrixAfterMultiple(j+1));
		m_GivensTrans.getGivensMatrixAfterMultiple(j+1,p_QZMatrix_Implicit_Step);
		//更新H-T矩阵对
		this->updateHTMatrixByZ();
		//更新综合转换矩阵Q\Z Total
		updateZMatrix_Total();
	}

};

//单值位移QZ迭代 隐式
void SingleShiftQZIteration::implicit_QZIteration(double input_shiftValue)
{
	//generateHessenTriangleOpMatrix();
	int i=0;
	while (i<10)
	{
		implicit_QZIteration_Step(input_shiftValue);
		//cout << "SingleShiftQZIteration--implicit_QZIteration----OP Hessenberg Matrix after:" << i << "iteration" <<endl;
		//p_OpMatrix_A->printMatrix();
		i++;
	}
};
//单值rayleigh商位移QZ迭代 隐式
void SingleShiftQZIteration::rayleigh_Quotient_IM_QZIteration(int iterateNum)
{
	p_QMatrix_Implicit_Total->resetMatrixToI();
	p_ZMatrix_Implicit_Total->resetMatrixToI();

	this->generateHessenTriangleOpMatrix();

	//int rayleighValueIndex = p_OpMatrix_A->rowNum - 1;
	//double rayleighValue;
	int i=0;
	while (i<iterateNum)
	{
		this->rayleigh_Quotient_IM_QZIteration_Step();
		//rayleighValue = p_OpMatrix_A->getMatrixElement(rayleighValueIndex,rayleighValueIndex);
		//implicit_QZIteration_Step(rayleighValue);
		//cout << "SingleShiftQZIteration--rayleigh_Quotient_IM_QRIteration----OP Hessenberg Matrix after:" << i << "iteration" <<endl;
		//p_OpMatrix_Hessenberg->printMatrix();
		i++;
	}
};

//单值rayleigh商位移QZ迭代 隐式 单步
void SingleShiftQZIteration::rayleigh_Quotient_IM_QZIteration_Step()
{
	//计算A*B^-1 右下角元素 作为单步位移值
	m_ABInvCalc.generateABinvLastOne(p_OpMatrix_A, p_OpMatrix_B);
	double rayleighValue = m_ABInvCalc.getABinv_N_N();

	implicit_QZIteration_Step(rayleighValue);
};

//获取Q 总体转换矩阵
BasicMatrix* SingleShiftQZIteration::getQMatrix_Total()
{
	return this->p_QMatrix_Implicit_Total;
};

//获取Z 总体转换矩阵
BasicMatrix* SingleShiftQZIteration::getZMatrix_Total()
{
	return this->p_ZMatrix_Implicit_Total;
};
