/*
 * GeneraliedEigenVectorCalcForReal.cpp
 *
 *  Created on: 2017年6月22日
 *      Author: looke
 */

#include "GeneralizedEigenVectorCalcForReal.h"

GeneralizedEigenVectorCalcForReal::GeneralizedEigenVectorCalcForReal(BasicMatrix* p_input_UpTriangleMatrix_A, BasicMatrix* p_input_UpTriangleMatrix_B, BasicMatrix* p_input_UpTriangleMatrix_BinvA, BasicMatrix* p_input_SubMatrix, BasicMatrix* p_input_SubInverMatrix, BasicMatrix* p_input_TempMatrix)
:m_Multiplier(p_input_UpTriangleMatrix_A,p_input_UpTriangleMatrix_B,p_input_TempMatrix),
 m_Inverser(p_input_SubMatrix,p_input_SubInverMatrix),
 m_EigenVectorCalcForReal(p_input_UpTriangleMatrix_BinvA,p_input_SubMatrix,p_input_SubInverMatrix,p_input_TempMatrix)
{
	this->init(
			p_input_UpTriangleMatrix_A,
			p_input_UpTriangleMatrix_B,
			p_input_UpTriangleMatrix_BinvA,
			p_input_SubMatrix,
			p_input_SubInverMatrix,
			p_input_TempMatrix);
};
void GeneralizedEigenVectorCalcForReal::init(BasicMatrix* p_input_UpTriangleMatrix_A, BasicMatrix* p_input_UpTriangleMatrix_B, BasicMatrix* p_input_UpTriangleMatrix_BinvA, BasicMatrix* p_input_SubMatrix, BasicMatrix* p_input_SubInverMatrix, BasicMatrix* p_input_TempMatrix)
{
	this->p_UpTriangleMatrix_A = p_input_UpTriangleMatrix_A;
	this->p_UpTriangleMatrix_B = p_input_UpTriangleMatrix_B;
	this->p_UpTriangleMatrix_BinvA = p_input_UpTriangleMatrix_BinvA;
	this->p_SubMatrix = p_input_SubMatrix;
	this->p_SubInvserMatrix = p_input_SubInverMatrix;
	this->p_TempMatrix = p_input_TempMatrix;

	this->preIndex = 0;
	this->deflationStart = 0;
	this->deflationEnd = 0;
};


void GeneralizedEigenVectorCalcForReal::reload(BasicMatrix* p_input_UpTriangleMatrix_A, BasicMatrix* p_input_UpTriangleMatrix_B, BasicMatrix* p_input_UpTriangleMatrix_BinvA, BasicMatrix* p_input_SubMatrix, BasicMatrix* p_input_SubInverMatrix, BasicMatrix* p_input_TempMatrix)
{
	this->init(p_input_UpTriangleMatrix_A, p_input_UpTriangleMatrix_B, p_input_UpTriangleMatrix_BinvA, p_input_SubMatrix, p_input_SubInverMatrix, p_input_TempMatrix);
};


//根据指定的对角线索引，计算对应对角线元素(特征值)的特征向量，对前驱元素要进行排查，如果有重复特征值，需要考虑重新规划索引位置
bool GeneralizedEigenVectorCalcForReal::getEigenVector(int diagonalIndex, BasicVector* p_resultVector)
{
	bool result = true;
	p_resultVector->resetVectorElementToZero();

	//首先根据B矩阵主对角线进行降阶，计算B矩阵的降阶点
	this->generateDeflationEnd();

	//生成已降阶的B矩阵
	this->generateDeflationMatrix_B(p_SubMatrix);

	//首先计算B的逆矩阵
	this->generateDeflationMatrix_Temp(p_SubInvserMatrix);
	m_Inverser.reload(p_SubMatrix, p_SubInvserMatrix);
	m_Inverser.generateInverseMatrix();

	this->generateDeflationMatrix_A(p_SubMatrix);

	//计算Binv * A
	this->generateDeflationMatrix_Temp(p_UpTriangleMatrix_BinvA);
	m_Multiplier.reload(p_SubInvserMatrix, p_SubMatrix, p_UpTriangleMatrix_BinvA);
	m_Multiplier.multiplyCalc();

	//以Binv * A 计算特征向量
	m_EigenVectorCalcForReal.reload(p_UpTriangleMatrix_BinvA,p_SubMatrix,p_SubInvserMatrix,p_TempMatrix);
	result = m_EigenVectorCalcForReal.getEigenVector(diagonalIndex,p_resultVector);

	p_resultVector->resetDimension(this->p_UpTriangleMatrix_A->columnNum);
	return result;
};

/*
 * 根据B矩阵的主对角线0元，计算新的降阶点
 */
void GeneralizedEigenVectorCalcForReal::generateDeflationEnd()
{
	this->deflationEnd = p_UpTriangleMatrix_B->maxIndexOfZeroOnDiagonal() - 1;

	//如果不存在0元，则将deflationEnd设置到右下角元素索引
	if(-1 == deflationEnd)
	{
		deflationEnd = this->p_UpTriangleMatrix_A->columnNum - 1;
	}
};

void GeneralizedEigenVectorCalcForReal::generateDeflationMatrix_B(BasicMatrix* deflated_Matrix_B)
{
	int matrixSize = this->deflationEnd - this->deflationStart + 1;

	deflated_Matrix_B->resetMatrixToI();
	deflated_Matrix_B->resizeMatrix(matrixSize,matrixSize);
	deflated_Matrix_B->resetMatrixToI();

	for(int i=deflationStart, m=0; i<=deflationEnd; i++,m++)
	{
		for(int j=deflationStart, n=0; j<=deflationEnd; j++,n++)
		{
			deflated_Matrix_B->setMatrixElement(m,n,p_UpTriangleMatrix_B->getMatrixElement(i,j));
		}
	}
};

void GeneralizedEigenVectorCalcForReal::generateDeflationMatrix_A(BasicMatrix* deflated_Matrix_A)
{
	int matrixSize = this->deflationEnd - this->deflationStart + 1;

	deflated_Matrix_A->resetMatrixToI();
	deflated_Matrix_A->resizeMatrix(matrixSize,matrixSize);
	deflated_Matrix_A->resetMatrixToI();

	for(int i=deflationStart, m=0; i<=deflationEnd; i++,m++)
	{
		for(int j=deflationStart, n=0; j<=deflationEnd; j++,n++)
		{
			deflated_Matrix_A->setMatrixElement(m,n,p_UpTriangleMatrix_A->getMatrixElement(i,j));
		}
	}
};

void GeneralizedEigenVectorCalcForReal::generateDeflationMatrix_Temp(BasicMatrix* deflated_Matrix_Temp)
{
	int matrixSize = this->deflationEnd - this->deflationStart + 1;

	deflated_Matrix_Temp->resetMatrixToI();
	deflated_Matrix_Temp->resizeMatrix(matrixSize,matrixSize);
	deflated_Matrix_Temp->resetMatrixToI();
};
