/*
 * NormalEigenVectorCalcForReal.cpp
 *
 *  Created on: 2017年6月19日
 *      Author: looke
 */

#include "NormalEigenVectorCalcForReal.h"

NormalEigenVectorCalcForReal::NormalEigenVectorCalcForReal(BasicMatrix* p_input_UpTriangleMatrix, BasicMatrix* p_input_SubMatrix, BasicMatrix* p_input_SubInverMatrix, BasicMatrix* p_input_TempMatrix)
:m_Multiplier(p_input_UpTriangleMatrix,p_input_SubMatrix,p_input_TempMatrix),m_Inverser(p_input_SubMatrix,p_input_SubInverMatrix)
{
	this->init(p_input_UpTriangleMatrix, p_input_SubMatrix, p_input_SubInverMatrix, p_input_TempMatrix);
};

void NormalEigenVectorCalcForReal::init(BasicMatrix* p_input_UpTriangleMatrix, BasicMatrix* p_input_SubMatrix, BasicMatrix* p_input_SubInverMatrix, BasicMatrix* p_input_TempMatrix)
{
	this->p_UpTriangleMatrix = p_input_UpTriangleMatrix;
	this->p_SubMatrix = p_input_SubMatrix;
	this->p_SubInvserMatrix = p_input_SubInverMatrix;
	this->p_TempMatrix = p_input_TempMatrix;

	this->preIndex = 0;
};


void NormalEigenVectorCalcForReal::reload(BasicMatrix* p_input_UpTriangleMatrix, BasicMatrix* p_input_SubMatrix, BasicMatrix* p_input_SubInverMatrix, BasicMatrix* p_input_TempMatrix)
{
	this->init(p_input_UpTriangleMatrix, p_input_SubMatrix, p_input_SubInverMatrix, p_input_TempMatrix);
};

/*
 * 根据指定对角线索引，查找上三角矩阵对角线前驱元素上是否还存在相同取值的特征值，如果存在，则返回对应元素的对角线索引
 * 2017.6.20:需要考虑前驱存在复数特征值的情况，避开复数特征值区域
 */
bool NormalEigenVectorCalcForReal::findSameEigenValuePreIndex(int diagonalIndex)
{
	double diagonalElement = p_UpTriangleMatrix->getMatrixElement(diagonalIndex,diagonalIndex);
	double temp;
	bool result = false;
	//查找相同的前驱特征值
	for(int i=diagonalIndex-1; i>=0; i--)
	{
		if(0 != i)
		{
			//判断次对角元是否为0，如果不为0，表明当前位于复数特征值对角块，需要避开
			temp = p_UpTriangleMatrix->getMatrixElement(i,i-1);
			if(0 != temp)
			{
				i--;
				continue;
			}
		}
		temp = p_UpTriangleMatrix->getMatrixElement(i,i);

		if(diagonalElement == temp)
		{
			this->preIndex = i;
			result = true;
		}
	}
	return result;
};

//根据指定的对角线索引，计算对应对角线元素(特征值)的特征向量
bool NormalEigenVectorCalcForReal::getEigenVector(int diagonalIndex, BasicVector* p_resultVector)
{
	if(diagonalIndex < 0 || diagonalIndex >= p_UpTriangleMatrix->columnNum)
	{
		return false;
	}
	int targetIndex = diagonalIndex;
	bool hasPreIndex = findSameEigenValuePreIndex(diagonalIndex);

	if(hasPreIndex)
	{
		targetIndex = this->preIndex;
	}

	getEigenVectorWithoutPreCheck(targetIndex, p_resultVector);

	return true;
}
//根据指定的对角线索引，计算对应对角线元素(特征值)的特征向量
void NormalEigenVectorCalcForReal::getEigenVectorWithoutPreCheck(int diagonalIndex, BasicVector* p_resultVector)
{
	double diagonalElement,temp;
	diagonalElement = p_UpTriangleMatrix->getMatrixElement(diagonalIndex,diagonalIndex);

	p_resultVector->resetDimension(p_UpTriangleMatrix->columnNum);

	p_SubMatrix->resizeMatrix(diagonalIndex,diagonalIndex);
	p_SubInvserMatrix->resizeMatrix(diagonalIndex,diagonalIndex);
	p_SubInvserMatrix->resetMatrixToI();

	//根据index初始化左上子矩阵
	for(int i=0; i < diagonalIndex; i++)
	{
		for(int j=0; j < diagonalIndex; j++)
		{
			temp = p_UpTriangleMatrix->getMatrixElement(i,j);
			if(i!= j)
			{
				p_SubMatrix->setMatrixElement(i,j,temp);
			}
			else
			{
				p_SubMatrix->setMatrixElement(i,j,temp-diagonalElement);
			}

		}
	}

	//计算子矩阵的逆矩阵
	m_Inverser.reload(p_SubMatrix, p_SubInvserMatrix);
	m_Inverser.generateInverseMatrix();

	//获取diagonalIndex对应的子向量
	p_SubMatrix->resizeMatrix(diagonalIndex,1);
	for(int i=0; i < diagonalIndex; i++)
	{
		p_SubMatrix->setMatrixElement(i,0,p_UpTriangleMatrix->getMatrixElement(i,diagonalIndex));
	}

	p_TempMatrix->resizeMatrix(diagonalIndex,1);

	m_Multiplier.reload(p_SubInvserMatrix, p_SubMatrix, p_TempMatrix);
	m_Multiplier.multiplyCalc();


	for(int i=0; i < diagonalIndex; i++)
	{
		temp = 0 - p_TempMatrix->getMatrixElement(i,0);
		p_resultVector->setElement(i,temp);
	}

	p_resultVector->setElement(diagonalIndex,1);

	for(int i=diagonalIndex+1; i < p_UpTriangleMatrix->columnNum; i++)
	{
		p_resultVector->setElement(i,0);
	}
	p_resultVector->normalizationVector();
};
