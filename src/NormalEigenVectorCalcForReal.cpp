/*
 * NormalEigenVectorCalcForReal.cpp
 *
 *  Created on: 2017��6��19��
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
 * ����ָ���Խ������������������Ǿ���Խ���ǰ��Ԫ�����Ƿ񻹴�����ͬȡֵ������ֵ��������ڣ��򷵻ض�ӦԪ�صĶԽ�������
 * 2017.6.20:��Ҫ����ǰ�����ڸ�������ֵ��������ܿ���������ֵ����
 */
bool NormalEigenVectorCalcForReal::findSameEigenValuePreIndex(int diagonalIndex)
{
	double diagonalElement = p_UpTriangleMatrix->getMatrixElement(diagonalIndex,diagonalIndex);
	double temp;
	bool result = false;
	//������ͬ��ǰ������ֵ
	for(int i=diagonalIndex-1; i>=0; i--)
	{
		if(0 != i)
		{
			//�жϴζԽ�Ԫ�Ƿ�Ϊ0�������Ϊ0��������ǰλ�ڸ�������ֵ�Խǿ飬��Ҫ�ܿ�
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

//����ָ���ĶԽ��������������Ӧ�Խ���Ԫ��(����ֵ)����������
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
//����ָ���ĶԽ��������������Ӧ�Խ���Ԫ��(����ֵ)����������
void NormalEigenVectorCalcForReal::getEigenVectorWithoutPreCheck(int diagonalIndex, BasicVector* p_resultVector)
{
	double diagonalElement,temp;
	diagonalElement = p_UpTriangleMatrix->getMatrixElement(diagonalIndex,diagonalIndex);

	p_resultVector->resetDimension(p_UpTriangleMatrix->columnNum);

	p_SubMatrix->resizeMatrix(diagonalIndex,diagonalIndex);
	p_SubInvserMatrix->resizeMatrix(diagonalIndex,diagonalIndex);
	p_SubInvserMatrix->resetMatrixToI();

	//����index��ʼ�������Ӿ���
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

	//�����Ӿ���������
	m_Inverser.reload(p_SubMatrix, p_SubInvserMatrix);
	m_Inverser.generateInverseMatrix();

	//��ȡdiagonalIndex��Ӧ��������
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
	//p_resultVector->normalizationVector();
};
