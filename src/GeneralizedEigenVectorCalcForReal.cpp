/*
 * GeneraliedEigenVectorCalcForReal.cpp
 *
 *  Created on: 2017��6��22��
 *      Author: looke
 */

#include "GeneralizedEigenVectorCalcForReal.h"

GeneralizedEigenVectorCalcForReal::GeneralizedEigenVectorCalcForReal(BasicMatrix* p_input_UpTriangleMatrix_A, BasicMatrix* p_input_UpTriangleMatrix_B, BasicMatrix* p_input_UpTriangleMatrix_BinvA, BasicMatrix* p_input_SubMatrix, BasicMatrix* p_input_SubInverMatrix, BasicMatrix* p_input_TempMatrix)
:m_Multiplier(p_input_UpTriangleMatrix_A,p_input_UpTriangleMatrix_B,p_input_TempMatrix),m_Inverser(p_input_SubMatrix,p_input_SubInverMatrix),m_EigenVectorCalcForReal(p_input_UpTriangleMatrix_BinvA,p_input_SubMatrix,p_input_SubInverMatrix,p_input_TempMatrix)
{
	this->init(p_input_UpTriangleMatrix_A, p_input_UpTriangleMatrix_B, p_input_UpTriangleMatrix_BinvA, p_input_SubMatrix, p_input_SubInverMatrix, p_input_TempMatrix);
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


//����ָ���ĶԽ��������������Ӧ�Խ���Ԫ��(����ֵ)��������������ǰ��Ԫ��Ҫ�����Ų�
bool GeneralizedEigenVectorCalcForReal::getEigenVector(int diagonalIndex, BasicVector* p_resultVector)
{
	bool result = true;
	p_resultVector->resetVectorElementToZero();

	//���ȸ���B�������Խ��߽��н��ף�����B����Ľ��׵�
	this->generateDeflationEnd();

	//�����ѽ��׵�B����
	this->generateDeflationMatrix_B(p_SubMatrix);

	//���ȼ���B�������
	this->generateDeflationMatrix_Temp(p_SubInvserMatrix);
	m_Inverser.reload(p_SubMatrix, p_SubInvserMatrix);
	m_Inverser.generateInverseMatrix();

	this->generateDeflationMatrix_A(p_SubMatrix);

	//����Binv * A
	this->generateDeflationMatrix_Temp(p_UpTriangleMatrix_BinvA);
	m_Multiplier.reload(p_SubInvserMatrix, p_SubMatrix, p_UpTriangleMatrix_BinvA);
	m_Multiplier.multiplyCalc();

	//��Binv * A ������������
	m_EigenVectorCalcForReal.reload(p_UpTriangleMatrix_BinvA,p_SubMatrix,p_SubInvserMatrix,p_TempMatrix);
	result = m_EigenVectorCalcForReal.getEigenVector(diagonalIndex,p_resultVector);

	p_resultVector->resetDimension(this->p_UpTriangleMatrix_A->columnNum);
	return result;
};

/*
 * ����B��������Խ���0Ԫ�������µĽ��׵�
 */
void GeneralizedEigenVectorCalcForReal::generateDeflationEnd()
{
	this->deflationEnd = p_UpTriangleMatrix_B->maxIndexOfZeroOnDiagonal() - 1;

	//���������0Ԫ����deflationEnd���õ����½�Ԫ������
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
};
