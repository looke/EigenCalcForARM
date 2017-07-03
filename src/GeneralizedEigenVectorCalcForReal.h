/*
 * GeneralizedEigenVectorCalcForReal.h
 *
 *  Created on: 2017��6��22��
 *      Author: looke
 */

#ifndef GENERALIZEDEIGENVECTORCALCFORREAL_H_
#define GENERALIZEDEIGENVECTORCALCFORREAL_H_

#include "BasicMatrix.h"
#include "BasicVector.h"
#include "MatrixMultiplier.h"
#include "MatrixInverser.h"
#include "NormalEigenVectorCalcForReal.h"
class GeneralizedEigenVectorCalcForReal
{
public:
	GeneralizedEigenVectorCalcForReal(BasicMatrix* p_input_UpTriangleMatrix_A, BasicMatrix* p_input_UpTriangleMatrix_B, BasicMatrix* p_input_UpTriangleMatrix_BinvA, BasicMatrix* p_input_SubMatrix, BasicMatrix* p_input_SubInverMatrix, BasicMatrix* p_input_TempMatrix);

	void init(BasicMatrix* p_input_UpTriangleMatrix_A, BasicMatrix* p_input_UpTriangleMatrix_B, BasicMatrix* p_input_UpTriangleMatrix_BinvA, BasicMatrix* p_input_SubMatrix, BasicMatrix* p_input_SubInverMatrix, BasicMatrix* p_input_TempMatrix);
	void reload(BasicMatrix* p_input_UpTriangleMatrix_A, BasicMatrix* p_input_UpTriangleMatrix_B, BasicMatrix* p_input_UpTriangleMatrix_BinvA, BasicMatrix* p_input_SubMatrix, BasicMatrix* p_input_SubInverMatrix, BasicMatrix* p_input_TempMatrix);

	void generateDeflationEnd();

	void generateDeflationMatrix_B(BasicMatrix* deflated_Matrix_B);
	void generateDeflationMatrix_A(BasicMatrix* deflated_Matrix_A);
	void generateDeflationMatrix_Temp(BasicMatrix* deflated_Matrix_Temp);
	//����ָ���ĶԽ��������������Ӧ�Խ���Ԫ��(����ֵ)��������������ǰ��Ԫ��Ҫ�����Ų飬������ظ�����ֵ����Ҫ�������¹滮����λ��
	bool getEigenVector(int diagonalIndex, BasicVector* p_resultVector);

protected:
	int preIndex;
	int deflationStart;
	int deflationEnd;
	BasicMatrix* p_UpTriangleMatrix_A;
	BasicMatrix* p_UpTriangleMatrix_B;

	BasicMatrix* p_UpTriangleMatrix_BinvA;
	BasicMatrix* p_SubMatrix;
	BasicMatrix* p_SubInvserMatrix;
	BasicMatrix* p_TempMatrix;

	MatrixMultiplier m_Multiplier;
	MatrixInverser m_Inverser;
	NormalEigenVectorCalcForReal m_EigenVectorCalcForReal;

};

#endif /* GENERALIZEDEIGENVECTORCALCFORREAL_H_ */
