/*
 * GeneralizedEigenVectorCalcForReal.h
 *
 *  Created on: 2017年6月22日
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
	//根据指定的对角线索引，计算对应对角线元素(特征值)的特征向量，对前驱元素要进行排查，如果有重复特征值，需要考虑重新规划索引位置
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
