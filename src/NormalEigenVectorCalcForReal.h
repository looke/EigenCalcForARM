/*
 * NormalEigenVectorCalcForReal.h
 *
 *  Created on: 2017年6月19日
 *      Author: looke
 */

#ifndef NORMALEIGENVECTORCALCFORREAL_H_
#define NORMALEIGENVECTORCALCFORREAL_H_


#include "BasicMatrix.h"
#include "BasicVector.h"
#include "MatrixMultiplier.h"
#include "MatrixInverser.h"

class NormalEigenVectorCalcForReal
{
public:
	NormalEigenVectorCalcForReal(BasicMatrix* p_input_UpTriangleMatrix, BasicMatrix* p_input_SubMatrix, BasicMatrix* p_input_SubInverMatrix, BasicMatrix* p_input_TempMatrix);

	void init(BasicMatrix* p_input_UpTriangleMatrix, BasicMatrix* p_input_SubMatrix, BasicMatrix* p_input_SubInverMatrix, BasicMatrix* p_input_TempMatrix);
	void reload(BasicMatrix* p_input_UpTriangleMatrix, BasicMatrix* p_input_SubMatrix, BasicMatrix* p_input_SubInverMatrix, BasicMatrix* p_input_TempMatrix);

	//根据指定的对角线索引，计算对应对角线元素(特征值)的特征向量，对前驱元素要进行排查
	bool getEigenVector(int diagonalIndex, BasicVector* p_resultVector);

	//根据指定对角线索引，查找上三角矩阵对角线前驱元素上是否还存在相同取值的特征值，如果存在，则返回对应元素的对角线索引
	bool findSameEigenValuePreIndex(int diagonalIndex);
protected:
	int preIndex;
	BasicMatrix* p_UpTriangleMatrix;
	BasicMatrix* p_SubMatrix;
	BasicMatrix* p_SubInvserMatrix;
	BasicMatrix* p_TempMatrix;

	MatrixMultiplier m_Multiplier;
	MatrixInverser m_Inverser;

	//根据指定的对角线索引，计算对应对角线元素(特征值)的特征向量
	void getEigenVectorWithoutPreCheck(int diagonalIndex, BasicVector* p_resultVector);
};

#endif /* NORMALEIGENVECTORCALCFORREAL_H_ */
