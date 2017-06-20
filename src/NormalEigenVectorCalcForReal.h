/*
 * NormalEigenVectorCalcForReal.h
 *
 *  Created on: 2017��6��19��
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

	//����ָ���ĶԽ��������������Ӧ�Խ���Ԫ��(����ֵ)��������������ǰ��Ԫ��Ҫ�����Ų�
	bool getEigenVector(int diagonalIndex, BasicVector* p_resultVector);

	//����ָ���Խ������������������Ǿ���Խ���ǰ��Ԫ�����Ƿ񻹴�����ͬȡֵ������ֵ��������ڣ��򷵻ض�ӦԪ�صĶԽ�������
	bool findSameEigenValuePreIndex(int diagonalIndex);
protected:
	int preIndex;
	BasicMatrix* p_UpTriangleMatrix;
	BasicMatrix* p_SubMatrix;
	BasicMatrix* p_SubInvserMatrix;
	BasicMatrix* p_TempMatrix;

	MatrixMultiplier m_Multiplier;
	MatrixInverser m_Inverser;

	//����ָ���ĶԽ��������������Ӧ�Խ���Ԫ��(����ֵ)����������
	void getEigenVectorWithoutPreCheck(int diagonalIndex, BasicVector* p_resultVector);
};

#endif /* NORMALEIGENVECTORCALCFORREAL_H_ */
