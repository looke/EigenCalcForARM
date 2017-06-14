/*
 * HessenbergTriangleFormular.h
 *
 *  Created on: 2017��4��28��
 *      Author: looke
 */

#ifndef TRANSFORMATION_BASIC_HESSENBERGTRIANGLEFORMULAR_H_
#define TRANSFORMATION_BASIC_HESSENBERGTRIANGLEFORMULAR_H_

#include "BasicMatrix.h"
#include "MatrixMultiplier.h"
#include "QRDecomposition.h"
#include "GivensTransformation.h"


class HessenbergTriangleFormular
{
public:
	//HessenbergTriangleFormular();
	HessenbergTriangleFormular(BasicMatrix* input_Matrix_A, BasicMatrix* input_Matrix_B,BasicMatrix* input_QMatrix_Total,BasicMatrix* input_ZMatrix_Total,BasicMatrix* input_QZMatrix_Step,BasicMatrix* input_TempMatrix_Trans,BasicMatrix* input_TempMatrix);

	void initABMatrix();
	void formularABMatrix();
	void formularColumnVector(int columnIndex);

	void updateOpMatrix_A_ByQ();
	void updateOpMatrix_B_ByQ();
	void updateOpMatrixByZ();

	void updateQMatrix_Total();
	void updateZMatrix_Total();

	BasicMatrix* getMatrixA();
	BasicMatrix* getMatrixB();
	BasicMatrix* getHessenbergMatrixA();
	BasicMatrix* getTriangleMatrixB();
	//BasicMatrix* getMatrixQ_Step();
	//BasicMatrix* getMatrixZ_Step();
	BasicMatrix* getMatrixQZ_Step();
	BasicMatrix* getMatrixQ_Total();
	BasicMatrix* getMatrixZ_Total();


	void init(BasicMatrix* input_Matrix_A, BasicMatrix* input_Matrix_B,BasicMatrix* input_QMatrix_Total,BasicMatrix* input_ZMatrix_Total,BasicMatrix* input_QZMatrix_Step,BasicMatrix* input_TempMatrix_Trans,BasicMatrix* input_TempMatrix);
	void reload(BasicMatrix* input_Matrix_A, BasicMatrix* input_Matrix_B,BasicMatrix* input_QMatrix_Total,BasicMatrix* input_ZMatrix_Total,BasicMatrix* input_QZMatrix_Step,BasicMatrix* input_TempMatrix_Trans,BasicMatrix* input_TempMatrix);

	//virtual ~HessenbergTriangleFormular(){};
protected:

	//ԭʼ��������A
	BasicMatrix* p_OpMatrix_A;
	//ԭʼ��������B
	BasicMatrix* p_OpMatrix_B;

	//Hessenberg����A
	//BasicMatrix* p_OpHessenbergMatrix_A;

	//�����Ƿ���B
	//BasicMatrix* p_OpTriangleMatrix_B;

	//����Q�任����,���ڻ�����˲�������
	BasicMatrix* p_QMatrix_Total;

	//����Z�任����,���ڻ����ҳ˲�������
	BasicMatrix* p_ZMatrix_Total;

	//����Q�任����,������˲�������
	//BasicMatrix* p_QMatrix_Step;

	//����Z�任����,�����ҳ˲�������
	//BasicMatrix* p_ZMatrix_Step;

	//Q* A Z--->Hessenberg����
	//Q* B Z--->�����Ǿ���

	//����QZ�任����,������˲�������
	BasicMatrix* p_QZMatrix_Step;

	//�м�ת������
	BasicMatrix* p_TempMatrix_Trans;

	//�м���̾���
	BasicMatrix* p_TempMatrix;

	//QR�ֽ�
	QRDecomposition m_QRDecomp;

	//Givens�任
	GivensTransformation m_GivensTrans;

	//�˷���
	MatrixMultiplier m_Multiplier;

	//ת����
	MatrixTransposer m_Transposer;

private:

};

#endif /* TRANSFORMATION_BASIC_HESSENBERGTRIANGLEFORMULAR_H_ */
