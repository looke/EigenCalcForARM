/*
 * HessenbergFormular.h
 *
 *  Created on: 2017��5��8��
 *      Author: looke
 */

#ifndef TRANSFORMATION_BASIC_HESSENBERGFORMULAR_H_
#define TRANSFORMATION_BASIC_HESSENBERGFORMULAR_H_

#include "BasicMatrix.h"
#include "BasicVector.h"
#include "HouseholderTransformation.h"
#include "MatrixMultiplier.h"
#include "MatrixTransposer.h"

class HessenbergFormular
{
public:
	//HessenbergFormular();
	HessenbergFormular(BasicMatrix* p_Input_OpMatrix, BasicMatrix* p_Input_preTransMatrix,BasicMatrix* p_Input_transMatrix, BasicMatrix* p_Input_tempMatrix);

	void init(BasicMatrix* p_Input_OpMatrix, BasicMatrix* p_Input_preTransMatrix,BasicMatrix* p_Input_transMatrix, BasicMatrix* p_Input_tempMatrix);
	void reload(BasicMatrix* p_Input_OpMatrix, BasicMatrix* p_Input_preTransMatrix,BasicMatrix* p_Input_transMatrix, BasicMatrix* p_Input_tempMatrix);
	//void resizeSubMatrix(int rowAndColumnNumber);

	//���ݵ�ǰ����  ��ʼ���ӱ任����
	//void initSubHouseholderTrans(int iterateNum);

	//Hessenberg��ʽ��
	void formularUpperHessnbergMatrix();

	//���ݵ�ǰ���� ���������ӱ任����
	void generateSubHouseholderTrans(int iterateNum);

	//���ӱ任��������Ϊȫά�ȱ任��
	void upgradeSubHouseholderTrans(int iterateNum);

	void updatePreTransMatrix();
	//void updateAfterTransMatrix();

	void updateOpMatrix();

	BasicMatrix* getOpMatrix();
	BasicMatrix* getPreTransMatrix();
	//BasicMatrix* getAfterTransMatrix();
	BasicMatrix* getTransMatrix();
	//BasicMatrix* getSubTransMatrix();

	//virtual ~HessenbergFormular(){};
protected:
	BasicMatrix* p_OpMatrix;

	//����任��--���
	BasicMatrix* p_preTransMatrix;

	//�ҳ˱任�� Ϊ��˱任���ת��
	//BasicMatrix* p_afterTransMatrix;

	//Householder�任��
	BasicMatrix* p_transMatrix;

	//Householder�ӱ任��
	//BasicMatrix* p_subTransMatrix;

	//�м���̾���
	BasicMatrix* p_tempMatrix;

	//Householder�任
	HouseholderTransformation m_HouseholderTrans;

	//�˷���
	MatrixMultiplier m_Multiplier;

	//ת����
	MatrixTransposer m_Transposer;
};



#endif /* TRANSFORMATION_BASIC_HESSENBERGFORMULAR_H_ */
