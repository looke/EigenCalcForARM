/*
 * QRDecomposition.h
 *
 *  Created on: 2017��4��25��
 *      Author: looke
 */

#ifndef TRANSFORMATION_BASIC_QRDECOMPOSITION_H_
#define TRANSFORMATION_BASIC_QRDECOMPOSITION_H_

#include "BasicMatrix.h"
#include "MatrixMultiplier.h"
#include "HouseholderTransformation.h"

class QRDecomposition
{

public:

	//QRDecomposition();
	QRDecomposition(BasicMatrix* p_input_OpMatrix, BasicMatrix* p_input_QMatrix, BasicMatrix* p_input_householderMatrix, BasicMatrix* p_input_TempMatrix);

	void init(BasicMatrix* p_input_OpMatrix, BasicMatrix* p_input_QMatrix, BasicMatrix* p_input_householderMatrix, BasicMatrix* p_input_TempMatrix);
	void reload(BasicMatrix* p_input_OpMatrix, BasicMatrix* p_input_QMatrix, BasicMatrix* p_input_householderMatrix, BasicMatrix* p_input_TempMatrix);

	void calcQRMatrix();

	//���²��������Q����
	void updateMatrix();

	//����HouseHolder�任����
	void generateHouseholderMatrix(int index);


	BasicMatrix* getQMatrix();
	BasicMatrix* getQTMatrix();
	BasicMatrix* getRMatrix();
	BasicMatrix* getHouseholderMatrix();


	//virtual ~QRDecomposition(){};

private:

protected:

	//�������󣬾��������󽫱���϶ԽǾ���
	BasicMatrix* p_OpMatrix;


	//Q����
	BasicMatrix* p_QMatrix;

	//Q�����ת��(�����)
	//BasicMatrix* p_QTMatrix;


	//householder����ÿ�ε���ʹ��
	BasicMatrix* p_householderMatrix;

	//�м���̾���
	BasicMatrix* p_TempMatrix;

	//�˷���
	MatrixMultiplier m_Multiplier;

	//householder�任
	HouseholderTransformation m_HouseholderTrans;

};


#endif /* TRANSFORMATION_BASIC_QRDECOMPOSITION_H_ */
