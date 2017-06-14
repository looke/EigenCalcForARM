/*
 * QRDecomposition.h
 *
 *  Created on: 2017年4月25日
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

	//更新操作矩阵和Q矩阵
	void updateMatrix();

	//生成HouseHolder变换矩阵
	void generateHouseholderMatrix(int index);


	BasicMatrix* getQMatrix();
	BasicMatrix* getQTMatrix();
	BasicMatrix* getRMatrix();
	BasicMatrix* getHouseholderMatrix();


	//virtual ~QRDecomposition(){};

private:

protected:

	//操作矩阵，经过迭代后将变成上对角矩阵
	BasicMatrix* p_OpMatrix;


	//Q矩阵
	BasicMatrix* p_QMatrix;

	//Q矩阵的转置(逆矩阵)
	//BasicMatrix* p_QTMatrix;


	//householder矩阵，每次迭代使用
	BasicMatrix* p_householderMatrix;

	//中间过程矩阵
	BasicMatrix* p_TempMatrix;

	//乘法器
	MatrixMultiplier m_Multiplier;

	//householder变换
	HouseholderTransformation m_HouseholderTrans;

};


#endif /* TRANSFORMATION_BASIC_QRDECOMPOSITION_H_ */
