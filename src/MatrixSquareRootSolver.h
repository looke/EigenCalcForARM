/*
 * MatrixSquareRootSolver.h
 * Solve square roots based on Denman-Beavers Iteration.
 *
 *  Created on: 2017��2��25��
 *      Author: looke
 */

#ifndef MATRIXSQUAREROOTSOLVER_H_
#define MATRIXSQUAREROOTSOLVER_H_
#include "BasicMatrix.h"
#include "MatrixInverser.h"

class MatrixSquareRootSolver
{
public:
	MatrixSquareRootSolver(BasicMatrix* input_opYMatrix, BasicMatrix* input_opZMatrix, BasicMatrix* input_YMatrix_Temp, BasicMatrix* input_ZMatrix_Temp, BasicMatrix* input_tempMatrix);

	virtual void init(BasicMatrix* input_opYMatrix, BasicMatrix* input_opZMatrix, BasicMatrix* input_YMatrix_Temp, BasicMatrix* input_ZMatrix_Temp, BasicMatrix* input_tempMatrix);
	virtual void reload(BasicMatrix* input_opYMatrix, BasicMatrix* input_opZMatrix, BasicMatrix* input_YMatrix_Temp, BasicMatrix* input_ZMatrix_Temp, BasicMatrix* input_tempMatrix);

	bool generateSquareRootMatrix();

	void DenmanBeaversIteration();

	//BasicMatrix* getSquareRootMatrix();


	//void printYMatrix();
	//void printZMatrix();


	int getIterationTime();


	virtual ~MatrixSquareRootSolver() {};

	//��������
	MatrixInverser m_inverser;


protected:
	int iterationTime;
	int maxIterationTime; //������������������������Ȼû�л��Ŀ�꾫�ȵ����ݣ�����Ϊ���ʧ��

	bool isSquareMatrix;
	int rowNumber;
	int columnNumber;
	double maxDiff;
	double threshold;

	BasicMatrix* p_YMatrix;
	BasicMatrix* p_ZMatrix;

	//�����м���̼������ľ���
	BasicMatrix* p_tempYMatrix;
	BasicMatrix* p_tempZMatrix;

	//�����м���̼������ľ���
	BasicMatrix* p_tempMatrix;

};

#endif /* MATRIXINVERSER_H_ */
