/*
 * test_SingleShiftQRIteration.cpp
 *
 *  Created on: 2017年6月13日
 *      Author: looke
 */
#include "..\gtest_src\gtest\gtest.h"
#include "StaticVector.h"
#include "StaticMatrix.h"
#include "MatrixMultiplier.h"
#include "MatrixTransposer.h"
#include "QRDecomposition.h"
#include "SingleShiftQRIteration.h"
#include "math.h"
#include <iostream>
using namespace std;

/*
 * 测试SingleShiftQRIteration 隐式迭代
 */
TEST(SingleShiftQRIterationIMITTest_Normal4x4, postive)
{
	MatrixTransposer m_Transposer = MatrixTransposer();
	double lowEdge, maxDiff;
	StaticMatrix test33 = StaticMatrix(3,3);

	//特征值 9 -1 0
	test33.setMatrixElement(0,0,1);
	test33.setMatrixElement(0,1,2);
	test33.setMatrixElement(0,2,3);

	test33.setMatrixElement(1,0,2);
	test33.setMatrixElement(1,1,1);
	test33.setMatrixElement(1,2,3);

	test33.setMatrixElement(2,0,3);
	test33.setMatrixElement(2,1,3);
	test33.setMatrixElement(2,2,6);
	StaticMatrix test33_Original = StaticMatrix(3,3);
	test33_Original.copyMatrixElementNoCheck(&test33);

	StaticMatrix test33_QTMatrix = StaticMatrix(3,3);
	StaticMatrix test33_QMatrix = StaticMatrix(3,3);
	StaticMatrix test33_QQTMatrix_Step = StaticMatrix(3,3);
	StaticMatrix test33_TempMatrix = StaticMatrix(3,3);

	SingleShiftQRIteration singleQRIt = SingleShiftQRIteration(&test33,&test33_QTMatrix,&test33_QQTMatrix_Step,&test33_TempMatrix);
	//singleQRIt.rayleigh_Quotient_IM_QRIteration(10);
	cout << "After 10 Iteration: OpMatrix" << endl;
	test33.printMatrix();
	lowEdge = test33.getLowEdge();
	EXPECT_GT(lowEdge, fabs(test33.getMatrixElement(0,0)-9));
	EXPECT_GT(lowEdge, fabs(test33.getMatrixElement(1,1)+1));
	EXPECT_GT(lowEdge, fabs(test33.getMatrixElement(2,2)-0));

	MatrixMultiplier m_Multi = MatrixMultiplier(&test33_QTMatrix, &test33_Original, &test33_TempMatrix);
	m_Multi.multiplyCalc();
	test33_Original.copyMatrixElementNoCheck(&test33_TempMatrix);

	test33_QMatrix.copyMatrixElementNoCheck(&test33_QTMatrix);
	m_Transposer.transposeSquareMatrix(&test33_QMatrix);
	m_Multi.reload(&test33_Original,&test33_QMatrix,&test33_TempMatrix);
	m_Multi.multiplyCalc();
	test33_Original.copyMatrixElementNoCheck(&test33_TempMatrix);

	cout << "QT * OP *Q" << endl;
	test33_Original.printMatrix();
	maxDiff = test33_Original.calcMaxDifferentialNoCheck(&test33);
	EXPECT_GE(lowEdge, maxDiff);

	StaticMatrix test44 = StaticMatrix(4,4);

	test44.setMatrixElement(0,0,54);
	test44.setMatrixElement(0,1,40);
	test44.setMatrixElement(0,2,10);
	test44.setMatrixElement(0,3,76);

	test44.setMatrixElement(1,0,47);
	test44.setMatrixElement(1,1,20);
	test44.setMatrixElement(1,2,94);
	test44.setMatrixElement(1,3,49);

	test44.setMatrixElement(2,0,26);
	test44.setMatrixElement(2,1,80);
	test44.setMatrixElement(2,2,94);
	test44.setMatrixElement(2,3,70);

	test44.setMatrixElement(3,0,3);
	test44.setMatrixElement(3,1,92);
	test44.setMatrixElement(3,2,83);
	test44.setMatrixElement(3,3,45);
	StaticMatrix test44_Original = StaticMatrix(4,4);
	test44_Original.copyMatrixElementNoCheck(&test44);
	StaticMatrix test44_QTMatrix = StaticMatrix(4,4);
	StaticMatrix test44_QMatrix = StaticMatrix(4,4);
	StaticMatrix test44_QQTMatrix_Step = StaticMatrix(4,4);
	StaticMatrix test44_TempMatrix = StaticMatrix(4,4);

	singleQRIt.reload(&test44, &test44_QTMatrix, &test44_QQTMatrix_Step, &test44_TempMatrix);
	singleQRIt.rayleigh_Quotient_IM_QRIteration(80);

	cout << "After 10 Iteration: OpMatrix" << endl;
	test44.printMatrix();

	m_Multi.reload(&test44_QTMatrix, &test44_Original, &test44_TempMatrix);
	m_Multi.multiplyCalc();
	test44_Original.copyMatrixElementNoCheck(&test44_TempMatrix);

	test44_QMatrix.copyMatrixElementNoCheck(&test44_QTMatrix);
	m_Transposer.transposeSquareMatrix(&test44_QMatrix);

	m_Multi.reload(&test44_Original,&test44_QMatrix,&test44_TempMatrix);
	m_Multi.multiplyCalc();
	test44_Original.copyMatrixElementNoCheck(&test44_TempMatrix);
	lowEdge = test44_Original.getLowEdge();
	maxDiff = test44_Original.calcMaxDifferentialNoCheck(&test44);
	EXPECT_GE(lowEdge, maxDiff);

}
