/*
 * test_DoubleShiftQRIteration.cpp
 *
 *  Created on: 2017年6月14日
 *      Author: looke
 */

#include "..\gtest_src\gtest\gtest.h"
#include "StaticVector.h"
#include "StaticMatrix.h"
#include "MatrixMultiplier.h"
#include "MatrixTransposer.h"
#include "QRDecomposition.h"
#include "DoubleShiftQRIteration.h"
#include "math.h"
#include <iostream>
using namespace std;

/*
 * 测试DoubleShiftQRIteration 隐式迭代
 */
TEST(DoubleShiftQRIterationIMITTest_Normal4x4, postive)
{
	double lowEdge, maxDiff;

	StaticMatrix test44 = StaticMatrix(4,4);

	test44.setMatrixElement(0,0,1);
	test44.setMatrixElement(0,1,2);
	test44.setMatrixElement(0,2,3);
	test44.setMatrixElement(0,3,4);

	test44.setMatrixElement(1,0,5);
	test44.setMatrixElement(1,1,6);
	test44.setMatrixElement(1,2,7);
	test44.setMatrixElement(1,3,8);

	test44.setMatrixElement(2,0,0);
	test44.setMatrixElement(2,1,9);
	test44.setMatrixElement(2,2,10);
	test44.setMatrixElement(2,3,11);

	test44.setMatrixElement(3,0,0);
	test44.setMatrixElement(3,1,0);
	test44.setMatrixElement(3,2,-12);
	test44.setMatrixElement(3,3,13);
	StaticMatrix test44_Original = StaticMatrix(4,4);
	test44_Original.copyMatrixElementNoCheck(&test44);

	StaticVector test44_vector  = StaticVector(4);
	StaticMatrix test44_QMatrix = StaticMatrix(4,4);
	StaticMatrix test44_QTMatrix = StaticMatrix(4,4);
	StaticMatrix test44_QQTMatrix_Step = StaticMatrix(4,4);
	StaticMatrix test44_TempMatrix = StaticMatrix(4,4);

	DoubleShiftQRIteration dsIteration = DoubleShiftQRIteration(&test44,&test44_vector,&test44_QMatrix,&test44_QTMatrix,&test44_QQTMatrix_Step,&test44_TempMatrix);
	//dsIteration.wilkinson_IM_QRIteration();


	MatrixMultiplier m_Multiplier = MatrixMultiplier(&test44_QTMatrix,&test44_Original, &test44_TempMatrix);
	m_Multiplier.multiplyCalc();
	test44_Original.copyMatrixElementNoCheck(&test44_TempMatrix);
	m_Multiplier.reload(&test44_Original, &test44_QMatrix, &test44_TempMatrix);
	m_Multiplier.multiplyCalc();
	test44_Original.copyMatrixElementNoCheck(&test44_TempMatrix);

	lowEdge = test44_Original.getLowEdge();
	maxDiff = test44_Original.calcMaxDifferentialNoCheck(&test44);

	EXPECT_GT(lowEdge,maxDiff);
}
