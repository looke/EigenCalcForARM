/*
 * test_MatrixSquareRootSolver.cpp
 *
 *  Created on: 2017年6月27日
 *      Author: looke
 */
#include "..\gtest_src\gtest\gtest.h"
#include "StaticMatrix.h"
#include "MatrixMultiplier.h"
#include "MatrixSquareRootSolver.h"
#include "math.h"

/*
 * 测试MatrixSquareRootSolver 3x3单位方阵 开方运算
 */

TEST(MatrixSquareRootSolverTest_Square_I, postive)
{
	StaticMatrix testMatrix33 = StaticMatrix(3,3);
	testMatrix33.setMatrixElement(0,0,1);
	testMatrix33.setMatrixElement(0,1,0);
	testMatrix33.setMatrixElement(0,2,0);

	testMatrix33.setMatrixElement(1,0,0);
	testMatrix33.setMatrixElement(1,1,1);
	testMatrix33.setMatrixElement(1,2,0);

	testMatrix33.setMatrixElement(2,0,0);
	testMatrix33.setMatrixElement(2,1,0);
	testMatrix33.setMatrixElement(2,2,1);

	StaticMatrix testMatrix33_Original = StaticMatrix(3,3);
	testMatrix33_Original.copyMatrixElementNoCheck(&testMatrix33);

	StaticMatrix testMatrix33_Z = StaticMatrix(3,3);
	StaticMatrix testMatrix33_Y_Temp = StaticMatrix(3,3);
	StaticMatrix testMatrix33_Z_Temp = StaticMatrix(3,3);
	StaticMatrix testMatrix33_Temp = StaticMatrix(3,3);

	StaticMatrix testMatrix33_Verify = StaticMatrix(3,3);

	MatrixSquareRootSolver squareCalc = MatrixSquareRootSolver(&testMatrix33,&testMatrix33_Z,&testMatrix33_Y_Temp,&testMatrix33_Z_Temp,&testMatrix33_Temp);
	bool result = squareCalc.generateSquareRootMatrix();

	EXPECT_EQ(true, result);

	MatrixMultiplier m_MultiCalc = MatrixMultiplier(&testMatrix33,&testMatrix33, &testMatrix33_Verify);
	m_MultiCalc.multiplyCalc();

	double maxDiff = testMatrix33_Original.calcMaxDifferentialNoCheck(&testMatrix33_Verify);
	double lowEdge = testMatrix33_Original.getLowEdge();

	EXPECT_GT(lowEdge, maxDiff);
}

/*
 * 测试MatrixSquareRootSolver 3x3普通方阵 开方运算
 */

TEST(MatrixSquareRootSolverTest_Square_Norm, postive)
{
	StaticMatrix testMatrix33 = StaticMatrix(3,3);
	testMatrix33.setMatrixElement(0,0,1);
	testMatrix33.setMatrixElement(0,1,0);
	testMatrix33.setMatrixElement(0,2,2);

	testMatrix33.setMatrixElement(1,0,0);
	testMatrix33.setMatrixElement(1,1,1);
	testMatrix33.setMatrixElement(1,2,0);

	testMatrix33.setMatrixElement(2,0,0);
	testMatrix33.setMatrixElement(2,1,0);
	testMatrix33.setMatrixElement(2,2,1);

	StaticMatrix testMatrix33_Original = StaticMatrix(3,3);
	testMatrix33_Original.copyMatrixElementNoCheck(&testMatrix33);

	StaticMatrix testMatrix33_Z = StaticMatrix(3,3);
	StaticMatrix testMatrix33_Y_Temp = StaticMatrix(3,3);
	StaticMatrix testMatrix33_Z_Temp = StaticMatrix(3,3);
	StaticMatrix testMatrix33_Temp = StaticMatrix(3,3);

	StaticMatrix testMatrix33_Verify = StaticMatrix(3,3);

	MatrixSquareRootSolver squareCalc = MatrixSquareRootSolver(&testMatrix33,&testMatrix33_Z,&testMatrix33_Y_Temp,&testMatrix33_Z_Temp, &testMatrix33_Temp);
	bool result = squareCalc.generateSquareRootMatrix();

	EXPECT_EQ(true, result);

	MatrixMultiplier m_MultiCalc = MatrixMultiplier(&testMatrix33,&testMatrix33, &testMatrix33_Verify);
	m_MultiCalc.multiplyCalc();

	double maxDiff = testMatrix33_Original.calcMaxDifferentialNoCheck(&testMatrix33_Verify);
	double lowEdge = testMatrix33_Original.getLowEdge();

	EXPECT_GT(lowEdge, maxDiff);
}
