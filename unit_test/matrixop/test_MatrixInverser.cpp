/*
 * test_MatrixInverser.cpp
 *
 *  Created on: 2017年6月19日
 *      Author: looke
 */
#include "..\gtest_src\gtest\gtest.h"
#include "StaticMatrix.h"
#include "MatrixMultiplier.h"
#include "MatrixInverser.h"
#include "math.h"

/*
 * 测试MatrixInverser创建
 */
TEST(MatrixInverserInversTest, postive)
{
	StaticMatrix testMatrix33 = StaticMatrix(3,3);
	testMatrix33.setMatrixElement(0,0,1);
	testMatrix33.setMatrixElement(0,1,1);
	testMatrix33.setMatrixElement(0,2,1);

	testMatrix33.setMatrixElement(1,0,0);
	testMatrix33.setMatrixElement(1,1,1);
	testMatrix33.setMatrixElement(1,2,1);

	testMatrix33.setMatrixElement(2,0,0);
	testMatrix33.setMatrixElement(2,1,1);
	testMatrix33.setMatrixElement(2,2,0);

	StaticMatrix testMatrix33_Original = StaticMatrix(3,3);
	testMatrix33_Original.copyMatrixElementNoCheck(&testMatrix33);

	StaticMatrix testMatrix33_inv = StaticMatrix(3,3);
	testMatrix33_inv.resetMatrixToI();
	StaticMatrix testMatrix33_multi = StaticMatrix(3,3);

	MatrixInverser inverser = MatrixInverser(&testMatrix33,&testMatrix33_inv);
	inverser.generateInverseMatrix();

	MatrixMultiplier multiplier = MatrixMultiplier(&testMatrix33_Original,&testMatrix33_inv,&testMatrix33_multi);
	multiplier.multiplyCalc();

	//testMatrix33_multi.printMatrix();

	EXPECT_EQ(1, testMatrix33_multi.getMatrixElement(0,0));
	EXPECT_EQ(1, testMatrix33_multi.getMatrixElement(1,1));
	EXPECT_EQ(1, testMatrix33_multi.getMatrixElement(2,2));
}



