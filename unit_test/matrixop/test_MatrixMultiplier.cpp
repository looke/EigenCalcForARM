/*
 * test_StaticMatrixMultiplier.cpp
 *
 *  Created on: 2017年6月8日
 *      Author: looke
 */

#include "..\gtest_src\gtest\gtest.h"
#include "StaticMatrix.h"
#include "MatrixMultiplier.h"
#include "math.h"

/*
 * 测试StaticMatrixMultiplier创建异常处理
 */
/*
TEST(StaticMatrixMultiplierCreateExceptionTest, postive)
{
	StaticMatrix leftMatrix32 = StaticMatrix(3,2);
	StaticMatrix rightMatrix33 = StaticMatrix(3,3);
	StaticMatrix resultMatrix33 = StaticMatrix(3,3);
	EXPECT_EQ(0, MatrixMultiplier(&leftMatrix32, &rightMatrix33, &resultMatrix33));

	StaticMatrix leftMatrix42 = StaticMatrix(4,2);
	StaticMatrix rightMatrix22 = StaticMatrix(2,2);
	StaticMatrix resultMatrix44 = StaticMatrix(4,4);
	EXPECT_EQ(0,MatrixMultiplier(&leftMatrix42, &rightMatrix22, &resultMatrix44));

	StaticMatrix leftMatrix53 = StaticMatrix(5,3);
	StaticMatrix rightMatrix34 = StaticMatrix(3,4);
	StaticMatrix resultMatrix24 = StaticMatrix(2,4);
	EXPECT_EQ(0,MatrixMultiplier(&leftMatrix53, &rightMatrix34, &resultMatrix24));

	StaticMatrix leftMatrix65 = StaticMatrix(6,5);
	StaticMatrix rightMatrix45 = StaticMatrix(4,5);
	StaticMatrix resultMatrix65 = StaticMatrix(6,5);
	EXPECT_EQ(0,MatrixMultiplier(&leftMatrix65, &rightMatrix45, &resultMatrix65));
}
*/
/*
 * 测试MatrixMultiplier创建
 */
TEST(MatrixMultiplierCreateTest, postive)
{
	StaticMatrix leftMatrix33 = StaticMatrix(3,3);
	StaticMatrix rightMatrix33 = StaticMatrix(3,3);
	StaticMatrix resultMatrix33 = StaticMatrix(3,3);
	EXPECT_NO_THROW(MatrixMultiplier(&leftMatrix33, &rightMatrix33, &resultMatrix33));

	StaticMatrix leftMatrix42 = StaticMatrix(4,2);
	StaticMatrix rightMatrix22 = StaticMatrix(2,2);
	StaticMatrix resultMatrix42 = StaticMatrix(4,2);
	EXPECT_NO_THROW(MatrixMultiplier(&leftMatrix42, &rightMatrix22, &resultMatrix42));

	StaticMatrix leftMatrix53 = StaticMatrix(5,3);
	StaticMatrix rightMatrix34 = StaticMatrix(3,4);
	StaticMatrix resultMatrix54 = StaticMatrix(5,4);
	EXPECT_NO_THROW(MatrixMultiplier(&leftMatrix53, &rightMatrix34, &resultMatrix54));

	StaticMatrix leftMatrix65 = StaticMatrix(6,5);
	StaticMatrix rightMatrix55 = StaticMatrix(5,5);
	StaticMatrix resultMatrix65 = StaticMatrix(6,5);
	EXPECT_NO_THROW(MatrixMultiplier(&leftMatrix65, &rightMatrix55, &resultMatrix65));
}


/*
 * 测试MatrixMultiplier 矩阵乘法
 */
TEST(MatrixMultiplierMultiCalcTest_IxI, postive)
{
	StaticMatrix leftMatrix33 = StaticMatrix(3,3);
	leftMatrix33.resetMatrixToI();
	StaticMatrix rightMatrix33 = StaticMatrix(3,3);
	rightMatrix33.resetMatrixToI();

	StaticMatrix resultMatrix33 = StaticMatrix(3,3);

	MatrixMultiplier(&leftMatrix33, &rightMatrix33, &resultMatrix33).multiplyCalc();

	EXPECT_EQ(1, rightMatrix33.getMatrixElement(0,0));
	EXPECT_EQ(0, rightMatrix33.getMatrixElement(0,1));
	EXPECT_EQ(0, rightMatrix33.getMatrixElement(0,2));

	EXPECT_EQ(0, rightMatrix33.getMatrixElement(1,0));
	EXPECT_EQ(1, rightMatrix33.getMatrixElement(1,1));
	EXPECT_EQ(0, rightMatrix33.getMatrixElement(1,2));

	EXPECT_EQ(0, rightMatrix33.getMatrixElement(2,0));
	EXPECT_EQ(0, rightMatrix33.getMatrixElement(2,1));
	EXPECT_EQ(1, rightMatrix33.getMatrixElement(2,2));
}

/*
 * 测试MatrixMultiplier 矩阵乘法
 */
TEST(MatrixMultiplierMultiCalcTest_NxN, postive)
{
	StaticMatrix leftMatrix44 = StaticMatrix(4,4);
	leftMatrix44.setMatrixElement(0,0,1);
	leftMatrix44.setMatrixElement(0,1,2);
	leftMatrix44.setMatrixElement(0,2,3);
	leftMatrix44.setMatrixElement(0,3,4);

	leftMatrix44.setMatrixElement(1,0,5);
	leftMatrix44.setMatrixElement(1,1,6);
	leftMatrix44.setMatrixElement(1,2,7);
	leftMatrix44.setMatrixElement(1,3,8);

	leftMatrix44.setMatrixElement(2,0,9);
	leftMatrix44.setMatrixElement(2,1,10);
	leftMatrix44.setMatrixElement(2,2,11);
	leftMatrix44.setMatrixElement(2,3,12);

	leftMatrix44.setMatrixElement(3,0,13);
	leftMatrix44.setMatrixElement(3,1,14);
	leftMatrix44.setMatrixElement(3,2,15);
	leftMatrix44.setMatrixElement(3,3,16);

	StaticMatrix rightMatrix43 = StaticMatrix(4,3);
	rightMatrix43.setMatrixElement(0,0,1);
	rightMatrix43.setMatrixElement(0,1,2);
	rightMatrix43.setMatrixElement(0,2,3);

	rightMatrix43.setMatrixElement(1,0,1);
	rightMatrix43.setMatrixElement(1,1,2);
	rightMatrix43.setMatrixElement(1,2,3);

	rightMatrix43.setMatrixElement(2,0,1);
	rightMatrix43.setMatrixElement(2,1,2);
	rightMatrix43.setMatrixElement(2,2,3);

	rightMatrix43.setMatrixElement(3,0,1);
	rightMatrix43.setMatrixElement(3,1,2);
	rightMatrix43.setMatrixElement(3,2,3);

	StaticMatrix resultMatrix43 = StaticMatrix(4,3);

	EXPECT_EQ(true,MatrixMultiplier(&leftMatrix44, &rightMatrix43, &resultMatrix43).multiplyCalc());

	EXPECT_EQ(10, resultMatrix43.getMatrixElement(0,0));
	EXPECT_EQ(20, resultMatrix43.getMatrixElement(0,1));
	EXPECT_EQ(30, resultMatrix43.getMatrixElement(0,2));

	EXPECT_EQ(26, resultMatrix43.getMatrixElement(1,0));
	EXPECT_EQ(52, resultMatrix43.getMatrixElement(1,1));
	EXPECT_EQ(78, resultMatrix43.getMatrixElement(1,2));

	EXPECT_EQ(42, resultMatrix43.getMatrixElement(2,0));
	EXPECT_EQ(84, resultMatrix43.getMatrixElement(2,1));
	EXPECT_EQ(126, resultMatrix43.getMatrixElement(2,2));

	EXPECT_EQ(58, resultMatrix43.getMatrixElement(3,0));
	EXPECT_EQ(116, resultMatrix43.getMatrixElement(3,1));
	EXPECT_EQ(174, resultMatrix43.getMatrixElement(3,2));
}

/*
 * 测试MatrixMultiplier 矩阵乘法 维度不匹配
 */
TEST(MatrixMultiplierMultiCalcExceptionTest, postive)
{
	StaticMatrix leftMatrix32 = StaticMatrix(3,2);
	StaticMatrix rightMatrix33 = StaticMatrix(3,3);
	StaticMatrix resultMatrix33 = StaticMatrix(3,3);
	bool result = false;
	EXPECT_EQ(result, MatrixMultiplier(&leftMatrix32, &rightMatrix33, &resultMatrix33).multiplyCalc());
}
