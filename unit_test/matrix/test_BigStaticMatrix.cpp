/*
 * test_BigStaticMatrix.cpp
 *
 *  Created on: 2017年6月24日
 *      Author: looke
 */

#include "..\gtest_src\gtest\gtest.h"
#include "BigStaticMatrix.h"
#include "math.h"

/*
 * 测试Matrix初始化 越界处理
 */
TEST(BigMatrixCreateExceptionTest, postive)
{
	BigStaticMatrix testMatrix_13 = BigStaticMatrix(-1,3);
	BigStaticMatrix testMatrix223 = BigStaticMatrix(2,23);
	BigStaticMatrix testMatrix2223 = BigStaticMatrix(22,23);
	BigStaticMatrix testMatrix213 = BigStaticMatrix(21,3);
	EXPECT_EQ(0,testMatrix_13.rowNum);
	EXPECT_EQ(0,testMatrix_13.columnNum);

	EXPECT_EQ(0,testMatrix223.rowNum);
	EXPECT_EQ(0,testMatrix223.columnNum);

	EXPECT_EQ(0,testMatrix2223.rowNum);
	EXPECT_EQ(0,testMatrix2223.columnNum);

	EXPECT_EQ(0,testMatrix213.rowNum);
	EXPECT_EQ(0,testMatrix213.columnNum);
}

/*
 * 测试Matrix初始化
 */
TEST(BigMatrixCreateTest, postive)
{
	BigStaticMatrix testMatrix33 = BigStaticMatrix(3,3);
	BigStaticMatrix testMatrix53 = BigStaticMatrix(5,3);
	BigStaticMatrix testMatrix00 = BigStaticMatrix(0,0);
	BigStaticMatrix testMatrix2020 = BigStaticMatrix(20,20);
	EXPECT_EQ(3,testMatrix33.rowNum);
	EXPECT_EQ(3,testMatrix33.columnNum);

	EXPECT_EQ(1,testMatrix33.getMatrixElement(0,0));
	EXPECT_EQ(2,testMatrix33.getMatrixElement(0,1));
	EXPECT_EQ(3,testMatrix33.getMatrixElement(0,2));

	EXPECT_EQ(4,testMatrix33.getMatrixElement(1,0));
	EXPECT_EQ(5,testMatrix33.getMatrixElement(1,1));
	EXPECT_EQ(6,testMatrix33.getMatrixElement(1,2));

	EXPECT_EQ(7,testMatrix33.getMatrixElement(2,0));
	EXPECT_EQ(8,testMatrix33.getMatrixElement(2,1));
	EXPECT_EQ(9,testMatrix33.getMatrixElement(2,2));

	EXPECT_EQ(5,testMatrix53.rowNum);
	EXPECT_EQ(3,testMatrix53.columnNum);

	EXPECT_EQ(1,testMatrix53.getMatrixElement(0,0));
	EXPECT_EQ(2,testMatrix53.getMatrixElement(0,1));
	EXPECT_EQ(3,testMatrix53.getMatrixElement(0,2));

	EXPECT_EQ(4,testMatrix53.getMatrixElement(1,0));
	EXPECT_EQ(5,testMatrix53.getMatrixElement(1,1));
	EXPECT_EQ(6,testMatrix53.getMatrixElement(1,2));

	EXPECT_EQ(7,testMatrix53.getMatrixElement(2,0));
	EXPECT_EQ(8,testMatrix53.getMatrixElement(2,1));
	EXPECT_EQ(9,testMatrix53.getMatrixElement(2,2));

	EXPECT_EQ(10,testMatrix53.getMatrixElement(3,0));
	EXPECT_EQ(11,testMatrix53.getMatrixElement(3,1));
	EXPECT_EQ(12,testMatrix53.getMatrixElement(3,2));

	EXPECT_EQ(13,testMatrix53.getMatrixElement(4,0));
	EXPECT_EQ(14,testMatrix53.getMatrixElement(4,1));
	EXPECT_EQ(15,testMatrix53.getMatrixElement(4,2));

	EXPECT_EQ(20,testMatrix2020.rowNum);
	EXPECT_EQ(20,testMatrix2020.columnNum);
}

/*
 * 测试Matrix 指定元素值设置
 */
TEST(BigMatrixSetElementTest, postive)
{
	BigStaticMatrix testMatrix33 = BigStaticMatrix(3,3);
	bool result = false;
	EXPECT_EQ(result,testMatrix33.setMatrixElement(-1,0,10));
	EXPECT_EQ(result,testMatrix33.setMatrixElement(0,-1,10));
	EXPECT_EQ(result,testMatrix33.setMatrixElement(3,2,10));
	EXPECT_EQ(result,testMatrix33.setMatrixElement(0,3,10));
	EXPECT_EQ(result,testMatrix33.setMatrixElement(5,3,10));
	EXPECT_EQ(result,testMatrix33.setMatrixElement(6,2,10));
	EXPECT_EQ(result,testMatrix33.setMatrixElement(1,20,10));
	EXPECT_EQ(result,testMatrix33.setMatrixElement(6,33,10));

	EXPECT_EQ(1,testMatrix33.getMatrixElement(0,0));
	EXPECT_EQ(2,testMatrix33.getMatrixElement(0,1));
	EXPECT_EQ(3,testMatrix33.getMatrixElement(0,2));

	EXPECT_EQ(4,testMatrix33.getMatrixElement(1,0));
	EXPECT_EQ(5,testMatrix33.getMatrixElement(1,1));
	EXPECT_EQ(6,testMatrix33.getMatrixElement(1,2));

	EXPECT_EQ(7,testMatrix33.getMatrixElement(2,0));
	EXPECT_EQ(8,testMatrix33.getMatrixElement(2,1));
	EXPECT_EQ(9,testMatrix33.getMatrixElement(2,2));

	result = true;
	EXPECT_EQ(result,testMatrix33.setMatrixElement(0,0,10));
	EXPECT_EQ(result,testMatrix33.setMatrixElement(0,1,11));
	EXPECT_EQ(result,testMatrix33.setMatrixElement(0,2,12));

	EXPECT_EQ(result,testMatrix33.setMatrixElement(1,0,13));
	EXPECT_EQ(result,testMatrix33.setMatrixElement(1,1,14));
	EXPECT_EQ(result,testMatrix33.setMatrixElement(1,2,15));

	EXPECT_EQ(result,testMatrix33.setMatrixElement(2,0,16));
	EXPECT_EQ(result,testMatrix33.setMatrixElement(2,1,17));
	EXPECT_EQ(result,testMatrix33.setMatrixElement(2,2,18));

	EXPECT_EQ(10,testMatrix33.getMatrixElement(0,0));
	EXPECT_EQ(11,testMatrix33.getMatrixElement(0,1));
	EXPECT_EQ(12,testMatrix33.getMatrixElement(0,2));

	EXPECT_EQ(13,testMatrix33.getMatrixElement(1,0));
	EXPECT_EQ(14,testMatrix33.getMatrixElement(1,1));
	EXPECT_EQ(15,testMatrix33.getMatrixElement(1,2));

	EXPECT_EQ(16,testMatrix33.getMatrixElement(2,0));
	EXPECT_EQ(17,testMatrix33.getMatrixElement(2,1));
	EXPECT_EQ(18,testMatrix33.getMatrixElement(2,2));
}

/*
 * 测试Matrix 交换行
 */
TEST(BigMatrixSwapRowTest, postive)
{
	BigStaticMatrix testMatrix_23 = BigStaticMatrix(2,3);
	bool result = false;
	EXPECT_EQ(result,testMatrix_23.swapRow(-1,0));
	EXPECT_EQ(result,testMatrix_23.swapRow(1,-2));
	EXPECT_EQ(result,testMatrix_23.swapRow(2,1));
	EXPECT_EQ(result,testMatrix_23.swapRow(0,3));

	EXPECT_EQ(1,testMatrix_23.getMatrixElement(0,0));
	EXPECT_EQ(2,testMatrix_23.getMatrixElement(0,1));
	EXPECT_EQ(3,testMatrix_23.getMatrixElement(0,2));

	EXPECT_EQ(4,testMatrix_23.getMatrixElement(1,0));
	EXPECT_EQ(5,testMatrix_23.getMatrixElement(1,1));
	EXPECT_EQ(6,testMatrix_23.getMatrixElement(1,2));

	result = true;
	EXPECT_EQ(result,testMatrix_23.swapRow(0,1));

	EXPECT_EQ(4,testMatrix_23.getMatrixElement(0,0));
	EXPECT_EQ(5,testMatrix_23.getMatrixElement(0,1));
	EXPECT_EQ(6,testMatrix_23.getMatrixElement(0,2));

	EXPECT_EQ(1,testMatrix_23.getMatrixElement(1,0));
	EXPECT_EQ(2,testMatrix_23.getMatrixElement(1,1));
	EXPECT_EQ(3,testMatrix_23.getMatrixElement(1,2));

	EXPECT_EQ(result,testMatrix_23.swapRow(1,0));

	EXPECT_EQ(1,testMatrix_23.getMatrixElement(0,0));
	EXPECT_EQ(2,testMatrix_23.getMatrixElement(0,1));
	EXPECT_EQ(3,testMatrix_23.getMatrixElement(0,2));

	EXPECT_EQ(4,testMatrix_23.getMatrixElement(1,0));
	EXPECT_EQ(5,testMatrix_23.getMatrixElement(1,1));
	EXPECT_EQ(6,testMatrix_23.getMatrixElement(1,2));
}

/*
 * 测试Matrix 交换列
 */
TEST(BigMatrixSwapColumnTest, postive)
{
	BigStaticMatrix testMatrix_23 = BigStaticMatrix(2,3);
	bool result = false;
	EXPECT_EQ(result,testMatrix_23.swapColumn(-1,0));
	EXPECT_EQ(result,testMatrix_23.swapColumn(1,-2));
	EXPECT_EQ(result,testMatrix_23.swapColumn(3,1));
	EXPECT_EQ(result,testMatrix_23.swapColumn(0,3));

	EXPECT_EQ(1,testMatrix_23.getMatrixElement(0,0));
	EXPECT_EQ(2,testMatrix_23.getMatrixElement(0,1));
	EXPECT_EQ(3,testMatrix_23.getMatrixElement(0,2));

	EXPECT_EQ(4,testMatrix_23.getMatrixElement(1,0));
	EXPECT_EQ(5,testMatrix_23.getMatrixElement(1,1));
	EXPECT_EQ(6,testMatrix_23.getMatrixElement(1,2));

	result = true;
	EXPECT_EQ(result,testMatrix_23.swapColumn(0,1));

	EXPECT_EQ(2,testMatrix_23.getMatrixElement(0,0));
	EXPECT_EQ(1,testMatrix_23.getMatrixElement(0,1));
	EXPECT_EQ(3,testMatrix_23.getMatrixElement(0,2));

	EXPECT_EQ(5,testMatrix_23.getMatrixElement(1,0));
	EXPECT_EQ(4,testMatrix_23.getMatrixElement(1,1));
	EXPECT_EQ(6,testMatrix_23.getMatrixElement(1,2));

	EXPECT_EQ(result,testMatrix_23.swapColumn(1,0));

	EXPECT_EQ(1,testMatrix_23.getMatrixElement(0,0));
	EXPECT_EQ(2,testMatrix_23.getMatrixElement(0,1));
	EXPECT_EQ(3,testMatrix_23.getMatrixElement(0,2));

	EXPECT_EQ(4,testMatrix_23.getMatrixElement(1,0));
	EXPECT_EQ(5,testMatrix_23.getMatrixElement(1,1));
	EXPECT_EQ(6,testMatrix_23.getMatrixElement(1,2));

	EXPECT_EQ(result,testMatrix_23.swapColumn(2,0));

	EXPECT_EQ(3,testMatrix_23.getMatrixElement(0,0));
	EXPECT_EQ(2,testMatrix_23.getMatrixElement(0,1));
	EXPECT_EQ(1,testMatrix_23.getMatrixElement(0,2));

	EXPECT_EQ(6,testMatrix_23.getMatrixElement(1,0));
	EXPECT_EQ(5,testMatrix_23.getMatrixElement(1,1));
	EXPECT_EQ(4,testMatrix_23.getMatrixElement(1,2));
}

/*
 * 测试Matrix 交换对角元
 */
TEST(BigMatrixSwapDiagonalTest, postive)
{
	BigStaticMatrix testMatrix_23 = BigStaticMatrix(2,3);
	bool result = false;
	EXPECT_EQ(result,testMatrix_23.swapDiagElement(0,2));
	EXPECT_EQ(result,testMatrix_23.swapDiagElement(1,2));

	EXPECT_EQ(1,testMatrix_23.getMatrixElement(0,0));
	EXPECT_EQ(2,testMatrix_23.getMatrixElement(0,1));
	EXPECT_EQ(3,testMatrix_23.getMatrixElement(0,2));

	EXPECT_EQ(4,testMatrix_23.getMatrixElement(1,0));
	EXPECT_EQ(5,testMatrix_23.getMatrixElement(1,1));
	EXPECT_EQ(6,testMatrix_23.getMatrixElement(1,2));

	BigStaticMatrix testMatrix44 = BigStaticMatrix(4,4);
	EXPECT_EQ(result,testMatrix44.swapDiagElement(-1,2));
	EXPECT_EQ(result,testMatrix44.swapDiagElement(4,0));
	EXPECT_EQ(result,testMatrix44.swapDiagElement(0,4));
	EXPECT_EQ(result,testMatrix44.swapDiagElement(6,3));

	EXPECT_EQ(1,testMatrix44.getMatrixElement(0,0));
	EXPECT_EQ(2,testMatrix44.getMatrixElement(0,1));
	EXPECT_EQ(3,testMatrix44.getMatrixElement(0,2));
	EXPECT_EQ(4,testMatrix44.getMatrixElement(0,3));

	EXPECT_EQ(5,testMatrix44.getMatrixElement(1,0));
	EXPECT_EQ(6,testMatrix44.getMatrixElement(1,1));
	EXPECT_EQ(7,testMatrix44.getMatrixElement(1,2));
	EXPECT_EQ(8,testMatrix44.getMatrixElement(1,3));

	EXPECT_EQ(9,testMatrix44.getMatrixElement(2,0));
	EXPECT_EQ(10,testMatrix44.getMatrixElement(2,1));
	EXPECT_EQ(11,testMatrix44.getMatrixElement(2,2));
	EXPECT_EQ(12,testMatrix44.getMatrixElement(2,3));

	EXPECT_EQ(13,testMatrix44.getMatrixElement(3,0));
	EXPECT_EQ(14,testMatrix44.getMatrixElement(3,1));
	EXPECT_EQ(15,testMatrix44.getMatrixElement(3,2));
	EXPECT_EQ(16,testMatrix44.getMatrixElement(3,3));

	result = true;
	EXPECT_EQ(result,testMatrix44.swapDiagElement(0,2));
	EXPECT_EQ(11,testMatrix44.getMatrixElement(0,0));
	EXPECT_EQ(2,testMatrix44.getMatrixElement(0,1));
	EXPECT_EQ(3,testMatrix44.getMatrixElement(0,2));
	EXPECT_EQ(4,testMatrix44.getMatrixElement(0,3));

	EXPECT_EQ(5,testMatrix44.getMatrixElement(1,0));
	EXPECT_EQ(6,testMatrix44.getMatrixElement(1,1));
	EXPECT_EQ(7,testMatrix44.getMatrixElement(1,2));
	EXPECT_EQ(8,testMatrix44.getMatrixElement(1,3));

	EXPECT_EQ(9,testMatrix44.getMatrixElement(2,0));
	EXPECT_EQ(10,testMatrix44.getMatrixElement(2,1));
	EXPECT_EQ(1,testMatrix44.getMatrixElement(2,2));
	EXPECT_EQ(12,testMatrix44.getMatrixElement(2,3));

	EXPECT_EQ(13,testMatrix44.getMatrixElement(3,0));
	EXPECT_EQ(14,testMatrix44.getMatrixElement(3,1));
	EXPECT_EQ(15,testMatrix44.getMatrixElement(3,2));
	EXPECT_EQ(16,testMatrix44.getMatrixElement(3,3));

	EXPECT_EQ(result,testMatrix44.swapDiagElement(0,2));
	EXPECT_EQ(1,testMatrix44.getMatrixElement(0,0));
	EXPECT_EQ(2,testMatrix44.getMatrixElement(0,1));
	EXPECT_EQ(3,testMatrix44.getMatrixElement(0,2));
	EXPECT_EQ(4,testMatrix44.getMatrixElement(0,3));

	EXPECT_EQ(5,testMatrix44.getMatrixElement(1,0));
	EXPECT_EQ(6,testMatrix44.getMatrixElement(1,1));
	EXPECT_EQ(7,testMatrix44.getMatrixElement(1,2));
	EXPECT_EQ(8,testMatrix44.getMatrixElement(1,3));

	EXPECT_EQ(9,testMatrix44.getMatrixElement(2,0));
	EXPECT_EQ(10,testMatrix44.getMatrixElement(2,1));
	EXPECT_EQ(11,testMatrix44.getMatrixElement(2,2));
	EXPECT_EQ(12,testMatrix44.getMatrixElement(2,3));

	EXPECT_EQ(13,testMatrix44.getMatrixElement(3,0));
	EXPECT_EQ(14,testMatrix44.getMatrixElement(3,1));
	EXPECT_EQ(15,testMatrix44.getMatrixElement(3,2));
	EXPECT_EQ(16,testMatrix44.getMatrixElement(3,3));

	EXPECT_EQ(result,testMatrix44.swapDiagElement(3,2));
	EXPECT_EQ(1,testMatrix44.getMatrixElement(0,0));
	EXPECT_EQ(2,testMatrix44.getMatrixElement(0,1));
	EXPECT_EQ(3,testMatrix44.getMatrixElement(0,2));
	EXPECT_EQ(4,testMatrix44.getMatrixElement(0,3));

	EXPECT_EQ(5,testMatrix44.getMatrixElement(1,0));
	EXPECT_EQ(6,testMatrix44.getMatrixElement(1,1));
	EXPECT_EQ(7,testMatrix44.getMatrixElement(1,2));
	EXPECT_EQ(8,testMatrix44.getMatrixElement(1,3));

	EXPECT_EQ(9,testMatrix44.getMatrixElement(2,0));
	EXPECT_EQ(10,testMatrix44.getMatrixElement(2,1));
	EXPECT_EQ(16,testMatrix44.getMatrixElement(2,2));
	EXPECT_EQ(12,testMatrix44.getMatrixElement(2,3));

	EXPECT_EQ(13,testMatrix44.getMatrixElement(3,0));
	EXPECT_EQ(14,testMatrix44.getMatrixElement(3,1));
	EXPECT_EQ(15,testMatrix44.getMatrixElement(3,2));
	EXPECT_EQ(11,testMatrix44.getMatrixElement(3,3));

	EXPECT_EQ(result,testMatrix44.swapDiagElement(2,3));
	EXPECT_EQ(1,testMatrix44.getMatrixElement(0,0));
	EXPECT_EQ(2,testMatrix44.getMatrixElement(0,1));
	EXPECT_EQ(3,testMatrix44.getMatrixElement(0,2));
	EXPECT_EQ(4,testMatrix44.getMatrixElement(0,3));

	EXPECT_EQ(5,testMatrix44.getMatrixElement(1,0));
	EXPECT_EQ(6,testMatrix44.getMatrixElement(1,1));
	EXPECT_EQ(7,testMatrix44.getMatrixElement(1,2));
	EXPECT_EQ(8,testMatrix44.getMatrixElement(1,3));

	EXPECT_EQ(9,testMatrix44.getMatrixElement(2,0));
	EXPECT_EQ(10,testMatrix44.getMatrixElement(2,1));
	EXPECT_EQ(11,testMatrix44.getMatrixElement(2,2));
	EXPECT_EQ(12,testMatrix44.getMatrixElement(2,3));

	EXPECT_EQ(13,testMatrix44.getMatrixElement(3,0));
	EXPECT_EQ(14,testMatrix44.getMatrixElement(3,1));
	EXPECT_EQ(15,testMatrix44.getMatrixElement(3,2));
	EXPECT_EQ(16,testMatrix44.getMatrixElement(3,3));

	EXPECT_EQ(result,testMatrix44.swapDiagElement(1,3));
	EXPECT_EQ(1,testMatrix44.getMatrixElement(0,0));
	EXPECT_EQ(2,testMatrix44.getMatrixElement(0,1));
	EXPECT_EQ(3,testMatrix44.getMatrixElement(0,2));
	EXPECT_EQ(4,testMatrix44.getMatrixElement(0,3));

	EXPECT_EQ(5,testMatrix44.getMatrixElement(1,0));
	EXPECT_EQ(16,testMatrix44.getMatrixElement(1,1));
	EXPECT_EQ(7,testMatrix44.getMatrixElement(1,2));
	EXPECT_EQ(8,testMatrix44.getMatrixElement(1,3));

	EXPECT_EQ(9,testMatrix44.getMatrixElement(2,0));
	EXPECT_EQ(10,testMatrix44.getMatrixElement(2,1));
	EXPECT_EQ(11,testMatrix44.getMatrixElement(2,2));
	EXPECT_EQ(12,testMatrix44.getMatrixElement(2,3));

	EXPECT_EQ(13,testMatrix44.getMatrixElement(3,0));
	EXPECT_EQ(14,testMatrix44.getMatrixElement(3,1));
	EXPECT_EQ(15,testMatrix44.getMatrixElement(3,2));
	EXPECT_EQ(6,testMatrix44.getMatrixElement(3,3));

	EXPECT_EQ(result,testMatrix44.swapDiagElement(3,1));
	EXPECT_EQ(1,testMatrix44.getMatrixElement(0,0));
	EXPECT_EQ(2,testMatrix44.getMatrixElement(0,1));
	EXPECT_EQ(3,testMatrix44.getMatrixElement(0,2));
	EXPECT_EQ(4,testMatrix44.getMatrixElement(0,3));

	EXPECT_EQ(5,testMatrix44.getMatrixElement(1,0));
	EXPECT_EQ(6,testMatrix44.getMatrixElement(1,1));
	EXPECT_EQ(7,testMatrix44.getMatrixElement(1,2));
	EXPECT_EQ(8,testMatrix44.getMatrixElement(1,3));

	EXPECT_EQ(9,testMatrix44.getMatrixElement(2,0));
	EXPECT_EQ(10,testMatrix44.getMatrixElement(2,1));
	EXPECT_EQ(11,testMatrix44.getMatrixElement(2,2));
	EXPECT_EQ(12,testMatrix44.getMatrixElement(2,3));

	EXPECT_EQ(13,testMatrix44.getMatrixElement(3,0));
	EXPECT_EQ(14,testMatrix44.getMatrixElement(3,1));
	EXPECT_EQ(15,testMatrix44.getMatrixElement(3,2));
	EXPECT_EQ(16,testMatrix44.getMatrixElement(3,3));
}

/*
 * 测试Matrix Reset 单位矩阵
 */
TEST(BigMatrixResetToE1Test, postive)
{
	BigStaticMatrix testMatrix_23 = BigStaticMatrix(2,3);
	bool result = false;
	EXPECT_EQ(result,testMatrix_23.resetMatrixToI());

	EXPECT_EQ(1,testMatrix_23.getMatrixElement(0,0));
	EXPECT_EQ(2,testMatrix_23.getMatrixElement(0,1));
	EXPECT_EQ(3,testMatrix_23.getMatrixElement(0,2));

	EXPECT_EQ(4,testMatrix_23.getMatrixElement(1,0));
	EXPECT_EQ(5,testMatrix_23.getMatrixElement(1,1));
	EXPECT_EQ(6,testMatrix_23.getMatrixElement(1,2));


	BigStaticMatrix testMatrix_33 = BigStaticMatrix(3,3);
	EXPECT_EQ(1,testMatrix_33.getMatrixElement(0,0));
	EXPECT_EQ(2,testMatrix_33.getMatrixElement(0,1));
	EXPECT_EQ(3,testMatrix_33.getMatrixElement(0,2));

	EXPECT_EQ(4,testMatrix_33.getMatrixElement(1,0));
	EXPECT_EQ(5,testMatrix_33.getMatrixElement(1,1));
	EXPECT_EQ(6,testMatrix_33.getMatrixElement(1,2));

	EXPECT_EQ(7,testMatrix_33.getMatrixElement(2,0));
	EXPECT_EQ(8,testMatrix_33.getMatrixElement(2,1));
	EXPECT_EQ(9,testMatrix_33.getMatrixElement(2,2));

	result = true;
	EXPECT_EQ(result,testMatrix_33.resetMatrixToI());

	EXPECT_EQ(1,testMatrix_33.getMatrixElement(0,0));
	EXPECT_EQ(0,testMatrix_33.getMatrixElement(0,1));
	EXPECT_EQ(0,testMatrix_33.getMatrixElement(0,2));

	EXPECT_EQ(0,testMatrix_33.getMatrixElement(1,0));
	EXPECT_EQ(1,testMatrix_33.getMatrixElement(1,1));
	EXPECT_EQ(0,testMatrix_33.getMatrixElement(1,2));

	EXPECT_EQ(0,testMatrix_33.getMatrixElement(2,0));
	EXPECT_EQ(0,testMatrix_33.getMatrixElement(2,1));
	EXPECT_EQ(1,testMatrix_33.getMatrixElement(2,2));
}

/*
 * 测试Matrix Reset 0矩阵
 */
TEST(BigMatrixResetToZeroTest, postive)
{
	BigStaticMatrix testMatrix_23 = BigStaticMatrix(2,3);

	testMatrix_23.resetMatrixToZero();

	EXPECT_EQ(0,testMatrix_23.getMatrixElement(0,0));
	EXPECT_EQ(0,testMatrix_23.getMatrixElement(0,1));
	EXPECT_EQ(0,testMatrix_23.getMatrixElement(0,2));

	EXPECT_EQ(0,testMatrix_23.getMatrixElement(1,0));
	EXPECT_EQ(0,testMatrix_23.getMatrixElement(1,1));
	EXPECT_EQ(0,testMatrix_23.getMatrixElement(1,2));


	BigStaticMatrix testMatrix_33 = BigStaticMatrix(3,3);
	EXPECT_EQ(1,testMatrix_33.getMatrixElement(0,0));
	EXPECT_EQ(2,testMatrix_33.getMatrixElement(0,1));
	EXPECT_EQ(3,testMatrix_33.getMatrixElement(0,2));

	EXPECT_EQ(4,testMatrix_33.getMatrixElement(1,0));
	EXPECT_EQ(5,testMatrix_33.getMatrixElement(1,1));
	EXPECT_EQ(6,testMatrix_33.getMatrixElement(1,2));

	EXPECT_EQ(7,testMatrix_33.getMatrixElement(2,0));
	EXPECT_EQ(8,testMatrix_33.getMatrixElement(2,1));
	EXPECT_EQ(9,testMatrix_33.getMatrixElement(2,2));

	testMatrix_33.resetMatrixToZero();

	EXPECT_EQ(0,testMatrix_33.getMatrixElement(0,0));
	EXPECT_EQ(0,testMatrix_33.getMatrixElement(0,1));
	EXPECT_EQ(0,testMatrix_33.getMatrixElement(0,2));

	EXPECT_EQ(0,testMatrix_33.getMatrixElement(1,0));
	EXPECT_EQ(0,testMatrix_33.getMatrixElement(1,1));
	EXPECT_EQ(0,testMatrix_33.getMatrixElement(1,2));

	EXPECT_EQ(0,testMatrix_33.getMatrixElement(2,0));
	EXPECT_EQ(0,testMatrix_33.getMatrixElement(2,1));
	EXPECT_EQ(0,testMatrix_33.getMatrixElement(2,2));
}


/*
 * 测试Matrix 主对角线子矩阵 向下移动
 */
TEST(BigMatrixDiagonalSubMoveDownTest, postive)
{
	BigStaticMatrix testMatrix33 = BigStaticMatrix(3,3);

	EXPECT_EQ(1,testMatrix33.getMatrixElement(0,0));
	EXPECT_EQ(2,testMatrix33.getMatrixElement(0,1));
	EXPECT_EQ(3,testMatrix33.getMatrixElement(0,2));

	EXPECT_EQ(4,testMatrix33.getMatrixElement(1,0));
	EXPECT_EQ(5,testMatrix33.getMatrixElement(1,1));
	EXPECT_EQ(6,testMatrix33.getMatrixElement(1,2));

	EXPECT_EQ(7,testMatrix33.getMatrixElement(2,0));
	EXPECT_EQ(8,testMatrix33.getMatrixElement(2,1));
	EXPECT_EQ(9,testMatrix33.getMatrixElement(2,2));

	testMatrix33.moveDiagonalSubMatrixDown(0,1,0);

	EXPECT_EQ(1,testMatrix33.getMatrixElement(0,0));
	EXPECT_EQ(2,testMatrix33.getMatrixElement(0,1));
	EXPECT_EQ(3,testMatrix33.getMatrixElement(0,2));

	EXPECT_EQ(4,testMatrix33.getMatrixElement(1,0));
	EXPECT_EQ(5,testMatrix33.getMatrixElement(1,1));
	EXPECT_EQ(6,testMatrix33.getMatrixElement(1,2));

	EXPECT_EQ(7,testMatrix33.getMatrixElement(2,0));
	EXPECT_EQ(8,testMatrix33.getMatrixElement(2,1));
	EXPECT_EQ(9,testMatrix33.getMatrixElement(2,2));

	testMatrix33.moveDiagonalSubMatrixDown(0,2,0);

	EXPECT_EQ(1,testMatrix33.getMatrixElement(0,0));
	EXPECT_EQ(2,testMatrix33.getMatrixElement(0,1));
	EXPECT_EQ(3,testMatrix33.getMatrixElement(0,2));

	EXPECT_EQ(4,testMatrix33.getMatrixElement(1,0));
	EXPECT_EQ(5,testMatrix33.getMatrixElement(1,1));
	EXPECT_EQ(6,testMatrix33.getMatrixElement(1,2));

	EXPECT_EQ(7,testMatrix33.getMatrixElement(2,0));
	EXPECT_EQ(8,testMatrix33.getMatrixElement(2,1));
	EXPECT_EQ(9,testMatrix33.getMatrixElement(2,2));

	testMatrix33.moveDiagonalSubMatrixDown(0,1,1);

	EXPECT_EQ(1,testMatrix33.getMatrixElement(0,0));
	EXPECT_EQ(0,testMatrix33.getMatrixElement(0,1));
	EXPECT_EQ(3,testMatrix33.getMatrixElement(0,2));

	EXPECT_EQ(0,testMatrix33.getMatrixElement(1,0));
	EXPECT_EQ(1,testMatrix33.getMatrixElement(1,1));
	EXPECT_EQ(2,testMatrix33.getMatrixElement(1,2));

	EXPECT_EQ(7,testMatrix33.getMatrixElement(2,0));
	EXPECT_EQ(4,testMatrix33.getMatrixElement(2,1));
	EXPECT_EQ(5,testMatrix33.getMatrixElement(2,2));

	testMatrix33.resetMatrixToI();

	testMatrix33.setMatrixElement(0,0,5);
	testMatrix33.setMatrixElement(0,1,6);

	testMatrix33.setMatrixElement(1,0,7);
	testMatrix33.setMatrixElement(1,1,8);

	EXPECT_EQ(5,testMatrix33.getMatrixElement(0,0));
	EXPECT_EQ(6,testMatrix33.getMatrixElement(0,1));
	EXPECT_EQ(0,testMatrix33.getMatrixElement(0,2));

	EXPECT_EQ(7,testMatrix33.getMatrixElement(1,0));
	EXPECT_EQ(8,testMatrix33.getMatrixElement(1,1));
	EXPECT_EQ(0,testMatrix33.getMatrixElement(1,2));

	EXPECT_EQ(0,testMatrix33.getMatrixElement(2,0));
	EXPECT_EQ(0,testMatrix33.getMatrixElement(2,1));
	EXPECT_EQ(1,testMatrix33.getMatrixElement(2,2));

	testMatrix33.moveDiagonalSubMatrixDown(0,1,1);

	EXPECT_EQ(1,testMatrix33.getMatrixElement(0,0));
	EXPECT_EQ(0,testMatrix33.getMatrixElement(0,1));
	EXPECT_EQ(0,testMatrix33.getMatrixElement(0,2));

	EXPECT_EQ(0,testMatrix33.getMatrixElement(1,0));
	EXPECT_EQ(5,testMatrix33.getMatrixElement(1,1));
	EXPECT_EQ(6,testMatrix33.getMatrixElement(1,2));

	EXPECT_EQ(0,testMatrix33.getMatrixElement(2,0));
	EXPECT_EQ(7,testMatrix33.getMatrixElement(2,1));
	EXPECT_EQ(8,testMatrix33.getMatrixElement(2,2));


	BigStaticMatrix testMatrix44 = BigStaticMatrix(4,4);

	EXPECT_EQ(1,testMatrix44.getMatrixElement(0,0));
	EXPECT_EQ(2,testMatrix44.getMatrixElement(0,1));
	EXPECT_EQ(3,testMatrix44.getMatrixElement(0,2));
	EXPECT_EQ(4,testMatrix44.getMatrixElement(0,3));

	EXPECT_EQ(5,testMatrix44.getMatrixElement(1,0));
	EXPECT_EQ(6,testMatrix44.getMatrixElement(1,1));
	EXPECT_EQ(7,testMatrix44.getMatrixElement(1,2));
	EXPECT_EQ(8,testMatrix44.getMatrixElement(1,3));

	EXPECT_EQ(9,testMatrix44.getMatrixElement(2,0));
	EXPECT_EQ(10,testMatrix44.getMatrixElement(2,1));
	EXPECT_EQ(11,testMatrix44.getMatrixElement(2,2));
	EXPECT_EQ(12,testMatrix44.getMatrixElement(2,3));

	EXPECT_EQ(13,testMatrix44.getMatrixElement(3,0));
	EXPECT_EQ(14,testMatrix44.getMatrixElement(3,1));
	EXPECT_EQ(15,testMatrix44.getMatrixElement(3,2));
	EXPECT_EQ(16,testMatrix44.getMatrixElement(3,3));

	testMatrix44.moveDiagonalSubMatrixDown(0,1,2);

	EXPECT_EQ(1,testMatrix44.getMatrixElement(0,0));
	EXPECT_EQ(0,testMatrix44.getMatrixElement(0,1));
	EXPECT_EQ(3,testMatrix44.getMatrixElement(0,2));
	EXPECT_EQ(4,testMatrix44.getMatrixElement(0,3));

	EXPECT_EQ(0,testMatrix44.getMatrixElement(1,0));
	EXPECT_EQ(1,testMatrix44.getMatrixElement(1,1));
	EXPECT_EQ(7,testMatrix44.getMatrixElement(1,2));
	EXPECT_EQ(8,testMatrix44.getMatrixElement(1,3));

	EXPECT_EQ(9,testMatrix44.getMatrixElement(2,0));
	EXPECT_EQ(10,testMatrix44.getMatrixElement(2,1));
	EXPECT_EQ(1,testMatrix44.getMatrixElement(2,2));
	EXPECT_EQ(2,testMatrix44.getMatrixElement(2,3));

	EXPECT_EQ(13,testMatrix44.getMatrixElement(3,0));
	EXPECT_EQ(14,testMatrix44.getMatrixElement(3,1));
	EXPECT_EQ(5,testMatrix44.getMatrixElement(3,2));
	EXPECT_EQ(6,testMatrix44.getMatrixElement(3,3));
}

/*
 * 测试Matrix 主对角线子矩阵 向上移动
 */
TEST(BigMatrixDiagonalSubMoveUpTest, postive)
{
	BigStaticMatrix testMatrix33 = BigStaticMatrix(3,3);

	EXPECT_EQ(1,testMatrix33.getMatrixElement(0,0));
	EXPECT_EQ(2,testMatrix33.getMatrixElement(0,1));
	EXPECT_EQ(3,testMatrix33.getMatrixElement(0,2));

	EXPECT_EQ(4,testMatrix33.getMatrixElement(1,0));
	EXPECT_EQ(5,testMatrix33.getMatrixElement(1,1));
	EXPECT_EQ(6,testMatrix33.getMatrixElement(1,2));

	EXPECT_EQ(7,testMatrix33.getMatrixElement(2,0));
	EXPECT_EQ(8,testMatrix33.getMatrixElement(2,1));
	EXPECT_EQ(9,testMatrix33.getMatrixElement(2,2));

	testMatrix33.moveDiagonalSubMatrixUp(1,2,0);
	EXPECT_EQ(1,testMatrix33.getMatrixElement(0,0));
	EXPECT_EQ(2,testMatrix33.getMatrixElement(0,1));
	EXPECT_EQ(3,testMatrix33.getMatrixElement(0,2));

	EXPECT_EQ(4,testMatrix33.getMatrixElement(1,0));
	EXPECT_EQ(5,testMatrix33.getMatrixElement(1,1));
	EXPECT_EQ(6,testMatrix33.getMatrixElement(1,2));

	EXPECT_EQ(7,testMatrix33.getMatrixElement(2,0));
	EXPECT_EQ(8,testMatrix33.getMatrixElement(2,1));
	EXPECT_EQ(9,testMatrix33.getMatrixElement(2,2));

	testMatrix33.moveDiagonalSubMatrixUp(0,2,0);
	EXPECT_EQ(1,testMatrix33.getMatrixElement(0,0));
	EXPECT_EQ(2,testMatrix33.getMatrixElement(0,1));
	EXPECT_EQ(3,testMatrix33.getMatrixElement(0,2));

	EXPECT_EQ(4,testMatrix33.getMatrixElement(1,0));
	EXPECT_EQ(5,testMatrix33.getMatrixElement(1,1));
	EXPECT_EQ(6,testMatrix33.getMatrixElement(1,2));

	EXPECT_EQ(7,testMatrix33.getMatrixElement(2,0));
	EXPECT_EQ(8,testMatrix33.getMatrixElement(2,1));
	EXPECT_EQ(9,testMatrix33.getMatrixElement(2,2));

	testMatrix33.moveDiagonalSubMatrixUp(1,2,1);

	EXPECT_EQ(5,testMatrix33.getMatrixElement(0,0));
	EXPECT_EQ(6,testMatrix33.getMatrixElement(0,1));
	EXPECT_EQ(3,testMatrix33.getMatrixElement(0,2));

	EXPECT_EQ(8,testMatrix33.getMatrixElement(1,0));
	EXPECT_EQ(9,testMatrix33.getMatrixElement(1,1));
	EXPECT_EQ(0,testMatrix33.getMatrixElement(1,2));

	EXPECT_EQ(7,testMatrix33.getMatrixElement(2,0));
	EXPECT_EQ(0,testMatrix33.getMatrixElement(2,1));
	EXPECT_EQ(1,testMatrix33.getMatrixElement(2,2));

	BigStaticMatrix testMatrix44 = BigStaticMatrix(4,4);

	EXPECT_EQ(1,testMatrix44.getMatrixElement(0,0));
	EXPECT_EQ(2,testMatrix44.getMatrixElement(0,1));
	EXPECT_EQ(3,testMatrix44.getMatrixElement(0,2));
	EXPECT_EQ(4,testMatrix44.getMatrixElement(0,3));

	EXPECT_EQ(5,testMatrix44.getMatrixElement(1,0));
	EXPECT_EQ(6,testMatrix44.getMatrixElement(1,1));
	EXPECT_EQ(7,testMatrix44.getMatrixElement(1,2));
	EXPECT_EQ(8,testMatrix44.getMatrixElement(1,3));

	EXPECT_EQ(9,testMatrix44.getMatrixElement(2,0));
	EXPECT_EQ(10,testMatrix44.getMatrixElement(2,1));
	EXPECT_EQ(11,testMatrix44.getMatrixElement(2,2));
	EXPECT_EQ(12,testMatrix44.getMatrixElement(2,3));

	EXPECT_EQ(13,testMatrix44.getMatrixElement(3,0));
	EXPECT_EQ(14,testMatrix44.getMatrixElement(3,1));
	EXPECT_EQ(15,testMatrix44.getMatrixElement(3,2));
	EXPECT_EQ(16,testMatrix44.getMatrixElement(3,3));

	testMatrix44.moveDiagonalSubMatrixUp(2,3,2);

	EXPECT_EQ(11,testMatrix44.getMatrixElement(0,0));
	EXPECT_EQ(12,testMatrix44.getMatrixElement(0,1));
	EXPECT_EQ(3,testMatrix44.getMatrixElement(0,2));
	EXPECT_EQ(4,testMatrix44.getMatrixElement(0,3));

	EXPECT_EQ(15,testMatrix44.getMatrixElement(1,0));
	EXPECT_EQ(16,testMatrix44.getMatrixElement(1,1));
	EXPECT_EQ(7,testMatrix44.getMatrixElement(1,2));
	EXPECT_EQ(8,testMatrix44.getMatrixElement(1,3));

	EXPECT_EQ(9,testMatrix44.getMatrixElement(2,0));
	EXPECT_EQ(10,testMatrix44.getMatrixElement(2,1));
	EXPECT_EQ(1,testMatrix44.getMatrixElement(2,2));
	EXPECT_EQ(0,testMatrix44.getMatrixElement(2,3));

	EXPECT_EQ(13,testMatrix44.getMatrixElement(3,0));
	EXPECT_EQ(14,testMatrix44.getMatrixElement(3,1));
	EXPECT_EQ(0,testMatrix44.getMatrixElement(3,2));
	EXPECT_EQ(1,testMatrix44.getMatrixElement(3,3));
}

