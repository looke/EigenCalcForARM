/*
 * test_StaticMatrix_RowColumnVector.cpp
 *
 *  Created on: 2017年6月11日
 *      Author: looke
 */
#include "..\gtest_src\gtest\gtest.h"
#include "StaticMatrix.h"
#include "StaticVector.h"
#include "math.h"

/*
 * 测试Matrix 获取行向量
 */
TEST(MatrixGetRowVectorTest, postive)
{
	StaticMatrix testMatrix33 = StaticMatrix(3,3);

	EXPECT_EQ(3,testMatrix33.getRowVector(0)->getDimension());

	EXPECT_EQ(1,testMatrix33.getRowVector(0)->getElement(0));
	EXPECT_EQ(2,testMatrix33.getRowVector(0)->getElement(1));
	EXPECT_EQ(3,testMatrix33.getRowVector(0)->getElement(2));

	EXPECT_EQ(4,testMatrix33.getRowVector(1)->getElement(0));
	EXPECT_EQ(5,testMatrix33.getRowVector(1)->getElement(1));
	EXPECT_EQ(6,testMatrix33.getRowVector(1)->getElement(2));

	EXPECT_EQ(7,testMatrix33.getRowVector(2)->getElement(0));
	EXPECT_EQ(8,testMatrix33.getRowVector(2)->getElement(1));
	EXPECT_EQ(9,testMatrix33.getRowVector(2)->getElement(2));

	StaticMatrix testMatrix24 = StaticMatrix(2,4);
	EXPECT_EQ(4,testMatrix24.getRowVector(0)->getDimension());

	EXPECT_EQ(1,testMatrix24.getRowVector(0)->getElement(0));
	EXPECT_EQ(2,testMatrix24.getRowVector(0)->getElement(1));
	EXPECT_EQ(3,testMatrix24.getRowVector(0)->getElement(2));
	EXPECT_EQ(4,testMatrix24.getRowVector(0)->getElement(3));

	EXPECT_EQ(5,testMatrix24.getRowVector(1)->getElement(0));
	EXPECT_EQ(6,testMatrix24.getRowVector(1)->getElement(1));
	EXPECT_EQ(7,testMatrix24.getRowVector(1)->getElement(2));
	EXPECT_EQ(8,testMatrix24.getRowVector(1)->getElement(3));
}

/*
 * 测试Matrix 获取行向量维度异常
 */
TEST(MatrixGetRowVectorExceptionTest, postive)
{
	StaticMatrix testMatrix33 = StaticMatrix(3,3);

	EXPECT_EQ(0,testMatrix33.getRowVector(-1)->getDimension());
	EXPECT_EQ(0,testMatrix33.getRowVector(3)->getDimension());
	EXPECT_EQ(0,testMatrix33.getRowVector(4)->getDimension());
}

/*
 * 测试Matrix 获取列向量
 */
TEST(MatrixGetColumnVectorTest, postive)
{
	StaticMatrix testMatrix33 = StaticMatrix(3,3);

	EXPECT_EQ(3,testMatrix33.getColumnVector(0)->getDimension());

	EXPECT_EQ(1,testMatrix33.getColumnVector(0)->getElement(0));
	EXPECT_EQ(4,testMatrix33.getColumnVector(0)->getElement(1));
	EXPECT_EQ(7,testMatrix33.getColumnVector(0)->getElement(2));

	EXPECT_EQ(2,testMatrix33.getColumnVector(1)->getElement(0));
	EXPECT_EQ(5,testMatrix33.getColumnVector(1)->getElement(1));
	EXPECT_EQ(8,testMatrix33.getColumnVector(1)->getElement(2));

	EXPECT_EQ(3,testMatrix33.getColumnVector(2)->getElement(0));
	EXPECT_EQ(6,testMatrix33.getColumnVector(2)->getElement(1));
	EXPECT_EQ(9,testMatrix33.getColumnVector(2)->getElement(2));

	StaticMatrix testMatrix24 = StaticMatrix(2,4);
	EXPECT_EQ(2,testMatrix24.getColumnVector(0)->getDimension());

	EXPECT_EQ(1,testMatrix24.getColumnVector(0)->getElement(0));
	EXPECT_EQ(5,testMatrix24.getColumnVector(0)->getElement(1));

	EXPECT_EQ(2,testMatrix24.getColumnVector(1)->getElement(0));
	EXPECT_EQ(6,testMatrix24.getColumnVector(1)->getElement(1));

	EXPECT_EQ(3,testMatrix24.getColumnVector(2)->getElement(0));
	EXPECT_EQ(7,testMatrix24.getColumnVector(2)->getElement(1));

	EXPECT_EQ(4,testMatrix24.getColumnVector(3)->getElement(0));
	EXPECT_EQ(8,testMatrix24.getColumnVector(3)->getElement(1));
}

/*
 * 测试Matrix 获取列向量维度异常
 */
TEST(MatrixGetColumnVectorExceptionTest, postive)
{
	StaticMatrix testMatrix33 = StaticMatrix(3,3);

	EXPECT_EQ(0,testMatrix33.getColumnVector(-1)->getDimension());
	EXPECT_EQ(0,testMatrix33.getColumnVector(3)->getDimension());
	EXPECT_EQ(0,testMatrix33.getColumnVector(4)->getDimension());
}

//获取指定对角子矩阵行向量
TEST(MatrixGetSubMatrixRowVectorTest, postive)
{
	StaticMatrix testMatrix33 = StaticMatrix(3,3);

	EXPECT_EQ(3,testMatrix33.getSubMatrixRowVector(0,0)->getDimension());

	EXPECT_EQ(1,testMatrix33.getSubMatrixRowVector(0,0)->getElement(0));
	EXPECT_EQ(2,testMatrix33.getSubMatrixRowVector(0,0)->getElement(1));
	EXPECT_EQ(3,testMatrix33.getSubMatrixRowVector(0,0)->getElement(2));

	EXPECT_EQ(4,testMatrix33.getSubMatrixRowVector(0,1)->getElement(0));
	EXPECT_EQ(5,testMatrix33.getSubMatrixRowVector(0,1)->getElement(1));
	EXPECT_EQ(6,testMatrix33.getSubMatrixRowVector(0,1)->getElement(2));

	EXPECT_EQ(7,testMatrix33.getSubMatrixRowVector(0,2)->getElement(0));
	EXPECT_EQ(8,testMatrix33.getSubMatrixRowVector(0,2)->getElement(1));
	EXPECT_EQ(9,testMatrix33.getSubMatrixRowVector(0,2)->getElement(2));

	EXPECT_EQ(2,testMatrix33.getSubMatrixRowVector(1,0)->getDimension());

	EXPECT_EQ(5,testMatrix33.getSubMatrixRowVector(1,0)->getElement(0));
	EXPECT_EQ(6,testMatrix33.getSubMatrixRowVector(1,0)->getElement(1));

	EXPECT_EQ(8,testMatrix33.getSubMatrixRowVector(1,1)->getElement(0));
	EXPECT_EQ(9,testMatrix33.getSubMatrixRowVector(1,1)->getElement(1));

	EXPECT_EQ(1,testMatrix33.getSubMatrixRowVector(2,0)->getDimension());

	EXPECT_EQ(9,testMatrix33.getSubMatrixRowVector(2,0)->getElement(0));
}

//获取指定对角子矩阵列向量
TEST(MatrixGetSubMatrixColumnVectorTest, postive)
{
	StaticMatrix testMatrix33 = StaticMatrix(3,3);

	EXPECT_EQ(3,testMatrix33.getSubMatrixColumnVector(0,0)->getDimension());

	EXPECT_EQ(1,testMatrix33.getSubMatrixColumnVector(0,0)->getElement(0));
	EXPECT_EQ(4,testMatrix33.getSubMatrixColumnVector(0,0)->getElement(1));
	EXPECT_EQ(7,testMatrix33.getSubMatrixColumnVector(0,0)->getElement(2));

	EXPECT_EQ(2,testMatrix33.getSubMatrixColumnVector(0,1)->getElement(0));
	EXPECT_EQ(5,testMatrix33.getSubMatrixColumnVector(0,1)->getElement(1));
	EXPECT_EQ(8,testMatrix33.getSubMatrixColumnVector(0,1)->getElement(2));

	EXPECT_EQ(3,testMatrix33.getSubMatrixColumnVector(0,2)->getElement(0));
	EXPECT_EQ(6,testMatrix33.getSubMatrixColumnVector(0,2)->getElement(1));
	EXPECT_EQ(9,testMatrix33.getSubMatrixColumnVector(0,2)->getElement(2));

	EXPECT_EQ(2,testMatrix33.getSubMatrixColumnVector(1,0)->getDimension());

	EXPECT_EQ(5,testMatrix33.getSubMatrixColumnVector(1,0)->getElement(0));
	EXPECT_EQ(8,testMatrix33.getSubMatrixColumnVector(1,0)->getElement(1));

	EXPECT_EQ(6,testMatrix33.getSubMatrixColumnVector(1,1)->getElement(0));
	EXPECT_EQ(9,testMatrix33.getSubMatrixColumnVector(1,1)->getElement(1));

	EXPECT_EQ(1,testMatrix33.getSubMatrixColumnVector(2,0)->getDimension());
	EXPECT_EQ(9,testMatrix33.getSubMatrixColumnVector(2,0)->getElement(0));
}

//获取指定对角子矩阵hessenberg列向量
TEST(MatrixGetSubMatrixColumnHessenbergVectorTest, postive)
{
	StaticMatrix testMatrix33 = StaticMatrix(3,3);
	EXPECT_EQ(2,testMatrix33.getSubMatrixHessenColumnVector(0)->getDimension());

	EXPECT_EQ(4,testMatrix33.getSubMatrixHessenColumnVector(0)->getElement(0));
	EXPECT_EQ(7,testMatrix33.getSubMatrixHessenColumnVector(0)->getElement(1));

	EXPECT_EQ(1,testMatrix33.getSubMatrixHessenColumnVector(1)->getDimension());

	EXPECT_EQ(8,testMatrix33.getSubMatrixHessenColumnVector(1)->getElement(0));
}
