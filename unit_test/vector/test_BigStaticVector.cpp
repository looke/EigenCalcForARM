/*
 * test_BigStaticVector.cpp
 *
 *  Created on: 2017年6月24日
 *      Author: looke
 */

#include "..\gtest_src\gtest\gtest.h"
#include "BigStaticVector.h"
#include "math.h"

/*
 * 测试vector创建以及异常处理
 */
TEST(BigStaticVectorCreateTest, postive)
{
	EXPECT_NO_THROW(BigStaticVector(0));
	EXPECT_NO_THROW(BigStaticVector(10));
	EXPECT_NO_THROW(BigStaticVector(20));
}

/*
 * 测试vector长度信息获取
 */
TEST(BigStaticVectorSizeTest, postive)
{
	BigStaticVector testVector = BigStaticVector(3);
	EXPECT_EQ(3, testVector.getDimension());

	BigStaticVector testVector2 = BigStaticVector(13);
	EXPECT_EQ(13, testVector2.getDimension());
}

/*
 * 测试vector长度重置
 */
TEST(BigStaticVectorSizeResetTest, postive)
{
	BigStaticVector testVector = BigStaticVector(3);
	testVector.resetDimension(13);
	EXPECT_EQ(13, testVector.getDimension());

	testVector.resetDimension(5);
	EXPECT_EQ(5, testVector.getDimension());

	testVector.resetDimension(20);
	EXPECT_EQ(20, testVector.getDimension());
}

/*
 * 测试vector长度重置异常处理
 */
TEST(BigStaticVectorSizeResetExceptionTest, postive)
{
	BigStaticVector testVector = BigStaticVector(3);
	EXPECT_NO_THROW(testVector.resetDimension(0));
	EXPECT_NO_THROW(testVector.resetDimension(20));
	EXPECT_EQ(0,testVector.resetDimension(21));
	EXPECT_EQ(0,testVector.resetDimension(-1));
}

/*
 * 测试vector元素查询
 */
TEST(BigStaticVectorGetElementTest, postive)
{
	BigStaticVector testVector = BigStaticVector(3);
	EXPECT_EQ(1, testVector.getElement(0));
	EXPECT_EQ(1, testVector.getElement(1));
	EXPECT_EQ(1, testVector.getElement(2));
}

/*
 * 测试vector元素查询 异常处理
 */
TEST(BigStaticVectorGetElementExceptionTest, negative)
{
	BigStaticVector testVector = BigStaticVector(5);
	//测试异常
	EXPECT_NO_THROW(testVector.getElement(5));
	//测试异常
	EXPECT_NO_THROW(testVector.getElement(-1));
	//测试异常
	EXPECT_NO_THROW(testVector.getElement(-3));
}

/*
 * 测试vector元素设置
 */
TEST(BigStaticVectorSetElementTest, postive)
{
	BigStaticVector testVector = BigStaticVector(3);
	testVector.setElement(0,23);
	testVector.setElement(1,33);
	testVector.setElement(2,13);
	EXPECT_EQ(23, testVector.getElement(0));
	EXPECT_EQ(33, testVector.getElement(1));
	EXPECT_EQ(13, testVector.getElement(2));
}

/*
 * 测试vector元素设置异常处理
 */
TEST(BigStaticVectorSetElementExceptionTest, postive)
{
	BigStaticVector testVector = BigStaticVector(3);
	testVector.setElement(0,23);
	EXPECT_EQ(0, testVector.setElement(3,23));
	EXPECT_EQ(0, testVector.setElement(-1,23));
}

/*
 * 测试vector计算模平方值
 */
TEST(BigStaticVectorNormPowerTest, postive)
{
	BigStaticVector testVector1 = BigStaticVector(1);
	EXPECT_EQ(1, testVector1.getNormPowerOfVector());

	BigStaticVector testVector2 = BigStaticVector(2);
	EXPECT_EQ(2, testVector2.getNormPowerOfVector());

	BigStaticVector testVector3 = BigStaticVector(3);
	EXPECT_EQ(3, testVector3.getNormPowerOfVector());

	testVector3.setElement(0,2);
	testVector3.setElement(1,4);
	testVector3.setElement(2,1);
	EXPECT_EQ(21, testVector3.getNormPowerOfVector());
}

/*
 * 测试vector计算模
 */
TEST(BigStaticVectorNormTest, postive)
{
	BigStaticVector testVector1 = BigStaticVector(1);
	EXPECT_EQ(sqrt(1), testVector1.getNormOfVector());

	BigStaticVector testVector2 = BigStaticVector(2);
	EXPECT_EQ(sqrt(2), testVector2.getNormOfVector());

	BigStaticVector testVector3 = BigStaticVector(3);
	EXPECT_EQ(sqrt(3), testVector3.getNormOfVector());

	testVector3.setElement(0,2);
	testVector3.setElement(1,4);
	testVector3.setElement(2,1);
	EXPECT_EQ(sqrt(21), testVector3.getNormOfVector());
}

