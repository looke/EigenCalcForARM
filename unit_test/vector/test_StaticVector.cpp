/*
 * test_StaticVector.cpp
 *
 *  Created on: 2017年6月7日
 *      Author: looke
 */

#include "..\gtest_src\gtest\gtest.h"
#include "StaticVector.h"
#include "math.h"
/*
 * 测试vector创建以及异常处理
 */
TEST(StaticVectorCreateTest, postive)
{
	EXPECT_NO_THROW(StaticVector(0));
	EXPECT_NO_THROW(StaticVector(10));
}

/*
 * 测试vector长度信息获取
 */
TEST(StaticVectorSizeTest, postive)
{
	StaticVector testVector = StaticVector(3);
	EXPECT_EQ(3, testVector.getDimension());

	StaticVector testVector2 = StaticVector(10);
	EXPECT_EQ(10, testVector2.getDimension());
}

/*
 * 测试vector长度重置
 */
TEST(StaticVectorSizeResetTest, postive)
{
	StaticVector testVector = StaticVector(3);
	testVector.resetDimension(8);
	EXPECT_EQ(8, testVector.getDimension());

	testVector.resetDimension(5);
	EXPECT_EQ(5, testVector.getDimension());
}

/*
 * 测试vector长度重置异常处理
 */
TEST(StaticVectorSizeResetExceptionTest, postive)
{
	StaticVector testVector = StaticVector(3);
	EXPECT_NO_THROW(testVector.resetDimension(0));
	EXPECT_NO_THROW(testVector.resetDimension(10));
	EXPECT_EQ(0,testVector.resetDimension(11));
	EXPECT_EQ(0,testVector.resetDimension(-1));
}

/*
 * 测试vector元素查询
 */
TEST(StaticVectorGetElementTest, postive)
{
	StaticVector testVector = StaticVector(3);
	EXPECT_EQ(1, testVector.getElement(0));
	EXPECT_EQ(1, testVector.getElement(1));
	EXPECT_EQ(1, testVector.getElement(2));
}

/*
 * 测试vector元素查询 异常处理
 */
TEST(StaticVectorGetElementExceptionTest, negative)
{
	StaticVector testVector = StaticVector(5);
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
TEST(StaticVectorSetElementTest, postive)
{
	StaticVector testVector = StaticVector(3);
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
TEST(StaticVectorSetElementExceptionTest, postive)
{
	StaticVector testVector = StaticVector(3);
	testVector.setElement(0,23);
	EXPECT_EQ(0, testVector.setElement(3,23));
	EXPECT_EQ(0, testVector.setElement(-1,23));
}

/*
 * 测试vector计算模平方值
 */
TEST(StaticVectorNormPowerTest, postive)
{
	StaticVector testVector1 = StaticVector(1);
	EXPECT_EQ(1, testVector1.getNormPowerOfVector());

	StaticVector testVector2 = StaticVector(2);
	EXPECT_EQ(2, testVector2.getNormPowerOfVector());

	StaticVector testVector3 = StaticVector(3);
	EXPECT_EQ(3, testVector3.getNormPowerOfVector());

	testVector3.setElement(0,2);
	testVector3.setElement(1,4);
	testVector3.setElement(2,1);
	EXPECT_EQ(21, testVector3.getNormPowerOfVector());
}

/*
 * 测试vector计算模
 */
TEST(StaticVectorNormTest, postive)
{
	StaticVector testVector1 = StaticVector(1);
	EXPECT_EQ(sqrt(1), testVector1.getNormOfVector());

	StaticVector testVector2 = StaticVector(2);
	EXPECT_EQ(sqrt(2), testVector2.getNormOfVector());

	StaticVector testVector3 = StaticVector(3);
	EXPECT_EQ(sqrt(3), testVector3.getNormOfVector());

	testVector3.setElement(0,2);
	testVector3.setElement(1,4);
	testVector3.setElement(2,1);
	EXPECT_EQ(sqrt(21), testVector3.getNormOfVector());
}
