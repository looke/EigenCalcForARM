/*
 * test_factorial.cpp
 *
 *  Created on: 2017年6月7日
 *      Author: looke
 */

#include "..\gtest_src\gtest\gtest.h"
#include "..\include\vector\static\StaticVector.h"

/*
 * 测试vector创建以及异常处理
 */
TEST(StaticVectorCreateTest, postive)
{
	EXPECT_NO_THROW(StaticVector(0));
	EXPECT_NO_THROW(StaticVector(10));
	EXPECT_NO_THROW(StaticVector(20));
	EXPECT_THROW(StaticVector(21), length_error);
	EXPECT_THROW(StaticVector(-1), length_error);
}

/*
 * 测试vector长度信息获取
 */
TEST(StaticVectorSizeTest, postive)
{
	StaticVector testVector = StaticVector(3);
	EXPECT_EQ(3, testVector.getDimension());

	StaticVector testVector2 = StaticVector(13);
	EXPECT_EQ(13, testVector2.getDimension());
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
TEST(StaticVectorGetElementExceptionTest, postive)
{
	StaticVector testVector = StaticVector(5);
	//测试异常
	EXPECT_THROW(testVector.getElement(5), length_error);
	//测试异常
	EXPECT_THROW(testVector.getElement(-1), length_error);
	//测试异常
	EXPECT_THROW(testVector.getElement(-3), length_error);
}

