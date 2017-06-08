/*
 * test_StaticVector_afterReDimension.cpp
 *
 *  Created on: 2017年6月8日
 *      Author: looke
 */

#include "..\gtest_src\gtest\gtest.h"
#include "..\include\vector\static\StaticVector.h"
#include "math.h"


class AfterReDimentionTest:public testing::Test
{
protected:
    void SetUp()
    {
    	m_Svector = StaticVector(20);
    	m_Svector.setElement(12,3);
    	m_Svector.setElement(0,11);
    	m_Svector.setElement(1,13);


    	m_Svector.resetDimension(3);

    }
    void TearDown()
    {

    }

    StaticVector m_Svector;
};

/*
 * 测试vector长度信息获取
 */
TEST_F(AfterReDimentionTest,StaticVectorSizeTest)
{
	EXPECT_EQ(3, m_Svector.getDimension());
}

/*
 * 测试vector元素查询 异常处理
 */
TEST_F(AfterReDimentionTest,StaticVectorGetElementExceptionTest)
{
	EXPECT_THROW(m_Svector.getElement(5), out_of_range);
	EXPECT_THROW(m_Svector.getElement(-1), out_of_range);
	EXPECT_THROW(m_Svector.getElement(3), out_of_range);
}

/*
 * 测试vector元素设置异常处理
 */
TEST_F(AfterReDimentionTest,StaticVectorSetElementExceptionTest)
{
	EXPECT_THROW(m_Svector.setElement(3,23), out_of_range);
	EXPECT_THROW(m_Svector.setElement(-1,23), out_of_range);
}

