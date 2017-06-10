/*
 * test_StaticVector_afterReDimension.cpp
 *
 *  Created on: 2017��6��8��
 *      Author: looke
 */

#include "..\gtest_src\gtest\gtest.h"
#include "StaticVector.h"
#include "math.h"
#include <stdexcept>

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
 * ����vector������Ϣ��ȡ
 */
TEST_F(AfterReDimentionTest,StaticVectorSizeTest)
{
	EXPECT_EQ(3, m_Svector.getDimension());
}

/*
 * ����vectorԪ�ز�ѯ �쳣����
 */

TEST_F(AfterReDimentionTest,StaticVectorGetElementExceptionTest)
{
	//EXPECT_THROW(m_Svector.getElement(25),std::logic_error);
	//EXPECT_THROW(m_Svector.getElement(-1),std::logic_error);
	//EXPECT_THROW(m_Svector.getElement(3),std::logic_error);
}

/*
 * ����vectorԪ�������쳣����
 */
TEST_F(AfterReDimentionTest,StaticVectorSetElementExceptionTest)
{
	EXPECT_EQ(0, m_Svector.setElement(3,23));
	EXPECT_EQ(0, m_Svector.setElement(-1,23));
}

