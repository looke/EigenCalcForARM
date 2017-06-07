/*
 * test_StaticVector.cpp
 *
 *  Created on: 2017��6��7��
 *      Author: looke
 */

#include "..\gtest_src\gtest\gtest.h"
#include "..\include\vector\static\StaticVector.h"

/*
 * ����vector�����Լ��쳣����
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
 * ����vector������Ϣ��ȡ
 */
TEST(StaticVectorSizeTest, postive)
{
	StaticVector testVector = StaticVector(3);
	EXPECT_EQ(3, testVector.getDimension());

	StaticVector testVector2 = StaticVector(13);
	EXPECT_EQ(13, testVector2.getDimension());
}

/*
 * ����vectorԪ�ز�ѯ
 */
TEST(StaticVectorGetElementTest, postive)
{
	StaticVector testVector = StaticVector(3);
	EXPECT_EQ(1, testVector.getElement(0));
	EXPECT_EQ(1, testVector.getElement(1));
	EXPECT_EQ(1, testVector.getElement(2));

}

/*
 * ����vectorԪ�ز�ѯ �쳣����
 */
TEST(StaticVectorGetElementExceptionTest, postive)
{
	StaticVector testVector = StaticVector(5);
	//�����쳣
	EXPECT_THROW(testVector.getElement(5), length_error);
	//�����쳣
	EXPECT_THROW(testVector.getElement(-1), length_error);
	//�����쳣
	EXPECT_THROW(testVector.getElement(-3), length_error);
}

/*
 * ����vectorԪ������
 */
TEST(StaticVectorSetElementTest, postive)
{
	StaticVector testVector = StaticVector(3);
	testVector.setElement(0,23);
	EXPECT_EQ(23, testVector.getElement(0));
}

