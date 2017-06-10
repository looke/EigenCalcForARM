/*
 * test_StaticVector.cpp
 *
 *  Created on: 2017��6��7��
 *      Author: looke
 */

#include "..\gtest_src\gtest\gtest.h"
#include "StaticVector.h"
#include "math.h"
/*
 * ����vector�����Լ��쳣����
 */
TEST(StaticVectorCreateTest, postive)
{
	EXPECT_NO_THROW(StaticVector(0));
	EXPECT_NO_THROW(StaticVector(10));
	EXPECT_NO_THROW(StaticVector(20));
	//EXPECT_THROW(StaticVector(21), out_of_range);
	//EXPECT_THROW(StaticVector(-1), out_of_range);
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
 * ����vector��������
 */
TEST(StaticVectorSizeResetTest, postive)
{
	StaticVector testVector = StaticVector(3);
	testVector.resetDimension(13);
	EXPECT_EQ(13, testVector.getDimension());

	testVector.resetDimension(5);
	EXPECT_EQ(5, testVector.getDimension());

	testVector.resetDimension(20);
	EXPECT_EQ(20, testVector.getDimension());
}

/*
 * ����vector���������쳣����
 */
TEST(StaticVectorSizeResetExceptionTest, postive)
{
	StaticVector testVector = StaticVector(3);
	EXPECT_NO_THROW(testVector.resetDimension(0));
	EXPECT_NO_THROW(testVector.resetDimension(20));
	EXPECT_EQ(0,testVector.resetDimension(21));
	EXPECT_EQ(0,testVector.resetDimension(-1));
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
TEST(StaticVectorGetElementExceptionTest, negative)
{
	StaticVector testVector = StaticVector(5);
	//�����쳣
	EXPECT_NO_THROW(testVector.getElement(5));
	//�����쳣
	EXPECT_NO_THROW(testVector.getElement(-1));
	//�����쳣
	EXPECT_NO_THROW(testVector.getElement(-3));
}

/*
 * ����vectorԪ������
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
 * ����vectorԪ�������쳣����
 */
TEST(StaticVectorSetElementExceptionTest, postive)
{
	StaticVector testVector = StaticVector(3);
	testVector.setElement(0,23);
	EXPECT_EQ(0, testVector.setElement(3,23));
	EXPECT_EQ(0, testVector.setElement(-1,23));
}

/*
 * ����vector����ģƽ��ֵ
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
 * ����vector����ģ
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
