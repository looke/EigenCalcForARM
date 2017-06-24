/*
 * test_BigStaticVector.cpp
 *
 *  Created on: 2017��6��24��
 *      Author: looke
 */

#include "..\gtest_src\gtest\gtest.h"
#include "BigStaticVector.h"
#include "math.h"

/*
 * ����vector�����Լ��쳣����
 */
TEST(BigStaticVectorCreateTest, postive)
{
	EXPECT_NO_THROW(BigStaticVector(0));
	EXPECT_NO_THROW(BigStaticVector(10));
	EXPECT_NO_THROW(BigStaticVector(20));
}

/*
 * ����vector������Ϣ��ȡ
 */
TEST(BigStaticVectorSizeTest, postive)
{
	BigStaticVector testVector = BigStaticVector(3);
	EXPECT_EQ(3, testVector.getDimension());

	BigStaticVector testVector2 = BigStaticVector(13);
	EXPECT_EQ(13, testVector2.getDimension());
}

/*
 * ����vector��������
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
 * ����vector���������쳣����
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
 * ����vectorԪ�ز�ѯ
 */
TEST(BigStaticVectorGetElementTest, postive)
{
	BigStaticVector testVector = BigStaticVector(3);
	EXPECT_EQ(1, testVector.getElement(0));
	EXPECT_EQ(1, testVector.getElement(1));
	EXPECT_EQ(1, testVector.getElement(2));
}

/*
 * ����vectorԪ�ز�ѯ �쳣����
 */
TEST(BigStaticVectorGetElementExceptionTest, negative)
{
	BigStaticVector testVector = BigStaticVector(5);
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
 * ����vectorԪ�������쳣����
 */
TEST(BigStaticVectorSetElementExceptionTest, postive)
{
	BigStaticVector testVector = BigStaticVector(3);
	testVector.setElement(0,23);
	EXPECT_EQ(0, testVector.setElement(3,23));
	EXPECT_EQ(0, testVector.setElement(-1,23));
}

/*
 * ����vector����ģƽ��ֵ
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
 * ����vector����ģ
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

