/*
 * test_factorial.cpp
 *
 *  Created on: 2017Äê6ÔÂ7ÈÕ
 *      Author: looke
 */

#include "..\gtest_src\gtest\gtest.h"
#include "..\include\vector\static\StaticVector.h"

TEST(StaticVectorSizeTest, postive)
{
	StaticVector testVector = StaticVector(3);
	EXPECT_EQ(3, testVector.getDimension());

	StaticVector testVector2 = StaticVector(13);
	EXPECT_EQ(13, testVector2.getDimension());
}

TEST(StaticVectorGetElementTest, postive)
{
	StaticVector testVector = StaticVector(3);
	EXPECT_EQ(1, testVector.getElement(0));
	EXPECT_EQ(1, testVector.getElement(1));
	EXPECT_EQ(1, testVector.getElement(2));

	//
	EXPECT_EQ(1, testVector.getElement(3));
}
