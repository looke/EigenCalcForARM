/*
 * test_GivensTransformation.cpp
 *
 *  Created on: 2017ƒÍ6‘¬10»’
 *      Author: looke
 */
#include "..\gtest_src\gtest\gtest.h"
#include "StaticVector.h"
#include "StaticMatrix.h"
#include "MatrixMultiplier.h"
#include "GivensTransformation.h"
#include "math.h"

/*
 * ≤‚ ‘GivensTransformation  ◊Û≥ÀGivensæÿ’Û…˙≥… 3x3
 */
TEST(GivensTransformationLeftGivensMatrixTest3x3, postive)
{
	StaticVector inputVector = StaticVector(3);
	inputVector.setElement(0,1.234);
	inputVector.setElement(1,2);
	inputVector.setElement(2,5);
	StaticMatrix vectorMatrix31 = StaticMatrix(3,1);
	vectorMatrix31.setMatrixElement(0,0,1.234);
	vectorMatrix31.setMatrixElement(1,0,2);
	vectorMatrix31.setMatrixElement(2,0,5);

	StaticMatrix resultMatrix31 = StaticMatrix(3,1);

	StaticMatrix givensMatrix33 = StaticMatrix(3,3);


	GivensTransformation testGivensTrans = GivensTransformation(&inputVector);
	bool result = true;

	EXPECT_EQ(result, testGivensTrans.getGivensMatrixPreMultiple(0,&givensMatrix33));
	MatrixMultiplier multiVerify = MatrixMultiplier(&givensMatrix33,&vectorMatrix31,&resultMatrix31);
	multiVerify.multiplyCalc();
	resultMatrix31.regularZeroElement();
	EXPECT_EQ(0, resultMatrix31.getMatrixElement(0,0));
	EXPECT_NE(0, resultMatrix31.getMatrixElement(1,0));
	EXPECT_NE(0, resultMatrix31.getMatrixElement(2,0));

	EXPECT_EQ(result, testGivensTrans.getGivensMatrixPreMultiple(1,&givensMatrix33));
	multiVerify.reload(&givensMatrix33,&vectorMatrix31,&resultMatrix31);
	multiVerify.multiplyCalc();
	resultMatrix31.regularZeroElement();
	EXPECT_NE(0, resultMatrix31.getMatrixElement(0,0));
	EXPECT_EQ(0, resultMatrix31.getMatrixElement(1,0));
	EXPECT_NE(0, resultMatrix31.getMatrixElement(2,0));

	EXPECT_EQ(result, testGivensTrans.getGivensMatrixPreMultiple(2,&givensMatrix33));
	multiVerify.reload(&givensMatrix33,&vectorMatrix31,&resultMatrix31);
	multiVerify.multiplyCalc();
	resultMatrix31.regularZeroElement();
	EXPECT_NE(0, resultMatrix31.getMatrixElement(0,0));
	EXPECT_NE(0, resultMatrix31.getMatrixElement(1,0));
	EXPECT_EQ(0, resultMatrix31.getMatrixElement(2,0));
}

/*
 * ≤‚ ‘GivensTransformation  ”“≥ÀGivensæÿ’Û…˙≥… 3x3
 */
TEST(GivensTransformationRightGivensMatrixTest3x3, postive)
{
	StaticVector inputVector = StaticVector(3);
	inputVector.setElement(0,1.234);
	inputVector.setElement(1,2);
	inputVector.setElement(2,5);
	StaticMatrix vectorMatrix13 = StaticMatrix(1,3);
	vectorMatrix13.setMatrixElement(0,0,1.234);
	vectorMatrix13.setMatrixElement(0,1,2);
	vectorMatrix13.setMatrixElement(0,2,5);

	StaticMatrix resultMatrix13 = StaticMatrix(1,3);

	StaticMatrix givensMatrix33 = StaticMatrix(3,3);

	GivensTransformation testGivensTrans = GivensTransformation(&inputVector);

	bool result = true;

	EXPECT_EQ(result, testGivensTrans.getGivensMatrixAfterMultiple(0,&givensMatrix33));
	MatrixMultiplier multiVerify = MatrixMultiplier(&vectorMatrix13,&givensMatrix33,&resultMatrix13);
	multiVerify.multiplyCalc();
	resultMatrix13.regularZeroElement();
	EXPECT_EQ(0, resultMatrix13.getMatrixElement(0,0));
	EXPECT_NE(0, resultMatrix13.getMatrixElement(0,1));
	EXPECT_NE(0, resultMatrix13.getMatrixElement(0,2));

	EXPECT_EQ(result, testGivensTrans.getGivensMatrixAfterMultiple(1,&givensMatrix33));
	multiVerify.reload(&vectorMatrix13,&givensMatrix33,&resultMatrix13);
	multiVerify.multiplyCalc();
	resultMatrix13.regularZeroElement();
	EXPECT_NE(0, resultMatrix13.getMatrixElement(0,0));
	EXPECT_EQ(0, resultMatrix13.getMatrixElement(0,1));
	EXPECT_NE(0, resultMatrix13.getMatrixElement(0,2));

	EXPECT_EQ(result, testGivensTrans.getGivensMatrixAfterMultiple(2,&givensMatrix33));
	multiVerify.reload(&vectorMatrix13,&givensMatrix33,&resultMatrix13);
	multiVerify.multiplyCalc();
	resultMatrix13.regularZeroElement();
	EXPECT_NE(0, resultMatrix13.getMatrixElement(0,0));
	EXPECT_NE(0, resultMatrix13.getMatrixElement(0,1));
	EXPECT_EQ(0, resultMatrix13.getMatrixElement(0,2));
}
