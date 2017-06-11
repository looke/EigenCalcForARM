/*
 * test_HouseHolderTransformation_toEn.cpp
 *
 *  Created on: 2017ƒÍ6‘¬11»’
 *      Author: looke
 */

#include "..\gtest_src\gtest\gtest.h"
#include "StaticVector.h"
#include "StaticMatrix.h"
#include "MatrixMultiplier.h"
#include "HouseholderTransformation.h"
#include "math.h"

/*
 * ≤‚ ‘HouseholderTransformation En ◊Û≥ÀHouseholderæÿ’Û Œ¨∂»“Ï≥£
 */
TEST(HouseholderTransformationLeftMatrixTest3x3_En_Exception, postive)
{
	StaticVector inputVector = StaticVector(3);
	inputVector.setElement(0,1.234);
	inputVector.setElement(1,2);
	inputVector.setElement(2,5);

	StaticMatrix householderMatrix32 = StaticMatrix(3,2);
	StaticMatrix householderMatrix23 = StaticMatrix(2,3);
	StaticMatrix householderMatrix22 = StaticMatrix(2,2);
	StaticMatrix householderMatrix55 = StaticMatrix(5,5);
	HouseholderTransformation testHouseholderTrans = HouseholderTransformation(&inputVector);
	bool result = false;
	EXPECT_EQ(result, testHouseholderTrans.getHouseholderMatrixToEn_ReverseElement(&householderMatrix32));
	EXPECT_EQ(result, testHouseholderTrans.getHouseholderMatrixToEn_ReverseElement(&householderMatrix23));
	EXPECT_EQ(result, testHouseholderTrans.getHouseholderMatrixToEn_ReverseElement(&householderMatrix22));
	EXPECT_EQ(result, testHouseholderTrans.getHouseholderMatrixToEn_ReverseElement(&householderMatrix55));

}
/*
 * ≤‚ ‘HouseholderTransformation En ◊Û≥ÀHouseholderæÿ’Û 3x3
 */
TEST(HouseholderTransformationLeftMatrixTest3x3_pXXXtoEn, postive)
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

	StaticMatrix householderMatrix33 = StaticMatrix(3,3);

	HouseholderTransformation testHouseholderTrans = HouseholderTransformation(&inputVector);
	bool result = true;
	EXPECT_EQ(result, testHouseholderTrans.getHouseholderMatrixToEn_ReverseElement(&householderMatrix33));

	MatrixMultiplier multiVerify = MatrixMultiplier(&householderMatrix33,&vectorMatrix31,&resultMatrix31);
	multiVerify.multiplyCalc();
	resultMatrix31.regularZeroElement();

	EXPECT_EQ(0, resultMatrix31.getMatrixElement(0,0));
	EXPECT_EQ(0, resultMatrix31.getMatrixElement(1,0));
	EXPECT_NE(0, resultMatrix31.getMatrixElement(2,0));
}

/*
 * ≤‚ ‘HouseholderTransformation En ◊Û≥ÀHouseholderæÿ’Û 3x3
 */
TEST(HouseholderTransformationLeftMatrixTest3x3_nXXXtoEn, postive)
{
	StaticVector inputVector = StaticVector(3);
	inputVector.setElement(0,-12);
	inputVector.setElement(1,23);
	inputVector.setElement(2,-13);
	StaticMatrix vectorMatrix31 = StaticMatrix(3,1);
	vectorMatrix31.setMatrixElement(0,0,-12);
	vectorMatrix31.setMatrixElement(1,0,23);
	vectorMatrix31.setMatrixElement(2,0,-13);

	StaticMatrix resultMatrix31 = StaticMatrix(3,1);

	StaticMatrix householderMatrix33 = StaticMatrix(3,3);

	HouseholderTransformation testHouseholderTrans = HouseholderTransformation(&inputVector);
	bool result = true;
	EXPECT_EQ(result, testHouseholderTrans.getHouseholderMatrixToEn_ReverseElement(&householderMatrix33));

	MatrixMultiplier multiVerify = MatrixMultiplier(&householderMatrix33,&vectorMatrix31,&resultMatrix31);
	multiVerify.multiplyCalc();
	resultMatrix31.regularZeroElement();

	EXPECT_EQ(0, resultMatrix31.getMatrixElement(0,0));
	EXPECT_EQ(0, resultMatrix31.getMatrixElement(1,0));
	EXPECT_NE(0, resultMatrix31.getMatrixElement(2,0));
}

/*
 * ≤‚ ‘HouseholderTransformation En ◊Û≥ÀHouseholderæÿ’Û 3x3
 */
TEST(HouseholderTransformationLeftMatrixTest3x3_p0XXtoEn, postive)
{
	StaticVector inputVector = StaticVector(3);
	inputVector.setElement(0,0);
	inputVector.setElement(1,2);
	inputVector.setElement(2,5);
	StaticMatrix vectorMatrix31 = StaticMatrix(3,1);
	vectorMatrix31.setMatrixElement(0,0,0);
	vectorMatrix31.setMatrixElement(1,0,2);
	vectorMatrix31.setMatrixElement(2,0,5);

	StaticMatrix resultMatrix31 = StaticMatrix(3,1);

	StaticMatrix householderMatrix33 = StaticMatrix(3,3);

	HouseholderTransformation testHouseholderTrans = HouseholderTransformation(&inputVector);
	bool result = true;
	EXPECT_EQ(result, testHouseholderTrans.getHouseholderMatrixToEn_ReverseElement(&householderMatrix33));

	MatrixMultiplier multiVerify = MatrixMultiplier(&householderMatrix33,&vectorMatrix31,&resultMatrix31);
	multiVerify.multiplyCalc();
	resultMatrix31.regularZeroElement();

	EXPECT_EQ(0, resultMatrix31.getMatrixElement(0,0));
	EXPECT_EQ(0, resultMatrix31.getMatrixElement(1,0));
	EXPECT_NE(0, resultMatrix31.getMatrixElement(2,0));
}

/*
 * ≤‚ ‘HouseholderTransformation En ◊Û≥ÀHouseholderæÿ’Û 3x3
 */
TEST(HouseholderTransformationLeftMatrixTest3x3_n0XXtoEn, postive)
{
	StaticVector inputVector = StaticVector(3);
	inputVector.setElement(0,0);
	inputVector.setElement(1,-32);
	inputVector.setElement(2,5);
	StaticMatrix vectorMatrix31 = StaticMatrix(3,1);
	vectorMatrix31.setMatrixElement(0,0,0);
	vectorMatrix31.setMatrixElement(1,0,-32);
	vectorMatrix31.setMatrixElement(2,0,5);

	StaticMatrix resultMatrix31 = StaticMatrix(3,1);

	StaticMatrix householderMatrix33 = StaticMatrix(3,3);

	HouseholderTransformation testHouseholderTrans = HouseholderTransformation(&inputVector);
	bool result = true;
	EXPECT_EQ(result, testHouseholderTrans.getHouseholderMatrixToEn_ReverseElement(&householderMatrix33));

	MatrixMultiplier multiVerify = MatrixMultiplier(&householderMatrix33,&vectorMatrix31,&resultMatrix31);
	multiVerify.multiplyCalc();
	resultMatrix31.regularZeroElement();

	EXPECT_EQ(0, resultMatrix31.getMatrixElement(0,0));
	EXPECT_EQ(0, resultMatrix31.getMatrixElement(1,0));
	EXPECT_NE(0, resultMatrix31.getMatrixElement(2,0));
}

/*
 * ≤‚ ‘HouseholderTransformation En ◊Û≥ÀHouseholderæÿ’Û 3x3
 */
TEST(HouseholderTransformationLeftMatrixTest3x3_p00XtoEn, postive)
{
	StaticVector inputVector = StaticVector(3);
	inputVector.setElement(0,0);
	inputVector.setElement(1,0);
	inputVector.setElement(2,5);
	StaticMatrix vectorMatrix31 = StaticMatrix(3,1);
	vectorMatrix31.setMatrixElement(0,0,0);
	vectorMatrix31.setMatrixElement(1,0,0);
	vectorMatrix31.setMatrixElement(2,0,5);

	StaticMatrix resultMatrix31 = StaticMatrix(3,1);

	StaticMatrix householderMatrix33 = StaticMatrix(3,3);

	HouseholderTransformation testHouseholderTrans = HouseholderTransformation(&inputVector);
	bool result = true;
	EXPECT_EQ(result, testHouseholderTrans.getHouseholderMatrixToEn_ReverseElement(&householderMatrix33));

	MatrixMultiplier multiVerify = MatrixMultiplier(&householderMatrix33,&vectorMatrix31,&resultMatrix31);
	multiVerify.multiplyCalc();
	resultMatrix31.regularZeroElement();

	EXPECT_EQ(0, resultMatrix31.getMatrixElement(0,0));
	EXPECT_EQ(0, resultMatrix31.getMatrixElement(1,0));
	EXPECT_EQ(-5, resultMatrix31.getMatrixElement(2,0));

	EXPECT_EQ(1, householderMatrix33.getMatrixElement(0,0));
	EXPECT_EQ(0, householderMatrix33.getMatrixElement(0,1));
	EXPECT_EQ(0, householderMatrix33.getMatrixElement(0,2));

	EXPECT_EQ(0, householderMatrix33.getMatrixElement(1,0));
	EXPECT_EQ(1, householderMatrix33.getMatrixElement(1,1));
	EXPECT_EQ(0, householderMatrix33.getMatrixElement(1,2));

	EXPECT_EQ(0, householderMatrix33.getMatrixElement(2,0));
	EXPECT_EQ(0, householderMatrix33.getMatrixElement(2,1));
	EXPECT_EQ(-1, householderMatrix33.getMatrixElement(2,2));
}

/*
 * ≤‚ ‘HouseholderTransformation En ◊Û≥ÀHouseholderæÿ’Û 3x3
 */
TEST(HouseholderTransformationLeftMatrixTest3x3_n00XtoEn, postive)
{
	StaticVector inputVector = StaticVector(3);
	inputVector.setElement(0,0);
	inputVector.setElement(1,0);
	inputVector.setElement(2,-3.35);
	StaticMatrix vectorMatrix31 = StaticMatrix(3,1);
	vectorMatrix31.setMatrixElement(0,0,0);
	vectorMatrix31.setMatrixElement(1,0,0);
	vectorMatrix31.setMatrixElement(2,0,-3.35);

	StaticMatrix resultMatrix31 = StaticMatrix(3,1);

	StaticMatrix householderMatrix33 = StaticMatrix(3,3);

	HouseholderTransformation testHouseholderTrans = HouseholderTransformation(&inputVector);
	bool result = true;
	EXPECT_EQ(result, testHouseholderTrans.getHouseholderMatrixToEn_ReverseElement(&householderMatrix33));

	MatrixMultiplier multiVerify = MatrixMultiplier(&householderMatrix33,&vectorMatrix31,&resultMatrix31);
	multiVerify.multiplyCalc();
	resultMatrix31.regularZeroElement();

	EXPECT_EQ(0, resultMatrix31.getMatrixElement(0,0));
	EXPECT_EQ(0, resultMatrix31.getMatrixElement(1,0));
	EXPECT_EQ(3.35, resultMatrix31.getMatrixElement(2,0));

	EXPECT_EQ(1, householderMatrix33.getMatrixElement(0,0));
	EXPECT_EQ(0, householderMatrix33.getMatrixElement(0,1));
	EXPECT_EQ(0, householderMatrix33.getMatrixElement(0,2));

	EXPECT_EQ(0, householderMatrix33.getMatrixElement(1,0));
	EXPECT_EQ(1, householderMatrix33.getMatrixElement(1,1));
	EXPECT_EQ(0, householderMatrix33.getMatrixElement(1,2));

	EXPECT_EQ(0, householderMatrix33.getMatrixElement(2,0));
	EXPECT_EQ(0, householderMatrix33.getMatrixElement(2,1));
	EXPECT_EQ(-1, householderMatrix33.getMatrixElement(2,2));
}

/*
 * ≤‚ ‘HouseholderTransformation En ”“≥ÀHouseholderæÿ’Û 3x3
 */
TEST(HouseholderTransformationRightMatrixTest3x3_pXXXtoEn, postive)
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

	StaticMatrix householderMatrix33 = StaticMatrix(3,3);

	HouseholderTransformation testHouseholderTrans = HouseholderTransformation(&inputVector);
	bool result = true;
	EXPECT_EQ(result, testHouseholderTrans.getHouseholderMatrixToEn_ReverseElement(&householderMatrix33));

	MatrixMultiplier multiVerify = MatrixMultiplier(&vectorMatrix13,&householderMatrix33,&resultMatrix13);
	multiVerify.multiplyCalc();
	resultMatrix13.regularZeroElement();

	EXPECT_EQ(0, resultMatrix13.getMatrixElement(0,0));
	EXPECT_EQ(0, resultMatrix13.getMatrixElement(0,1));
	EXPECT_NE(0, resultMatrix13.getMatrixElement(0,2));
}

/*
 * ≤‚ ‘HouseholderTransformation En ”“≥ÀHouseholderæÿ’Û 3x3
 */
TEST(HouseholderTransformationRightMatrixTest3x3_nXXXtoEn, postive)
{
	StaticVector inputVector = StaticVector(3);
	inputVector.setElement(0,-1.234);
	inputVector.setElement(1,2);
	inputVector.setElement(2,-5);
	StaticMatrix vectorMatrix13 = StaticMatrix(1,3);
	vectorMatrix13.setMatrixElement(0,0,-1.234);
	vectorMatrix13.setMatrixElement(0,1,2);
	vectorMatrix13.setMatrixElement(0,2,-5);

	StaticMatrix resultMatrix13 = StaticMatrix(1,3);

	StaticMatrix householderMatrix33 = StaticMatrix(3,3);

	HouseholderTransformation testHouseholderTrans = HouseholderTransformation(&inputVector);
	bool result = true;
	EXPECT_EQ(result, testHouseholderTrans.getHouseholderMatrixToEn_ReverseElement(&householderMatrix33));

	MatrixMultiplier multiVerify = MatrixMultiplier(&vectorMatrix13,&householderMatrix33,&resultMatrix13);
	multiVerify.multiplyCalc();
	resultMatrix13.regularZeroElement();

	EXPECT_EQ(0, resultMatrix13.getMatrixElement(0,0));
	EXPECT_EQ(0, resultMatrix13.getMatrixElement(0,1));
	EXPECT_NE(0, resultMatrix13.getMatrixElement(0,2));
}

/*
 * ≤‚ ‘HouseholderTransformation En ”“≥ÀHouseholderæÿ’Û 3x3
 */
TEST(HouseholderTransformationRightMatrixTest3x3_p0XXtoEn, postive)
{
	StaticVector inputVector = StaticVector(3);
	inputVector.setElement(0,0);
	inputVector.setElement(1,2);
	inputVector.setElement(2,5);
	StaticMatrix vectorMatrix13 = StaticMatrix(1,3);
	vectorMatrix13.setMatrixElement(0,0,0);
	vectorMatrix13.setMatrixElement(0,1,2);
	vectorMatrix13.setMatrixElement(0,2,5);

	StaticMatrix resultMatrix13 = StaticMatrix(1,3);

	StaticMatrix householderMatrix33 = StaticMatrix(3,3);

	HouseholderTransformation testHouseholderTrans = HouseholderTransformation(&inputVector);
	bool result = true;
	EXPECT_EQ(result, testHouseholderTrans.getHouseholderMatrixToEn_ReverseElement(&householderMatrix33));

	MatrixMultiplier multiVerify = MatrixMultiplier(&vectorMatrix13,&householderMatrix33,&resultMatrix13);
	multiVerify.multiplyCalc();
	resultMatrix13.regularZeroElement();

	EXPECT_EQ(0, resultMatrix13.getMatrixElement(0,0));
	EXPECT_EQ(0, resultMatrix13.getMatrixElement(0,1));
	EXPECT_NE(0, resultMatrix13.getMatrixElement(0,2));
}


/*
 * ≤‚ ‘HouseholderTransformation En ”“≥ÀHouseholderæÿ’Û 3x3
 */
TEST(HouseholderTransformationRightMatrixTest3x3_n0XXtoEn, postive)
{
	StaticVector inputVector = StaticVector(3);
	inputVector.setElement(0,0);
	inputVector.setElement(1,2);
	inputVector.setElement(2,-2.5);
	StaticMatrix vectorMatrix13 = StaticMatrix(1,3);
	vectorMatrix13.setMatrixElement(0,0,0);
	vectorMatrix13.setMatrixElement(0,1,2);
	vectorMatrix13.setMatrixElement(0,2,-2.5);

	StaticMatrix resultMatrix13 = StaticMatrix(1,3);

	StaticMatrix householderMatrix33 = StaticMatrix(3,3);

	HouseholderTransformation testHouseholderTrans = HouseholderTransformation(&inputVector);
	bool result = true;
	EXPECT_EQ(result, testHouseholderTrans.getHouseholderMatrixToEn_ReverseElement(&householderMatrix33));

	MatrixMultiplier multiVerify = MatrixMultiplier(&vectorMatrix13,&householderMatrix33,&resultMatrix13);
	multiVerify.multiplyCalc();
	resultMatrix13.regularZeroElement();

	EXPECT_EQ(0, resultMatrix13.getMatrixElement(0,0));
	EXPECT_EQ(0, resultMatrix13.getMatrixElement(0,1));
	EXPECT_NE(0, resultMatrix13.getMatrixElement(0,2));
}

/*
 * ≤‚ ‘HouseholderTransformation En ”“≥ÀHouseholderæÿ’Û 3x3
 */
TEST(HouseholderTransformationRightMatrixTest3x3_p00XtoEn, postive)
{
	StaticVector inputVector = StaticVector(3);
	inputVector.setElement(0,0);
	inputVector.setElement(1,0);
	inputVector.setElement(2,0.15);
	StaticMatrix vectorMatrix13 = StaticMatrix(1,3);
	vectorMatrix13.setMatrixElement(0,0,0);
	vectorMatrix13.setMatrixElement(0,1,0);
	vectorMatrix13.setMatrixElement(0,2,0.15);

	StaticMatrix resultMatrix13 = StaticMatrix(1,3);

	StaticMatrix householderMatrix33 = StaticMatrix(3,3);

	HouseholderTransformation testHouseholderTrans = HouseholderTransformation(&inputVector);
	bool result = true;
	EXPECT_EQ(result, testHouseholderTrans.getHouseholderMatrixToEn_ReverseElement(&householderMatrix33));

	MatrixMultiplier multiVerify = MatrixMultiplier(&vectorMatrix13,&householderMatrix33,&resultMatrix13);
	multiVerify.multiplyCalc();
	resultMatrix13.regularZeroElement();

	EXPECT_EQ(0, resultMatrix13.getMatrixElement(0,0));
	EXPECT_EQ(0, resultMatrix13.getMatrixElement(0,1));
	EXPECT_EQ(-0.15, resultMatrix13.getMatrixElement(0,2));


	EXPECT_EQ(1, householderMatrix33.getMatrixElement(0,0));
	EXPECT_EQ(0, householderMatrix33.getMatrixElement(0,1));
	EXPECT_EQ(0, householderMatrix33.getMatrixElement(0,2));

	EXPECT_EQ(0, householderMatrix33.getMatrixElement(1,0));
	EXPECT_EQ(1, householderMatrix33.getMatrixElement(1,1));
	EXPECT_EQ(0, householderMatrix33.getMatrixElement(1,2));

	EXPECT_EQ(0, householderMatrix33.getMatrixElement(2,0));
	EXPECT_EQ(0, householderMatrix33.getMatrixElement(2,1));
	EXPECT_EQ(-1, householderMatrix33.getMatrixElement(2,2));
}

/*
 * ≤‚ ‘HouseholderTransformation En ”“≥ÀHouseholderæÿ’Û 3x3
 */
TEST(HouseholderTransformationRightMatrixTest3x3_n00XtoEn, postive)
{
	StaticVector inputVector = StaticVector(3);
	inputVector.setElement(0,0);
	inputVector.setElement(1,0);
	inputVector.setElement(2,-30.15);
	StaticMatrix vectorMatrix13 = StaticMatrix(1,3);
	vectorMatrix13.setMatrixElement(0,0,0);
	vectorMatrix13.setMatrixElement(0,1,0);
	vectorMatrix13.setMatrixElement(0,2,-30.15);

	StaticMatrix resultMatrix13 = StaticMatrix(1,3);

	StaticMatrix householderMatrix33 = StaticMatrix(3,3);

	HouseholderTransformation testHouseholderTrans = HouseholderTransformation(&inputVector);
	bool result = true;
	EXPECT_EQ(result, testHouseholderTrans.getHouseholderMatrixToEn_ReverseElement(&householderMatrix33));

	MatrixMultiplier multiVerify = MatrixMultiplier(&vectorMatrix13,&householderMatrix33,&resultMatrix13);
	multiVerify.multiplyCalc();
	resultMatrix13.regularZeroElement();

	EXPECT_EQ(0, resultMatrix13.getMatrixElement(0,0));
	EXPECT_EQ(0, resultMatrix13.getMatrixElement(0,1));
	EXPECT_EQ(30.15, resultMatrix13.getMatrixElement(0,2));


	EXPECT_EQ(1, householderMatrix33.getMatrixElement(0,0));
	EXPECT_EQ(0, householderMatrix33.getMatrixElement(0,1));
	EXPECT_EQ(0, householderMatrix33.getMatrixElement(0,2));

	EXPECT_EQ(0, householderMatrix33.getMatrixElement(1,0));
	EXPECT_EQ(1, householderMatrix33.getMatrixElement(1,1));
	EXPECT_EQ(0, householderMatrix33.getMatrixElement(1,2));

	EXPECT_EQ(0, householderMatrix33.getMatrixElement(2,0));
	EXPECT_EQ(0, householderMatrix33.getMatrixElement(2,1));
	EXPECT_EQ(-1, householderMatrix33.getMatrixElement(2,2));
}

