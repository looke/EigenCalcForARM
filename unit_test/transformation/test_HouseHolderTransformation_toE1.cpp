/*
 * test_HouseHolderTransformation_toE1.cpp
 *
 *  Created on: 2017ƒÍ6‘¬10»’
 *      Author: looke
 */

#include "..\gtest_src\gtest\gtest.h"
#include "StaticVector.h"
#include "StaticMatrix.h"
#include "MatrixMultiplier.h"
#include "HouseholderTransformation.h"
#include "math.h"

/*
 * ≤‚ ‘HouseholderTransformation E1 ◊Û≥ÀHouseholderæÿ’Û 3x3
 */
TEST(HouseholderTransformationLeftMatrixTest3x3_E1_Exception, postive)
{
	StaticVector inputVector = StaticVector(3);
	inputVector.setElement(0,1.234);
	inputVector.setElement(1,2);
	inputVector.setElement(2,5);

	StaticMatrix householderMatrix32 = StaticMatrix(3,2);
	StaticMatrix householderMatrix23 = StaticMatrix(2,3);
	StaticMatrix householderMatrix44 = StaticMatrix(4,4);
	StaticMatrix householderMatrix22 = StaticMatrix(2,2);

	HouseholderTransformation testHouseholderTrans = HouseholderTransformation(&inputVector);
	bool result = false;
	EXPECT_EQ(result, testHouseholderTrans.getHouseholderMatrixToE1_ReverseElement(&householderMatrix32));
	EXPECT_EQ(result, testHouseholderTrans.getHouseholderMatrixToE1_ReverseElement(&householderMatrix23));
	EXPECT_EQ(result, testHouseholderTrans.getHouseholderMatrixToE1_ReverseElement(&householderMatrix44));
	EXPECT_EQ(result, testHouseholderTrans.getHouseholderMatrixToE1_ReverseElement(&householderMatrix22));
}
/*
 * ≤‚ ‘HouseholderTransformation E1 ◊Û≥ÀHouseholderæÿ’Û 3x3
 */
TEST(HouseholderTransformationLeftMatrixTest3x3_pXXXtoE1, postive)
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
	EXPECT_EQ(result, testHouseholderTrans.getHouseholderMatrixToE1_ReverseElement(&householderMatrix33));

	MatrixMultiplier multiVerify = MatrixMultiplier(&householderMatrix33,&vectorMatrix31,&resultMatrix31);
	multiVerify.multiplyCalc();
	resultMatrix31.regularZeroElement();

	EXPECT_NE(0, resultMatrix31.getMatrixElement(0,0));
	EXPECT_EQ(0, resultMatrix31.getMatrixElement(1,0));
	EXPECT_EQ(0, resultMatrix31.getMatrixElement(2,0));
}

/*
 * ≤‚ ‘HouseholderTransformation E1 ◊Û≥ÀHouseholderæÿ’Û 3x3
 */
TEST(HouseholderTransformationLeftMatrixTest3x3_nXXXtoE1, postive)
{
	StaticVector inputVector = StaticVector(3);
	inputVector.setElement(0,-34);
	inputVector.setElement(1,2);
	inputVector.setElement(2,5);
	StaticMatrix vectorMatrix31 = StaticMatrix(3,1);
	vectorMatrix31.setMatrixElement(0,0,-34);
	vectorMatrix31.setMatrixElement(1,0,2);
	vectorMatrix31.setMatrixElement(2,0,5);

	StaticMatrix resultMatrix31 = StaticMatrix(3,1);

	StaticMatrix householderMatrix33 = StaticMatrix(3,3);

	HouseholderTransformation testHouseholderTrans = HouseholderTransformation(&inputVector);
	bool result = true;
	EXPECT_EQ(result, testHouseholderTrans.getHouseholderMatrixToE1_ReverseElement(&householderMatrix33));

	MatrixMultiplier multiVerify = MatrixMultiplier(&householderMatrix33,&vectorMatrix31,&resultMatrix31);
	multiVerify.multiplyCalc();
	resultMatrix31.regularZeroElement();

	EXPECT_NE(0, resultMatrix31.getMatrixElement(0,0));
	EXPECT_EQ(0, resultMatrix31.getMatrixElement(1,0));
	EXPECT_EQ(0, resultMatrix31.getMatrixElement(2,0));
}

/*
 * ≤‚ ‘HouseholderTransformation xx0 to E1 ◊Û≥ÀHouseholderæÿ’Û 3x3
 */
TEST(HouseholderTransformationLeftMatrixTest3x3_pXX0toE1, postive)
{
	StaticVector inputVector = StaticVector(3);
	inputVector.setElement(0,1.234);
	inputVector.setElement(1,2);
	inputVector.setElement(2,0);
	StaticMatrix vectorMatrix31 = StaticMatrix(3,1);
	vectorMatrix31.setMatrixElement(0,0,1.234);
	vectorMatrix31.setMatrixElement(1,0,2);
	vectorMatrix31.setMatrixElement(2,0,0);

	StaticMatrix resultMatrix31 = StaticMatrix(3,1);

	StaticMatrix householderMatrix33 = StaticMatrix(3,3);

	HouseholderTransformation testHouseholderTrans = HouseholderTransformation(&inputVector);
	bool result = true;
	EXPECT_EQ(result, testHouseholderTrans.getHouseholderMatrixToE1_ReverseElement(&householderMatrix33));

	MatrixMultiplier multiVerify = MatrixMultiplier(&householderMatrix33,&vectorMatrix31,&resultMatrix31);
	multiVerify.multiplyCalc();
	resultMatrix31.regularZeroElement();

	EXPECT_NE(0, resultMatrix31.getMatrixElement(0,0));
	EXPECT_EQ(0, resultMatrix31.getMatrixElement(1,0));
	EXPECT_EQ(0, resultMatrix31.getMatrixElement(2,0));
}

/*
 * ≤‚ ‘HouseholderTransformation xx0 to E1 ◊Û≥ÀHouseholderæÿ’Û 3x3
 */
TEST(HouseholderTransformationLeftMatrixTest3x3_nXX0toE1, postive)
{
	StaticVector inputVector = StaticVector(3);
	inputVector.setElement(0,-14);
	inputVector.setElement(1,2);
	inputVector.setElement(2,0);
	StaticMatrix vectorMatrix31 = StaticMatrix(3,1);
	vectorMatrix31.setMatrixElement(0,0,-14);
	vectorMatrix31.setMatrixElement(1,0,2);
	vectorMatrix31.setMatrixElement(2,0,0);

	StaticMatrix resultMatrix31 = StaticMatrix(3,1);

	StaticMatrix householderMatrix33 = StaticMatrix(3,3);

	HouseholderTransformation testHouseholderTrans = HouseholderTransformation(&inputVector);
	bool result = true;
	EXPECT_EQ(result, testHouseholderTrans.getHouseholderMatrixToE1_ReverseElement(&householderMatrix33));

	MatrixMultiplier multiVerify = MatrixMultiplier(&householderMatrix33,&vectorMatrix31,&resultMatrix31);
	multiVerify.multiplyCalc();
	resultMatrix31.regularZeroElement();

	EXPECT_NE(0, resultMatrix31.getMatrixElement(0,0));
	EXPECT_EQ(0, resultMatrix31.getMatrixElement(1,0));
	EXPECT_EQ(0, resultMatrix31.getMatrixElement(2,0));
}

/*
 * ≤‚ ‘HouseholderTransformation pX00 to E1 ◊Û≥ÀHouseholderæÿ’Û 3x3 pXŒ™’˝ ˝
 */
TEST(HouseholderTransformationLeftMatrixTest3x3_pX00toE1, postive)
{
	StaticVector inputVector = StaticVector(3);
	inputVector.setElement(0,1.234);
	inputVector.setElement(1,0);
	inputVector.setElement(2,0);
	StaticMatrix vectorMatrix31 = StaticMatrix(3,1);
	vectorMatrix31.setMatrixElement(0,0,1.234);
	vectorMatrix31.setMatrixElement(1,0,0);
	vectorMatrix31.setMatrixElement(2,0,0);

	StaticMatrix resultMatrix31 = StaticMatrix(3,1);

	StaticMatrix householderMatrix33 = StaticMatrix(3,3);

	HouseholderTransformation testHouseholderTrans = HouseholderTransformation(&inputVector);
	bool result = true;
	EXPECT_EQ(result, testHouseholderTrans.getHouseholderMatrixToE1_ReverseElement(&householderMatrix33));

	MatrixMultiplier multiVerify = MatrixMultiplier(&householderMatrix33,&vectorMatrix31,&resultMatrix31);
	multiVerify.multiplyCalc();
	resultMatrix31.regularZeroElement();

	EXPECT_EQ(-1.234, resultMatrix31.getMatrixElement(0,0));
	EXPECT_EQ(0, resultMatrix31.getMatrixElement(1,0));
	EXPECT_EQ(0, resultMatrix31.getMatrixElement(2,0));


	EXPECT_EQ(-1, householderMatrix33.getMatrixElement(0,0));
	EXPECT_EQ(0, householderMatrix33.getMatrixElement(0,1));
	EXPECT_EQ(0, householderMatrix33.getMatrixElement(0,2));

	EXPECT_EQ(0, householderMatrix33.getMatrixElement(1,0));
	EXPECT_EQ(1, householderMatrix33.getMatrixElement(1,1));
	EXPECT_EQ(0, householderMatrix33.getMatrixElement(1,2));

	EXPECT_EQ(0, householderMatrix33.getMatrixElement(2,0));
	EXPECT_EQ(0, householderMatrix33.getMatrixElement(2,1));
	EXPECT_EQ(1, householderMatrix33.getMatrixElement(2,2));
}

/*
 * ≤‚ ‘HouseholderTransformation nX00 to E1 ◊Û≥ÀHouseholderæÿ’Û 3x3 nXŒ™∏∫ ˝
 */
TEST(HouseholderTransformationLeftMatrixTest3x3_nX00toE1, postive)
{
	StaticVector inputVector = StaticVector(3);
	inputVector.setElement(0,-1.5);
	inputVector.setElement(1,0);
	inputVector.setElement(2,0);
	StaticMatrix vectorMatrix31 = StaticMatrix(3,1);
	vectorMatrix31.setMatrixElement(0,0,-1.5);
	vectorMatrix31.setMatrixElement(1,0,0);
	vectorMatrix31.setMatrixElement(2,0,0);

	StaticMatrix resultMatrix31 = StaticMatrix(3,1);

	StaticMatrix householderMatrix33 = StaticMatrix(3,3);

	HouseholderTransformation testHouseholderTrans = HouseholderTransformation(&inputVector);
	bool result = true;
	EXPECT_EQ(result, testHouseholderTrans.getHouseholderMatrixToE1_ReverseElement(&householderMatrix33));

	MatrixMultiplier multiVerify = MatrixMultiplier(&householderMatrix33,&vectorMatrix31,&resultMatrix31);
	multiVerify.multiplyCalc();
	resultMatrix31.regularZeroElement();

	EXPECT_EQ(1.5, resultMatrix31.getMatrixElement(0,0));
	EXPECT_EQ(0, resultMatrix31.getMatrixElement(1,0));
	EXPECT_EQ(0, resultMatrix31.getMatrixElement(2,0));


	EXPECT_EQ(-1, householderMatrix33.getMatrixElement(0,0));
	EXPECT_EQ(0, householderMatrix33.getMatrixElement(0,1));
	EXPECT_EQ(0, householderMatrix33.getMatrixElement(0,2));

	EXPECT_EQ(0, householderMatrix33.getMatrixElement(1,0));
	EXPECT_EQ(1, householderMatrix33.getMatrixElement(1,1));
	EXPECT_EQ(0, householderMatrix33.getMatrixElement(1,2));

	EXPECT_EQ(0, householderMatrix33.getMatrixElement(2,0));
	EXPECT_EQ(0, householderMatrix33.getMatrixElement(2,1));
	EXPECT_EQ(1, householderMatrix33.getMatrixElement(2,2));
}


/*
 * ≤‚ ‘HouseholderTransformation E1 ”“≥ÀHouseholderæÿ’Û 3x3
 */
TEST(HouseholderTransformationRightMatrixTest3x3_pXXXtoE1, postive)
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
	EXPECT_EQ(result, testHouseholderTrans.getHouseholderMatrixToE1_ReverseElement(&householderMatrix33));

	MatrixMultiplier multiVerify = MatrixMultiplier(&vectorMatrix13,&householderMatrix33,&resultMatrix13);
	multiVerify.multiplyCalc();
	resultMatrix13.regularZeroElement();

	EXPECT_NE(0, resultMatrix13.getMatrixElement(0,0));
	EXPECT_EQ(0, resultMatrix13.getMatrixElement(0,1));
	EXPECT_EQ(0, resultMatrix13.getMatrixElement(0,2));
}

/*
 * ≤‚ ‘HouseholderTransformation E1 ”“≥ÀHouseholderæÿ’Û 3x3
 */
TEST(HouseholderTransformationRightMatrixTest3x3_nXXXtoE1, postive)
{
	StaticVector inputVector = StaticVector(3);
	inputVector.setElement(0,-1.234);
	inputVector.setElement(1,2);
	inputVector.setElement(2,5);
	StaticMatrix vectorMatrix13 = StaticMatrix(1,3);
	vectorMatrix13.setMatrixElement(0,0,-1.234);
	vectorMatrix13.setMatrixElement(0,1,2);
	vectorMatrix13.setMatrixElement(0,2,5);

	StaticMatrix resultMatrix13 = StaticMatrix(1,3);

	StaticMatrix householderMatrix33 = StaticMatrix(3,3);

	HouseholderTransformation testHouseholderTrans = HouseholderTransformation(&inputVector);
	bool result = true;
	EXPECT_EQ(result, testHouseholderTrans.getHouseholderMatrixToE1_ReverseElement(&householderMatrix33));

	MatrixMultiplier multiVerify = MatrixMultiplier(&vectorMatrix13,&householderMatrix33,&resultMatrix13);
	multiVerify.multiplyCalc();
	resultMatrix13.regularZeroElement();

	EXPECT_NE(0, resultMatrix13.getMatrixElement(0,0));
	EXPECT_EQ(0, resultMatrix13.getMatrixElement(0,1));
	EXPECT_EQ(0, resultMatrix13.getMatrixElement(0,2));
}

/*
 * ≤‚ ‘HouseholderTransformation E1 ”“≥ÀHouseholderæÿ’Û 3x3
 */
TEST(HouseholderTransformationRightMatrixTest3x3_pXX0toE1, postive)
{
	StaticVector inputVector = StaticVector(3);
	inputVector.setElement(0,1.234);
	inputVector.setElement(1,2);
	inputVector.setElement(2,0);
	StaticMatrix vectorMatrix13 = StaticMatrix(1,3);
	vectorMatrix13.setMatrixElement(0,0,1.234);
	vectorMatrix13.setMatrixElement(0,1,2);
	vectorMatrix13.setMatrixElement(0,2,0);

	StaticMatrix resultMatrix13 = StaticMatrix(1,3);

	StaticMatrix householderMatrix33 = StaticMatrix(3,3);

	HouseholderTransformation testHouseholderTrans = HouseholderTransformation(&inputVector);
	bool result = true;
	EXPECT_EQ(result, testHouseholderTrans.getHouseholderMatrixToE1_ReverseElement(&householderMatrix33));

	MatrixMultiplier multiVerify = MatrixMultiplier(&vectorMatrix13,&householderMatrix33,&resultMatrix13);
	multiVerify.multiplyCalc();
	resultMatrix13.regularZeroElement();

	EXPECT_NE(0, resultMatrix13.getMatrixElement(0,0));
	EXPECT_EQ(0, resultMatrix13.getMatrixElement(0,1));
	EXPECT_EQ(0, resultMatrix13.getMatrixElement(0,2));
}


/*
 * ≤‚ ‘HouseholderTransformation E1 ”“≥ÀHouseholderæÿ’Û 3x3
 */
TEST(HouseholderTransformationRightMatrixTest3x3_nXX0toE1, postive)
{
	StaticVector inputVector = StaticVector(3);
	inputVector.setElement(0,-11);
	inputVector.setElement(1,2);
	inputVector.setElement(2,0);
	StaticMatrix vectorMatrix13 = StaticMatrix(1,3);
	vectorMatrix13.setMatrixElement(0,0,-11);
	vectorMatrix13.setMatrixElement(0,1,2);
	vectorMatrix13.setMatrixElement(0,2,0);

	StaticMatrix resultMatrix13 = StaticMatrix(1,3);

	StaticMatrix householderMatrix33 = StaticMatrix(3,3);

	HouseholderTransformation testHouseholderTrans = HouseholderTransformation(&inputVector);
	bool result = true;
	EXPECT_EQ(result, testHouseholderTrans.getHouseholderMatrixToE1_ReverseElement(&householderMatrix33));

	MatrixMultiplier multiVerify = MatrixMultiplier(&vectorMatrix13,&householderMatrix33,&resultMatrix13);
	multiVerify.multiplyCalc();
	resultMatrix13.regularZeroElement();

	EXPECT_NE(0, resultMatrix13.getMatrixElement(0,0));
	EXPECT_EQ(0, resultMatrix13.getMatrixElement(0,1));
	EXPECT_EQ(0, resultMatrix13.getMatrixElement(0,2));
}

/*
 * ≤‚ ‘HouseholderTransformation E1 ”“≥ÀHouseholderæÿ’Û 3x3
 */
TEST(HouseholderTransformationRightMatrixTest3x3_pX00toE1, postive)
{
	StaticVector inputVector = StaticVector(3);
	inputVector.setElement(0,1.234);
	inputVector.setElement(1,0);
	inputVector.setElement(2,0);
	StaticMatrix vectorMatrix13 = StaticMatrix(1,3);
	vectorMatrix13.setMatrixElement(0,0,1.234);
	vectorMatrix13.setMatrixElement(0,1,0);
	vectorMatrix13.setMatrixElement(0,2,0);

	StaticMatrix resultMatrix13 = StaticMatrix(1,3);

	StaticMatrix householderMatrix33 = StaticMatrix(3,3);

	HouseholderTransformation testHouseholderTrans = HouseholderTransformation(&inputVector);
	bool result = true;
	EXPECT_EQ(result, testHouseholderTrans.getHouseholderMatrixToE1_ReverseElement(&householderMatrix33));

	MatrixMultiplier multiVerify = MatrixMultiplier(&vectorMatrix13,&householderMatrix33,&resultMatrix13);
	multiVerify.multiplyCalc();
	resultMatrix13.regularZeroElement();

	EXPECT_EQ(-1.234, resultMatrix13.getMatrixElement(0,0));
	EXPECT_EQ(0, resultMatrix13.getMatrixElement(0,1));
	EXPECT_EQ(0, resultMatrix13.getMatrixElement(0,2));

	EXPECT_EQ(-1, householderMatrix33.getMatrixElement(0,0));
	EXPECT_EQ(0, householderMatrix33.getMatrixElement(0,1));
	EXPECT_EQ(0, householderMatrix33.getMatrixElement(0,2));

	EXPECT_EQ(0, householderMatrix33.getMatrixElement(1,0));
	EXPECT_EQ(1, householderMatrix33.getMatrixElement(1,1));
	EXPECT_EQ(0, householderMatrix33.getMatrixElement(1,2));

	EXPECT_EQ(0, householderMatrix33.getMatrixElement(2,0));
	EXPECT_EQ(0, householderMatrix33.getMatrixElement(2,1));
	EXPECT_EQ(1, householderMatrix33.getMatrixElement(2,2));
}

/*
 * ≤‚ ‘HouseholderTransformation E1 ”“≥ÀHouseholderæÿ’Û 3x3
 */
TEST(HouseholderTransformationRightMatrixTest3x3_nX00toE1, postive)
{
	StaticVector inputVector = StaticVector(3);
	inputVector.setElement(0,-1.234);
	inputVector.setElement(1,0);
	inputVector.setElement(2,0);
	StaticMatrix vectorMatrix13 = StaticMatrix(1,3);
	vectorMatrix13.setMatrixElement(0,0,-1.234);
	vectorMatrix13.setMatrixElement(0,1,0);
	vectorMatrix13.setMatrixElement(0,2,0);

	StaticMatrix resultMatrix13 = StaticMatrix(1,3);

	StaticMatrix householderMatrix33 = StaticMatrix(3,3);

	HouseholderTransformation testHouseholderTrans = HouseholderTransformation(&inputVector);
	bool result = true;
	EXPECT_EQ(result, testHouseholderTrans.getHouseholderMatrixToE1_ReverseElement(&householderMatrix33));

	MatrixMultiplier multiVerify = MatrixMultiplier(&vectorMatrix13,&householderMatrix33,&resultMatrix13);
	multiVerify.multiplyCalc();
	resultMatrix13.regularZeroElement();

	EXPECT_EQ(1.234, resultMatrix13.getMatrixElement(0,0));
	EXPECT_EQ(0, resultMatrix13.getMatrixElement(0,1));
	EXPECT_EQ(0, resultMatrix13.getMatrixElement(0,2));

	EXPECT_EQ(-1, householderMatrix33.getMatrixElement(0,0));
	EXPECT_EQ(0, householderMatrix33.getMatrixElement(0,1));
	EXPECT_EQ(0, householderMatrix33.getMatrixElement(0,2));

	EXPECT_EQ(0, householderMatrix33.getMatrixElement(1,0));
	EXPECT_EQ(1, householderMatrix33.getMatrixElement(1,1));
	EXPECT_EQ(0, householderMatrix33.getMatrixElement(1,2));

	EXPECT_EQ(0, householderMatrix33.getMatrixElement(2,0));
	EXPECT_EQ(0, householderMatrix33.getMatrixElement(2,1));
	EXPECT_EQ(1, householderMatrix33.getMatrixElement(2,2));
}
