/*
 * test_GivensTransformation.cpp
 *
 *  Created on: 2017年6月10日
 *      Author: looke
 */
#include "..\gtest_src\gtest\gtest.h"
#include "StaticVector.h"
#include "StaticMatrix.h"
#include "MatrixMultiplier.h"
#include "GivensTransformation.h"
#include "math.h"

/*
 * 测试GivensTransformation  左乘Givens矩阵生成 3x3
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
 * 测试GivensTransformation  右乘Givens矩阵生成 3x3
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


/*
 * 测试GivensTransformation Givens矩阵维度不匹配
 */
TEST(GivensTransformationExceptionTest, postive)
{
	StaticVector inputVector = StaticVector(3);
	StaticMatrix givensMatrix33 = StaticMatrix(3,3);
	StaticMatrix givensMatrix42 = StaticMatrix(4,2);

	StaticMatrix givensMatrix32 = StaticMatrix(3,2);
	StaticMatrix givensMatrix53 = StaticMatrix(5,3);

	GivensTransformation testGivensTrans = GivensTransformation(&inputVector);

	bool result = false;

	EXPECT_EQ(result, testGivensTrans.getGivensMatrixAfterMultiple(0,&givensMatrix42));
	EXPECT_EQ(result, testGivensTrans.getGivensMatrixAfterMultiple(1,&givensMatrix42));
	EXPECT_EQ(result, testGivensTrans.getGivensMatrixAfterMultiple(2,&givensMatrix42));

	EXPECT_EQ(result, testGivensTrans.getGivensMatrixAfterMultiple(0,&givensMatrix32));
	EXPECT_EQ(result, testGivensTrans.getGivensMatrixAfterMultiple(1,&givensMatrix32));
	EXPECT_EQ(result, testGivensTrans.getGivensMatrixAfterMultiple(2,&givensMatrix32));

	EXPECT_EQ(result, testGivensTrans.getGivensMatrixAfterMultiple(0,&givensMatrix53));
	EXPECT_EQ(result, testGivensTrans.getGivensMatrixAfterMultiple(1,&givensMatrix53));
	EXPECT_EQ(result, testGivensTrans.getGivensMatrixAfterMultiple(2,&givensMatrix53));

	EXPECT_EQ(result, testGivensTrans.getGivensMatrixAfterMultiple(3,&givensMatrix33));
	EXPECT_EQ(result, testGivensTrans.getGivensMatrixAfterMultiple(4,&givensMatrix33));

	EXPECT_EQ(result, testGivensTrans.getGivensMatrixPreMultiple(0,&givensMatrix42));
	EXPECT_EQ(result, testGivensTrans.getGivensMatrixPreMultiple(1,&givensMatrix42));
	EXPECT_EQ(result, testGivensTrans.getGivensMatrixPreMultiple(2,&givensMatrix42));

	EXPECT_EQ(result, testGivensTrans.getGivensMatrixPreMultiple(0,&givensMatrix32));
	EXPECT_EQ(result, testGivensTrans.getGivensMatrixPreMultiple(1,&givensMatrix32));
	EXPECT_EQ(result, testGivensTrans.getGivensMatrixPreMultiple(2,&givensMatrix32));

	EXPECT_EQ(result, testGivensTrans.getGivensMatrixPreMultiple(0,&givensMatrix53));
	EXPECT_EQ(result, testGivensTrans.getGivensMatrixPreMultiple(1,&givensMatrix53));
	EXPECT_EQ(result, testGivensTrans.getGivensMatrixPreMultiple(2,&givensMatrix53));

	EXPECT_EQ(result, testGivensTrans.getGivensMatrixPreMultiple(3,&givensMatrix33));
	EXPECT_EQ(result, testGivensTrans.getGivensMatrixPreMultiple(4,&givensMatrix33));
}

/*
 * 测试GivensTransformation Givens矩阵Reload后求Pre转换
 */
TEST(GivensTransformationReloadTestPre, postive)
{
	StaticVector inputVector_old = StaticVector(3);
	inputVector_old.setElement(0,1.234);
	inputVector_old.setElement(1,2);
	inputVector_old.setElement(2,5);
	StaticMatrix vectorMatrix31_old = StaticMatrix(3,1);
	vectorMatrix31_old.setMatrixElement(0,0,1.234);
	vectorMatrix31_old.setMatrixElement(1,0,2);
	vectorMatrix31_old.setMatrixElement(2,0,5);

	StaticMatrix resultMatrix31_old = StaticMatrix(3,1);
	StaticMatrix givensMatrix33_old = StaticMatrix(3,3);

	GivensTransformation testGivensTrans = GivensTransformation(&inputVector_old);
	bool result = true;

	EXPECT_EQ(result, testGivensTrans.getGivensMatrixPreMultiple(0,&givensMatrix33_old));
	MatrixMultiplier multiVerify = MatrixMultiplier(&givensMatrix33_old,&vectorMatrix31_old,&resultMatrix31_old);
	multiVerify.multiplyCalc();
	resultMatrix31_old.regularZeroElement();
	EXPECT_EQ(0, resultMatrix31_old.getMatrixElement(0,0));
	EXPECT_NE(0, resultMatrix31_old.getMatrixElement(1,0));
	EXPECT_NE(0, resultMatrix31_old.getMatrixElement(2,0));

	StaticVector inputVector_new = StaticVector(3);
	inputVector_new.setElement(0,1);
	inputVector_new.setElement(1,3);
	inputVector_new.setElement(2,8);
	StaticMatrix vectorMatrix31_new = StaticMatrix(3,1);
	vectorMatrix31_new.setMatrixElement(0,0,1);
	vectorMatrix31_new.setMatrixElement(1,0,3);
	vectorMatrix31_new.setMatrixElement(2,0,8);

	StaticMatrix resultMatrix31_new = StaticMatrix(3,1);

	StaticMatrix givensMatrix33_new = StaticMatrix(3,3);

	//reload givens
	testGivensTrans.reload(&inputVector_new);
	EXPECT_EQ(result, testGivensTrans.getGivensMatrixPreMultiple(0,&givensMatrix33_new));
	multiVerify.reload(&givensMatrix33_new,&vectorMatrix31_new,&resultMatrix31_new);
	multiVerify.multiplyCalc();
	resultMatrix31_new.regularZeroElement();
	EXPECT_EQ(0, resultMatrix31_new.getMatrixElement(0,0));
	EXPECT_NE(0, resultMatrix31_new.getMatrixElement(1,0));
	EXPECT_NE(0, resultMatrix31_new.getMatrixElement(2,0));

	EXPECT_EQ(result, testGivensTrans.getGivensMatrixPreMultiple(1,&givensMatrix33_new));
	multiVerify.reload(&givensMatrix33_new,&vectorMatrix31_new,&resultMatrix31_new);
	multiVerify.multiplyCalc();
	resultMatrix31_new.regularZeroElement();
	EXPECT_NE(0, resultMatrix31_new.getMatrixElement(0,0));
	EXPECT_EQ(0, resultMatrix31_new.getMatrixElement(1,0));
	EXPECT_NE(0, resultMatrix31_new.getMatrixElement(2,0));

	EXPECT_EQ(result, testGivensTrans.getGivensMatrixPreMultiple(2,&givensMatrix33_new));
	multiVerify.reload(&givensMatrix33_new,&vectorMatrix31_new,&resultMatrix31_new);
	multiVerify.multiplyCalc();
	resultMatrix31_new.regularZeroElement();
	EXPECT_NE(0, resultMatrix31_new.getMatrixElement(0,0));
	EXPECT_NE(0, resultMatrix31_new.getMatrixElement(1,0));
	EXPECT_EQ(0, resultMatrix31_new.getMatrixElement(2,0));
}

/*
 * 测试GivensTransformation Givens矩阵Reload后求After转换
 */
TEST(GivensTransformationReloadTestAfter, postive)
{
	StaticVector inputVector_old = StaticVector(3);
	inputVector_old.setElement(0,1.234);
	inputVector_old.setElement(1,2);
	inputVector_old.setElement(2,5);
	StaticMatrix vectorMatrix13_old = StaticMatrix(1,3);
	vectorMatrix13_old.setMatrixElement(0,0,1.234);
	vectorMatrix13_old.setMatrixElement(0,1,2);
	vectorMatrix13_old.setMatrixElement(0,2,5);

	StaticMatrix resultMatrix13_old = StaticMatrix(1,3);

	StaticMatrix givensMatrix33_old = StaticMatrix(3,3);

	GivensTransformation testGivensTrans = GivensTransformation(&inputVector_old);

	bool result = true;

	EXPECT_EQ(result, testGivensTrans.getGivensMatrixAfterMultiple(0,&givensMatrix33_old));
	MatrixMultiplier multiVerify = MatrixMultiplier(&vectorMatrix13_old,&givensMatrix33_old,&resultMatrix13_old);
	multiVerify.multiplyCalc();
	resultMatrix13_old.regularZeroElement();
	EXPECT_EQ(0, resultMatrix13_old.getMatrixElement(0,0));
	EXPECT_NE(0, resultMatrix13_old.getMatrixElement(0,1));
	EXPECT_NE(0, resultMatrix13_old.getMatrixElement(0,2));

	StaticVector inputVector_new = StaticVector(3);
	inputVector_new.setElement(0,3);
	inputVector_new.setElement(1,5);
	inputVector_new.setElement(2,13);
	StaticMatrix vectorMatrix13_new = StaticMatrix(1,3);
	vectorMatrix13_new.setMatrixElement(0,0,3);
	vectorMatrix13_new.setMatrixElement(0,1,5);
	vectorMatrix13_new.setMatrixElement(0,2,13);

	StaticMatrix resultMatrix13_new = StaticMatrix(1,3);

	StaticMatrix givensMatrix33_new = StaticMatrix(3,3);

	//reload givens
	testGivensTrans.reload(&inputVector_new);
	EXPECT_EQ(result, testGivensTrans.getGivensMatrixAfterMultiple(0,&givensMatrix33_new));
	multiVerify.reload(&vectorMatrix13_new,&givensMatrix33_new,&resultMatrix13_new);
	multiVerify.multiplyCalc();
	resultMatrix13_new.regularZeroElement();
	EXPECT_EQ(0, resultMatrix13_new.getMatrixElement(0,0));
	EXPECT_NE(0, resultMatrix13_new.getMatrixElement(0,1));
	EXPECT_NE(0, resultMatrix13_new.getMatrixElement(0,2));

	EXPECT_EQ(result, testGivensTrans.getGivensMatrixAfterMultiple(1,&givensMatrix33_new));
	multiVerify.reload(&vectorMatrix13_new,&givensMatrix33_new,&resultMatrix13_new);
	multiVerify.multiplyCalc();
	resultMatrix13_new.regularZeroElement();
	EXPECT_NE(0, resultMatrix13_new.getMatrixElement(0,0));
	EXPECT_EQ(0, resultMatrix13_new.getMatrixElement(0,1));
	EXPECT_NE(0, resultMatrix13_new.getMatrixElement(0,2));

	EXPECT_EQ(result, testGivensTrans.getGivensMatrixAfterMultiple(2,&givensMatrix33_new));
	multiVerify.reload(&vectorMatrix13_new,&givensMatrix33_new,&resultMatrix13_new);
	multiVerify.multiplyCalc();
	resultMatrix13_new.regularZeroElement();
	EXPECT_NE(0, resultMatrix13_new.getMatrixElement(0,0));
	EXPECT_NE(0, resultMatrix13_new.getMatrixElement(0,1));
	EXPECT_EQ(0, resultMatrix13_new.getMatrixElement(0,2));
}

/*
 * 测试GivensTransformation 特殊向量011的转换
 */
TEST(GivensTransformation011VectorAfterTrans, postive)
{
	StaticVector inputVector = StaticVector(3);
	inputVector.setElement(0,0);
	inputVector.setElement(1,1);
	inputVector.setElement(2,1);

	StaticMatrix vectorMatrix13 = StaticMatrix(1,3);
	vectorMatrix13.setMatrixElement(0,0,0);
	vectorMatrix13.setMatrixElement(0,1,1);
	vectorMatrix13.setMatrixElement(0,2,1);

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

	EXPECT_EQ(1, givensMatrix33.getMatrixElement(0,0));
	EXPECT_EQ(0, givensMatrix33.getMatrixElement(0,1));
	EXPECT_EQ(0, givensMatrix33.getMatrixElement(0,2));

	EXPECT_EQ(0, givensMatrix33.getMatrixElement(1,0));
	EXPECT_EQ(1, givensMatrix33.getMatrixElement(1,1));
	EXPECT_EQ(0, givensMatrix33.getMatrixElement(1,2));

	EXPECT_EQ(0, givensMatrix33.getMatrixElement(2,0));
	EXPECT_EQ(0, givensMatrix33.getMatrixElement(2,1));
	EXPECT_EQ(1, givensMatrix33.getMatrixElement(2,2));

	EXPECT_EQ(result, testGivensTrans.getGivensMatrixAfterMultiple(1,&givensMatrix33));
	multiVerify.reload(&vectorMatrix13,&givensMatrix33,&resultMatrix13);
	multiVerify.multiplyCalc();
	resultMatrix13.regularZeroElement();
	EXPECT_EQ(0, resultMatrix13.getMatrixElement(0,0));
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

/*
 * 测试GivensTransformation 特殊向量011的转换
 */
TEST(GivensTransformation011VectorPreTrans, postive)
{
	StaticVector inputVector = StaticVector(3);
	inputVector.setElement(0,0);
	inputVector.setElement(1,1);
	inputVector.setElement(2,1);

	StaticMatrix vectorMatrix31 = StaticMatrix(3,1);
	vectorMatrix31.setMatrixElement(0,0,0);
	vectorMatrix31.setMatrixElement(1,0,1);
	vectorMatrix31.setMatrixElement(2,0,1);

	StaticMatrix resultMatrix31 = StaticMatrix(3,1);

	StaticMatrix givensMatrix33 = StaticMatrix(3,3);

	GivensTransformation testGivensTrans = GivensTransformation(&inputVector);

	bool result = true;

	EXPECT_EQ(result, testGivensTrans.getGivensMatrixPreMultiple(0,&givensMatrix33));
	MatrixMultiplier multiVerify = MatrixMultiplier(&givensMatrix33,&vectorMatrix31,&resultMatrix31);
	multiVerify.multiplyCalc();
	resultMatrix31.regularZeroElement();
	EXPECT_EQ(0, resultMatrix31.getMatrixElement(0,0));
	EXPECT_EQ(1, resultMatrix31.getMatrixElement(1,0));
	EXPECT_EQ(1, resultMatrix31.getMatrixElement(2,0));

	EXPECT_EQ(1, givensMatrix33.getMatrixElement(0,0));
	EXPECT_EQ(0, givensMatrix33.getMatrixElement(0,1));
	EXPECT_EQ(0, givensMatrix33.getMatrixElement(0,2));

	EXPECT_EQ(0, givensMatrix33.getMatrixElement(1,0));
	EXPECT_EQ(1, givensMatrix33.getMatrixElement(1,1));
	EXPECT_EQ(0, givensMatrix33.getMatrixElement(1,2));

	EXPECT_EQ(0, givensMatrix33.getMatrixElement(2,0));
	EXPECT_EQ(0, givensMatrix33.getMatrixElement(2,1));
	EXPECT_EQ(1, givensMatrix33.getMatrixElement(2,2));

	EXPECT_EQ(result, testGivensTrans.getGivensMatrixPreMultiple(1,&givensMatrix33));
	multiVerify.reload(&givensMatrix33,&vectorMatrix31,&resultMatrix31);
	multiVerify.multiplyCalc();
	resultMatrix31.regularZeroElement();
	EXPECT_NE(0, resultMatrix31.getMatrixElement(0,0));
	EXPECT_EQ(0, resultMatrix31.getMatrixElement(1,0));
	EXPECT_EQ(1, resultMatrix31.getMatrixElement(2,0));

	EXPECT_EQ(result, testGivensTrans.getGivensMatrixPreMultiple(2,&givensMatrix33));
	multiVerify.reload(&givensMatrix33,&vectorMatrix31,&resultMatrix31);
	multiVerify.multiplyCalc();
	resultMatrix31.regularZeroElement();
	EXPECT_EQ(0, resultMatrix31.getMatrixElement(0,0));
	EXPECT_NE(0, resultMatrix31.getMatrixElement(1,0));
	EXPECT_EQ(0, resultMatrix31.getMatrixElement(2,0));

}

/*
 * 测试GivensTransformation 特殊向量001的转换
 */
TEST(GivensTransformation001VectorAfterTrans, postive)
{
	StaticVector inputVector = StaticVector(3);
	inputVector.setElement(0,0);
	inputVector.setElement(1,0);
	inputVector.setElement(2,1);

	StaticMatrix vectorMatrix13 = StaticMatrix(1,3);
	vectorMatrix13.setMatrixElement(0,0,0);
	vectorMatrix13.setMatrixElement(0,1,0);
	vectorMatrix13.setMatrixElement(0,2,1);

	StaticMatrix resultMatrix13 = StaticMatrix(1,3);

	StaticMatrix givensMatrix33 = StaticMatrix(3,3);

	GivensTransformation testGivensTrans = GivensTransformation(&inputVector);

	bool result = true;

	EXPECT_EQ(result, testGivensTrans.getGivensMatrixAfterMultiple(0,&givensMatrix33));
	MatrixMultiplier multiVerify = MatrixMultiplier(&vectorMatrix13,&givensMatrix33,&resultMatrix13);
	multiVerify.multiplyCalc();
	resultMatrix13.regularZeroElement();
	EXPECT_EQ(0, resultMatrix13.getMatrixElement(0,0));
	EXPECT_EQ(0, resultMatrix13.getMatrixElement(0,1));
	EXPECT_NE(0, resultMatrix13.getMatrixElement(0,2));

	EXPECT_EQ(1, givensMatrix33.getMatrixElement(0,0));
	EXPECT_EQ(0, givensMatrix33.getMatrixElement(0,1));
	EXPECT_EQ(0, givensMatrix33.getMatrixElement(0,2));

	EXPECT_EQ(0, givensMatrix33.getMatrixElement(1,0));
	EXPECT_EQ(1, givensMatrix33.getMatrixElement(1,1));
	EXPECT_EQ(0, givensMatrix33.getMatrixElement(1,2));

	EXPECT_EQ(0, givensMatrix33.getMatrixElement(2,0));
	EXPECT_EQ(0, givensMatrix33.getMatrixElement(2,1));
	EXPECT_EQ(1, givensMatrix33.getMatrixElement(2,2));

	EXPECT_EQ(result, testGivensTrans.getGivensMatrixAfterMultiple(1,&givensMatrix33));
	multiVerify.reload(&vectorMatrix13,&givensMatrix33,&resultMatrix13);
	multiVerify.multiplyCalc();
	resultMatrix13.regularZeroElement();
	EXPECT_EQ(0, resultMatrix13.getMatrixElement(0,0));
	EXPECT_EQ(0, resultMatrix13.getMatrixElement(0,1));
	EXPECT_EQ(1, resultMatrix13.getMatrixElement(0,2));

	EXPECT_EQ(1, givensMatrix33.getMatrixElement(0,0));
	EXPECT_EQ(0, givensMatrix33.getMatrixElement(0,1));
	EXPECT_EQ(0, givensMatrix33.getMatrixElement(0,2));

	EXPECT_EQ(0, givensMatrix33.getMatrixElement(1,0));
	EXPECT_EQ(1, givensMatrix33.getMatrixElement(1,1));
	EXPECT_EQ(0, givensMatrix33.getMatrixElement(1,2));

	EXPECT_EQ(0, givensMatrix33.getMatrixElement(2,0));
	EXPECT_EQ(0, givensMatrix33.getMatrixElement(2,1));
	EXPECT_EQ(1, givensMatrix33.getMatrixElement(2,2));

	EXPECT_EQ(result, testGivensTrans.getGivensMatrixAfterMultiple(2,&givensMatrix33));
	multiVerify.reload(&vectorMatrix13,&givensMatrix33,&resultMatrix13);
	multiVerify.multiplyCalc();
	resultMatrix13.regularZeroElement();
	EXPECT_NE(0, resultMatrix13.getMatrixElement(0,0));
	EXPECT_EQ(0, resultMatrix13.getMatrixElement(0,1));
	EXPECT_EQ(0, resultMatrix13.getMatrixElement(0,2));
}

/*
 * 测试GivensTransformation 特殊向量001的转换
 */
TEST(GivensTransformation001VectorPreTrans, postive)
{
	StaticVector inputVector = StaticVector(3);
	inputVector.setElement(0,0);
	inputVector.setElement(1,0);
	inputVector.setElement(2,1);

	StaticMatrix vectorMatrix31 = StaticMatrix(3,1);
	vectorMatrix31.setMatrixElement(0,0,0);
	vectorMatrix31.setMatrixElement(1,0,0);
	vectorMatrix31.setMatrixElement(2,0,1);

	StaticMatrix resultMatrix31 = StaticMatrix(3,1);

	StaticMatrix givensMatrix33 = StaticMatrix(3,3);

	GivensTransformation testGivensTrans = GivensTransformation(&inputVector);

	bool result = true;

	EXPECT_EQ(result, testGivensTrans.getGivensMatrixPreMultiple(0,&givensMatrix33));
	MatrixMultiplier multiVerify = MatrixMultiplier(&givensMatrix33,&vectorMatrix31,&resultMatrix31);
	multiVerify.multiplyCalc();
	resultMatrix31.regularZeroElement();
	EXPECT_EQ(0, resultMatrix31.getMatrixElement(0,0));
	EXPECT_EQ(0, resultMatrix31.getMatrixElement(1,0));
	EXPECT_EQ(1, resultMatrix31.getMatrixElement(2,0));

	EXPECT_EQ(1, givensMatrix33.getMatrixElement(0,0));
	EXPECT_EQ(0, givensMatrix33.getMatrixElement(0,1));
	EXPECT_EQ(0, givensMatrix33.getMatrixElement(0,2));

	EXPECT_EQ(0, givensMatrix33.getMatrixElement(1,0));
	EXPECT_EQ(1, givensMatrix33.getMatrixElement(1,1));
	EXPECT_EQ(0, givensMatrix33.getMatrixElement(1,2));

	EXPECT_EQ(0, givensMatrix33.getMatrixElement(2,0));
	EXPECT_EQ(0, givensMatrix33.getMatrixElement(2,1));
	EXPECT_EQ(1, givensMatrix33.getMatrixElement(2,2));

	EXPECT_EQ(result, testGivensTrans.getGivensMatrixPreMultiple(1,&givensMatrix33));
	multiVerify.reload(&givensMatrix33,&vectorMatrix31,&resultMatrix31);
	multiVerify.multiplyCalc();
	resultMatrix31.regularZeroElement();
	EXPECT_EQ(0, resultMatrix31.getMatrixElement(0,0));
	EXPECT_EQ(0, resultMatrix31.getMatrixElement(1,0));
	EXPECT_EQ(1, resultMatrix31.getMatrixElement(2,0));

	EXPECT_EQ(1, givensMatrix33.getMatrixElement(0,0));
	EXPECT_EQ(0, givensMatrix33.getMatrixElement(0,1));
	EXPECT_EQ(0, givensMatrix33.getMatrixElement(0,2));

	EXPECT_EQ(0, givensMatrix33.getMatrixElement(1,0));
	EXPECT_EQ(1, givensMatrix33.getMatrixElement(1,1));
	EXPECT_EQ(0, givensMatrix33.getMatrixElement(1,2));

	EXPECT_EQ(0, givensMatrix33.getMatrixElement(2,0));
	EXPECT_EQ(0, givensMatrix33.getMatrixElement(2,1));
	EXPECT_EQ(1, givensMatrix33.getMatrixElement(2,2));

	EXPECT_EQ(result, testGivensTrans.getGivensMatrixPreMultiple(2,&givensMatrix33));
	multiVerify.reload(&givensMatrix33,&vectorMatrix31,&resultMatrix31);
	multiVerify.multiplyCalc();
	resultMatrix31.regularZeroElement();
	EXPECT_EQ(0, resultMatrix31.getMatrixElement(0,0));
	EXPECT_NE(0, resultMatrix31.getMatrixElement(1,0));
	EXPECT_EQ(0, resultMatrix31.getMatrixElement(2,0));

}


/*
 * 测试GivensTransformation 特殊向量000的转换
 */
TEST(GivensTransformation000VectorAfterTrans, postive)
{
	StaticVector inputVector = StaticVector(3);
	inputVector.setElement(0,0);
	inputVector.setElement(1,0);
	inputVector.setElement(2,0);

	StaticMatrix vectorMatrix13 = StaticMatrix(1,3);
	vectorMatrix13.setMatrixElement(0,0,0);
	vectorMatrix13.setMatrixElement(0,1,0);
	vectorMatrix13.setMatrixElement(0,2,0);

	StaticMatrix resultMatrix13 = StaticMatrix(1,3);

	StaticMatrix givensMatrix33 = StaticMatrix(3,3);

	GivensTransformation testGivensTrans = GivensTransformation(&inputVector);

	bool result = true;

	EXPECT_EQ(result, testGivensTrans.getGivensMatrixAfterMultiple(0,&givensMatrix33));
	MatrixMultiplier multiVerify = MatrixMultiplier(&vectorMatrix13,&givensMatrix33,&resultMatrix13);
	multiVerify.multiplyCalc();
	resultMatrix13.regularZeroElement();
	EXPECT_EQ(0, resultMatrix13.getMatrixElement(0,0));
	EXPECT_EQ(0, resultMatrix13.getMatrixElement(0,1));
	EXPECT_EQ(0, resultMatrix13.getMatrixElement(0,2));

	EXPECT_EQ(1, givensMatrix33.getMatrixElement(0,0));
	EXPECT_EQ(0, givensMatrix33.getMatrixElement(0,1));
	EXPECT_EQ(0, givensMatrix33.getMatrixElement(0,2));

	EXPECT_EQ(0, givensMatrix33.getMatrixElement(1,0));
	EXPECT_EQ(1, givensMatrix33.getMatrixElement(1,1));
	EXPECT_EQ(0, givensMatrix33.getMatrixElement(1,2));

	EXPECT_EQ(0, givensMatrix33.getMatrixElement(2,0));
	EXPECT_EQ(0, givensMatrix33.getMatrixElement(2,1));
	EXPECT_EQ(1, givensMatrix33.getMatrixElement(2,2));

	EXPECT_EQ(result, testGivensTrans.getGivensMatrixAfterMultiple(1,&givensMatrix33));
	multiVerify.reload(&vectorMatrix13,&givensMatrix33,&resultMatrix13);
	multiVerify.multiplyCalc();
	resultMatrix13.regularZeroElement();
	EXPECT_EQ(0, resultMatrix13.getMatrixElement(0,0));
	EXPECT_EQ(0, resultMatrix13.getMatrixElement(0,1));
	EXPECT_EQ(0, resultMatrix13.getMatrixElement(0,2));

	EXPECT_EQ(1, givensMatrix33.getMatrixElement(0,0));
	EXPECT_EQ(0, givensMatrix33.getMatrixElement(0,1));
	EXPECT_EQ(0, givensMatrix33.getMatrixElement(0,2));

	EXPECT_EQ(0, givensMatrix33.getMatrixElement(1,0));
	EXPECT_EQ(1, givensMatrix33.getMatrixElement(1,1));
	EXPECT_EQ(0, givensMatrix33.getMatrixElement(1,2));

	EXPECT_EQ(0, givensMatrix33.getMatrixElement(2,0));
	EXPECT_EQ(0, givensMatrix33.getMatrixElement(2,1));
	EXPECT_EQ(1, givensMatrix33.getMatrixElement(2,2));

	EXPECT_EQ(result, testGivensTrans.getGivensMatrixAfterMultiple(2,&givensMatrix33));
	multiVerify.reload(&vectorMatrix13,&givensMatrix33,&resultMatrix13);
	multiVerify.multiplyCalc();
	resultMatrix13.regularZeroElement();
	EXPECT_EQ(0, resultMatrix13.getMatrixElement(0,0));
	EXPECT_EQ(0, resultMatrix13.getMatrixElement(0,1));
	EXPECT_EQ(0, resultMatrix13.getMatrixElement(0,2));

	EXPECT_EQ(1, givensMatrix33.getMatrixElement(0,0));
	EXPECT_EQ(0, givensMatrix33.getMatrixElement(0,1));
	EXPECT_EQ(0, givensMatrix33.getMatrixElement(0,2));

	EXPECT_EQ(0, givensMatrix33.getMatrixElement(1,0));
	EXPECT_EQ(1, givensMatrix33.getMatrixElement(1,1));
	EXPECT_EQ(0, givensMatrix33.getMatrixElement(1,2));

	EXPECT_EQ(0, givensMatrix33.getMatrixElement(2,0));
	EXPECT_EQ(0, givensMatrix33.getMatrixElement(2,1));
	EXPECT_EQ(1, givensMatrix33.getMatrixElement(2,2));
}

/*
 * 测试GivensTransformation 特殊向量000的转换
 */
TEST(GivensTransformation000VectorPreTrans, postive)
{
	StaticVector inputVector = StaticVector(3);
	inputVector.setElement(0,0);
	inputVector.setElement(1,0);
	inputVector.setElement(2,0);

	StaticMatrix vectorMatrix31 = StaticMatrix(3,1);
	vectorMatrix31.setMatrixElement(0,0,0);
	vectorMatrix31.setMatrixElement(1,0,0);
	vectorMatrix31.setMatrixElement(2,0,0);

	StaticMatrix resultMatrix31 = StaticMatrix(3,1);

	StaticMatrix givensMatrix33 = StaticMatrix(3,3);

	GivensTransformation testGivensTrans = GivensTransformation(&inputVector);

	bool result = true;

	EXPECT_EQ(result, testGivensTrans.getGivensMatrixPreMultiple(0,&givensMatrix33));
	MatrixMultiplier multiVerify = MatrixMultiplier(&givensMatrix33,&vectorMatrix31,&resultMatrix31);
	multiVerify.multiplyCalc();
	resultMatrix31.regularZeroElement();
	EXPECT_EQ(0, resultMatrix31.getMatrixElement(0,0));
	EXPECT_EQ(0, resultMatrix31.getMatrixElement(1,0));
	EXPECT_EQ(0, resultMatrix31.getMatrixElement(2,0));

	EXPECT_EQ(1, givensMatrix33.getMatrixElement(0,0));
	EXPECT_EQ(0, givensMatrix33.getMatrixElement(0,1));
	EXPECT_EQ(0, givensMatrix33.getMatrixElement(0,2));

	EXPECT_EQ(0, givensMatrix33.getMatrixElement(1,0));
	EXPECT_EQ(1, givensMatrix33.getMatrixElement(1,1));
	EXPECT_EQ(0, givensMatrix33.getMatrixElement(1,2));

	EXPECT_EQ(0, givensMatrix33.getMatrixElement(2,0));
	EXPECT_EQ(0, givensMatrix33.getMatrixElement(2,1));
	EXPECT_EQ(1, givensMatrix33.getMatrixElement(2,2));

	EXPECT_EQ(result, testGivensTrans.getGivensMatrixPreMultiple(1,&givensMatrix33));
	multiVerify.reload(&givensMatrix33,&vectorMatrix31,&resultMatrix31);
	multiVerify.multiplyCalc();
	resultMatrix31.regularZeroElement();
	EXPECT_EQ(0, resultMatrix31.getMatrixElement(0,0));
	EXPECT_EQ(0, resultMatrix31.getMatrixElement(1,0));
	EXPECT_EQ(0, resultMatrix31.getMatrixElement(2,0));

	EXPECT_EQ(1, givensMatrix33.getMatrixElement(0,0));
	EXPECT_EQ(0, givensMatrix33.getMatrixElement(0,1));
	EXPECT_EQ(0, givensMatrix33.getMatrixElement(0,2));

	EXPECT_EQ(0, givensMatrix33.getMatrixElement(1,0));
	EXPECT_EQ(1, givensMatrix33.getMatrixElement(1,1));
	EXPECT_EQ(0, givensMatrix33.getMatrixElement(1,2));

	EXPECT_EQ(0, givensMatrix33.getMatrixElement(2,0));
	EXPECT_EQ(0, givensMatrix33.getMatrixElement(2,1));
	EXPECT_EQ(1, givensMatrix33.getMatrixElement(2,2));

	EXPECT_EQ(result, testGivensTrans.getGivensMatrixPreMultiple(2,&givensMatrix33));
	multiVerify.reload(&givensMatrix33,&vectorMatrix31,&resultMatrix31);
	multiVerify.multiplyCalc();
	resultMatrix31.regularZeroElement();
	EXPECT_EQ(0, resultMatrix31.getMatrixElement(0,0));
	EXPECT_EQ(0, resultMatrix31.getMatrixElement(1,0));
	EXPECT_EQ(0, resultMatrix31.getMatrixElement(2,0));

	EXPECT_EQ(1, givensMatrix33.getMatrixElement(0,0));
	EXPECT_EQ(0, givensMatrix33.getMatrixElement(0,1));
	EXPECT_EQ(0, givensMatrix33.getMatrixElement(0,2));

	EXPECT_EQ(0, givensMatrix33.getMatrixElement(1,0));
	EXPECT_EQ(1, givensMatrix33.getMatrixElement(1,1));
	EXPECT_EQ(0, givensMatrix33.getMatrixElement(1,2));

	EXPECT_EQ(0, givensMatrix33.getMatrixElement(2,0));
	EXPECT_EQ(0, givensMatrix33.getMatrixElement(2,1));
	EXPECT_EQ(1, givensMatrix33.getMatrixElement(2,2));
}
