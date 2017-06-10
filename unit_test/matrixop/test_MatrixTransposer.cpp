/*
 * test_MatrixTransposer.cpp
 *
 *  Created on: 2017ƒÍ6‘¬8»’
 *      Author: looke
 */

#include "..\gtest_src\gtest\gtest.h"
#include "StaticMatrix.h"
#include "MatrixTransposer.h"
#include "math.h"

/*
 * ≤‚ ‘MatrixTransposer æÿ’Û◊™÷√∑Ω’Û
 */
TEST(MatrixTransposerTransposeSquareTest, postive)
{
	StaticMatrix testMatrix33 = StaticMatrix(3,3);
	testMatrix33.resetMatrixToI();

	MatrixTransposer testTransposer = MatrixTransposer();
	testTransposer.transposeSquareMatrix(&testMatrix33);

	EXPECT_EQ(1, testMatrix33.getMatrixElement(0,0));
	EXPECT_EQ(0, testMatrix33.getMatrixElement(0,1));
	EXPECT_EQ(0, testMatrix33.getMatrixElement(0,2));

	EXPECT_EQ(0, testMatrix33.getMatrixElement(1,0));
	EXPECT_EQ(1, testMatrix33.getMatrixElement(1,1));
	EXPECT_EQ(0, testMatrix33.getMatrixElement(1,2));

	EXPECT_EQ(0, testMatrix33.getMatrixElement(2,0));
	EXPECT_EQ(0, testMatrix33.getMatrixElement(2,1));
	EXPECT_EQ(1, testMatrix33.getMatrixElement(2,2));

	StaticMatrix testMatrix44 = StaticMatrix(4,4);
	testMatrix44.setMatrixElement(0,0,1);
	testMatrix44.setMatrixElement(0,1,2);
	testMatrix44.setMatrixElement(0,2,3);
	testMatrix44.setMatrixElement(0,3,4);

	testMatrix44.setMatrixElement(1,0,5);
	testMatrix44.setMatrixElement(1,1,6);
	testMatrix44.setMatrixElement(1,2,7);
	testMatrix44.setMatrixElement(1,3,8);

	testMatrix44.setMatrixElement(2,0,9);
	testMatrix44.setMatrixElement(2,1,10);
	testMatrix44.setMatrixElement(2,2,11);
	testMatrix44.setMatrixElement(2,3,12);

	testMatrix44.setMatrixElement(3,0,13);
	testMatrix44.setMatrixElement(3,1,14);
	testMatrix44.setMatrixElement(3,2,15);
	testMatrix44.setMatrixElement(3,3,16);

	testTransposer.transposeSquareMatrix(&testMatrix44);
	EXPECT_EQ(1, testMatrix44.getMatrixElement(0,0));
	EXPECT_EQ(2, testMatrix44.getMatrixElement(1,0));
	EXPECT_EQ(3, testMatrix44.getMatrixElement(2,0));
	EXPECT_EQ(4, testMatrix44.getMatrixElement(3,0));

	EXPECT_EQ(5, testMatrix44.getMatrixElement(0,1));
	EXPECT_EQ(6, testMatrix44.getMatrixElement(1,1));
	EXPECT_EQ(7, testMatrix44.getMatrixElement(2,1));
	EXPECT_EQ(8, testMatrix44.getMatrixElement(3,1));

	EXPECT_EQ(9, testMatrix44.getMatrixElement(0,2));
	EXPECT_EQ(10, testMatrix44.getMatrixElement(1,2));
	EXPECT_EQ(11, testMatrix44.getMatrixElement(2,2));
	EXPECT_EQ(12, testMatrix44.getMatrixElement(3,2));

	EXPECT_EQ(13, testMatrix44.getMatrixElement(0,3));
	EXPECT_EQ(14, testMatrix44.getMatrixElement(1,3));
	EXPECT_EQ(15, testMatrix44.getMatrixElement(2,3));
	EXPECT_EQ(16, testMatrix44.getMatrixElement(3,3));
}


/*
 * ≤‚ ‘MatrixTransposer æÿ’Û◊™÷√ ∑«∑Ω’Û“Ï≥£
 */
TEST(MatrixTransposerNotSquareTest, postive)
{
	StaticMatrix testMatrix45 = StaticMatrix(4,5);
	MatrixTransposer testTransposer = MatrixTransposer();

	EXPECT_EQ(0,testTransposer.transposeSquareMatrix(&testMatrix45));

	EXPECT_EQ(1, testMatrix45.getMatrixElement(0,0));
	EXPECT_EQ(2, testMatrix45.getMatrixElement(0,1));
	EXPECT_EQ(3, testMatrix45.getMatrixElement(0,2));
	EXPECT_EQ(4, testMatrix45.getMatrixElement(0,3));
	EXPECT_EQ(5, testMatrix45.getMatrixElement(0,4));


	EXPECT_EQ(6, testMatrix45.getMatrixElement(1,0));
	EXPECT_EQ(7, testMatrix45.getMatrixElement(1,1));
	EXPECT_EQ(8, testMatrix45.getMatrixElement(1,2));
	EXPECT_EQ(9, testMatrix45.getMatrixElement(1,3));
	EXPECT_EQ(10,testMatrix45.getMatrixElement(1,4));

	EXPECT_EQ(11, testMatrix45.getMatrixElement(2,0));
	EXPECT_EQ(12, testMatrix45.getMatrixElement(2,1));
	EXPECT_EQ(13, testMatrix45.getMatrixElement(2,2));
	EXPECT_EQ(14, testMatrix45.getMatrixElement(2,3));
	EXPECT_EQ(15, testMatrix45.getMatrixElement(2,4));

	EXPECT_EQ(16, testMatrix45.getMatrixElement(3,0));
	EXPECT_EQ(17, testMatrix45.getMatrixElement(3,1));
	EXPECT_EQ(18, testMatrix45.getMatrixElement(3,2));
	EXPECT_EQ(19, testMatrix45.getMatrixElement(3,3));
	EXPECT_EQ(20, testMatrix45.getMatrixElement(3,4));
}

/*
 * ≤‚ ‘MatrixTransposer æÿ’Û◊™÷√»Œ“‚æÿ’Û
 */
TEST(MatrixTransposerTransposeTest, postive)
{
	StaticMatrix testMatrix45 = StaticMatrix(4,5);
	StaticMatrix testMatrix54 = StaticMatrix(5,4);

	MatrixTransposer testTransposer = MatrixTransposer();
	testTransposer.transposeMatrix(&testMatrix45, &testMatrix54);

	EXPECT_EQ(1, testMatrix54.getMatrixElement(0,0));
	EXPECT_EQ(2, testMatrix54.getMatrixElement(1,0));
	EXPECT_EQ(3, testMatrix54.getMatrixElement(2,0));
	EXPECT_EQ(4, testMatrix54.getMatrixElement(3,0));
	EXPECT_EQ(5, testMatrix54.getMatrixElement(4,0));

	EXPECT_EQ(6, testMatrix54.getMatrixElement(0,1));
	EXPECT_EQ(7, testMatrix54.getMatrixElement(1,1));
	EXPECT_EQ(8, testMatrix54.getMatrixElement(2,1));
	EXPECT_EQ(9, testMatrix54.getMatrixElement(3,1));
	EXPECT_EQ(10,testMatrix54.getMatrixElement(4,1));

	EXPECT_EQ(11, testMatrix54.getMatrixElement(0,2));
	EXPECT_EQ(12, testMatrix54.getMatrixElement(1,2));
	EXPECT_EQ(13, testMatrix54.getMatrixElement(2,2));
	EXPECT_EQ(14, testMatrix54.getMatrixElement(3,2));
	EXPECT_EQ(15, testMatrix54.getMatrixElement(4,2));

	EXPECT_EQ(16, testMatrix54.getMatrixElement(0,3));
	EXPECT_EQ(17, testMatrix54.getMatrixElement(1,3));
	EXPECT_EQ(18, testMatrix54.getMatrixElement(2,3));
	EXPECT_EQ(19, testMatrix54.getMatrixElement(3,3));
	EXPECT_EQ(20, testMatrix54.getMatrixElement(4,3));
}

/*
 * ≤‚ ‘MatrixTransposer æÿ’ÛŒ¨∂»“Ï≥£
 */
TEST(MatrixTransposerMatrixNotMatchExceptionTest, postive)
{
	StaticMatrix testMatrix45 = StaticMatrix(4,5);
	StaticMatrix testMatrix44 = StaticMatrix(4,4);

	StaticMatrix testMatrix65 = StaticMatrix(6,5);
	StaticMatrix testMatrix54 = StaticMatrix(5,4);


	MatrixTransposer testTransposer = MatrixTransposer();

	EXPECT_EQ(0,testTransposer.transposeMatrix(&testMatrix45, &testMatrix44));
	EXPECT_EQ(0,testTransposer.transposeMatrix(&testMatrix65, &testMatrix54));
}
