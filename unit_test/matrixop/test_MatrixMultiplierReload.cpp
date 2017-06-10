/*
 * test_StaticMatrixMultiplierReload.cpp
 *
 *  Created on: 2017��6��8��
 *      Author: looke
 */


#include "..\gtest_src\gtest\gtest.h"
#include "StaticMatrix.h"
#include "MatrixMultiplier.h"
#include "math.h"


class MatrixMultiplierReloadTest:public testing::Test
{
public:
	MatrixMultiplierReloadTest():leftMatrix_old(2,3),rightMatrix_old(3,5),resultMatrix(2,5),sMatrixMultiplier(&leftMatrix_old, &rightMatrix_old, &resultMatrix)
	{}
protected:
    void SetUp()
    {
    	sMatrixMultiplier.multiplyCalc();
    }
    void TearDown()
    {

    }

    StaticMatrix leftMatrix_old;
    StaticMatrix rightMatrix_old;
    StaticMatrix resultMatrix;

    MatrixMultiplier sMatrixMultiplier;
};

/*
 * �������¼��غ�ľ���˷�
 */
TEST_F(MatrixMultiplierReloadTest,MultiplierMulticalcTest)
{
	StaticMatrix leftMatrix44 = StaticMatrix(4,4);
	leftMatrix44.setMatrixElement(0,0,1);
	leftMatrix44.setMatrixElement(0,1,2);
	leftMatrix44.setMatrixElement(0,2,3);
	leftMatrix44.setMatrixElement(0,3,4);

	leftMatrix44.setMatrixElement(1,0,5);
	leftMatrix44.setMatrixElement(1,1,6);
	leftMatrix44.setMatrixElement(1,2,7);
	leftMatrix44.setMatrixElement(1,3,8);

	leftMatrix44.setMatrixElement(2,0,9);
	leftMatrix44.setMatrixElement(2,1,10);
	leftMatrix44.setMatrixElement(2,2,11);
	leftMatrix44.setMatrixElement(2,3,12);

	leftMatrix44.setMatrixElement(3,0,13);
	leftMatrix44.setMatrixElement(3,1,14);
	leftMatrix44.setMatrixElement(3,2,15);
	leftMatrix44.setMatrixElement(3,3,16);

	StaticMatrix rightMatrix43 = StaticMatrix(4,3);
	rightMatrix43.setMatrixElement(0,0,1);
	rightMatrix43.setMatrixElement(0,1,2);
	rightMatrix43.setMatrixElement(0,2,3);

	rightMatrix43.setMatrixElement(1,0,1);
	rightMatrix43.setMatrixElement(1,1,2);
	rightMatrix43.setMatrixElement(1,2,3);

	rightMatrix43.setMatrixElement(2,0,1);
	rightMatrix43.setMatrixElement(2,1,2);
	rightMatrix43.setMatrixElement(2,2,3);

	rightMatrix43.setMatrixElement(3,0,1);
	rightMatrix43.setMatrixElement(3,1,2);
	rightMatrix43.setMatrixElement(3,2,3);

	StaticMatrix resultMatrix43 = StaticMatrix(4,3);

	sMatrixMultiplier.reload(&leftMatrix44, &rightMatrix43, &resultMatrix43);

	bool result = true;

	EXPECT_EQ(result, sMatrixMultiplier.multiplyCalc());
	EXPECT_EQ(10, resultMatrix43.getMatrixElement(0,0));
	EXPECT_EQ(20, resultMatrix43.getMatrixElement(0,1));
	EXPECT_EQ(30, resultMatrix43.getMatrixElement(0,2));

	EXPECT_EQ(26, resultMatrix43.getMatrixElement(1,0));
	EXPECT_EQ(52, resultMatrix43.getMatrixElement(1,1));
	EXPECT_EQ(78, resultMatrix43.getMatrixElement(1,2));

	EXPECT_EQ(42, resultMatrix43.getMatrixElement(2,0));
	EXPECT_EQ(84, resultMatrix43.getMatrixElement(2,1));
	EXPECT_EQ(126, resultMatrix43.getMatrixElement(2,2));

	EXPECT_EQ(58, resultMatrix43.getMatrixElement(3,0));
	EXPECT_EQ(116, resultMatrix43.getMatrixElement(3,1));
	EXPECT_EQ(174, resultMatrix43.getMatrixElement(3,2));
}


/*
 * �������¼��غ�ľ���˷��쳣����
 */
TEST_F(MatrixMultiplierReloadTest, MatrixMultiplierMultiCalcExceptionTest)
{
	StaticMatrix leftMatrix32 = StaticMatrix(3,2);
	StaticMatrix rightMatrix33 = StaticMatrix(3,3);
	StaticMatrix resultMatrix33 = StaticMatrix(3,3);
	sMatrixMultiplier.reload(&leftMatrix32, &rightMatrix33, &resultMatrix33);
	bool result = false;
	EXPECT_EQ(result, sMatrixMultiplier.multiplyCalc());
}
