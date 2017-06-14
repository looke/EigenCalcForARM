/*
 * test_QRDecomposition.cpp
 *
 *  Created on: 2017年6月12日
 *      Author: looke
 */
#include "..\gtest_src\gtest\gtest.h"
#include "StaticVector.h"
#include "StaticMatrix.h"
#include "MatrixMultiplier.h"
#include "QRDecomposition.h"
#include "math.h"

/*
 * 测试QRDecomposition QR分解
 */
TEST(QRDecompositionQRDecompTest, postive)
{
	bool result = true;
	double lowEdge;
	StaticMatrix test44_upTriangle = StaticMatrix(4,4);
	test44_upTriangle.setMatrixElement(0,0,1);
	test44_upTriangle.setMatrixElement(0,1,2);
	test44_upTriangle.setMatrixElement(0,2,3);
	test44_upTriangle.setMatrixElement(0,3,4);

	test44_upTriangle.setMatrixElement(1,0,0);
	test44_upTriangle.setMatrixElement(1,1,5);
	test44_upTriangle.setMatrixElement(1,2,6);
	test44_upTriangle.setMatrixElement(1,3,7);

	test44_upTriangle.setMatrixElement(2,0,0);
	test44_upTriangle.setMatrixElement(2,1,0);
	test44_upTriangle.setMatrixElement(2,2,8);
	test44_upTriangle.setMatrixElement(2,3,9);

	test44_upTriangle.setMatrixElement(3,0,0);
	test44_upTriangle.setMatrixElement(3,1,0);
	test44_upTriangle.setMatrixElement(3,2,0);
	test44_upTriangle.setMatrixElement(3,3,10);

	StaticMatrix test44_QMatrix = StaticMatrix(4,4);
	StaticMatrix test44_HouseholderTransMatrix = StaticMatrix(4,4);
	StaticMatrix test44_TempMatrix = StaticMatrix(4,4);

	QRDecomposition qRDecompTest = QRDecomposition(&test44_upTriangle, &test44_QMatrix, &test44_HouseholderTransMatrix, &test44_TempMatrix);

	qRDecompTest.calcQRMatrix();

	EXPECT_EQ(-1, test44_upTriangle.getMatrixElement(0,0));
	EXPECT_EQ(-2, test44_upTriangle.getMatrixElement(0,1));
	EXPECT_EQ(-3, test44_upTriangle.getMatrixElement(0,2));
	EXPECT_EQ(-4, test44_upTriangle.getMatrixElement(0,3));

	EXPECT_EQ(0, test44_upTriangle.getMatrixElement(1,0));
	EXPECT_EQ(-5, test44_upTriangle.getMatrixElement(1,1));
	EXPECT_EQ(-6, test44_upTriangle.getMatrixElement(1,2));
	EXPECT_EQ(-7, test44_upTriangle.getMatrixElement(1,3));

	EXPECT_EQ(0, test44_upTriangle.getMatrixElement(2,0));
	EXPECT_EQ(0, test44_upTriangle.getMatrixElement(2,1));
	EXPECT_EQ(-8, test44_upTriangle.getMatrixElement(2,2));
	EXPECT_EQ(-9, test44_upTriangle.getMatrixElement(2,3));

	EXPECT_EQ(0, test44_upTriangle.getMatrixElement(3,0));
	EXPECT_EQ(0, test44_upTriangle.getMatrixElement(3,1));
	EXPECT_EQ(0, test44_upTriangle.getMatrixElement(3,2));
	EXPECT_EQ(10, test44_upTriangle.getMatrixElement(3,3));

	MatrixMultiplier multi = MatrixMultiplier(&test44_QMatrix,&test44_upTriangle, &test44_TempMatrix);
	multi.multiplyCalc();
	lowEdge = test44_TempMatrix.getLowEdge();
	EXPECT_EQ(result, fabs(test44_TempMatrix.getMatrixElement(0,0)-(1))<lowEdge);
	EXPECT_EQ(result, fabs(test44_TempMatrix.getMatrixElement(0,1)-(2))<lowEdge);
	EXPECT_EQ(result, fabs(test44_TempMatrix.getMatrixElement(0,2)-(3))<lowEdge);
	EXPECT_EQ(result, fabs(test44_TempMatrix.getMatrixElement(0,3)-(4))<lowEdge);

	EXPECT_EQ(result, fabs(test44_TempMatrix.getMatrixElement(1,0)-(0))<lowEdge);
	EXPECT_EQ(result, fabs(test44_TempMatrix.getMatrixElement(1,1)-(5))<lowEdge);
	EXPECT_EQ(result, fabs(test44_TempMatrix.getMatrixElement(1,2)-(6))<lowEdge);
	EXPECT_EQ(result, fabs(test44_TempMatrix.getMatrixElement(1,3)-(7))<lowEdge);

	EXPECT_EQ(result, fabs(test44_TempMatrix.getMatrixElement(2,0)-(0))<lowEdge);
	EXPECT_EQ(result, fabs(test44_TempMatrix.getMatrixElement(2,1)-(0))<lowEdge);
	EXPECT_EQ(result, fabs(test44_TempMatrix.getMatrixElement(2,2)-(8))<lowEdge);
	EXPECT_EQ(result, fabs(test44_TempMatrix.getMatrixElement(2,3)-(9))<lowEdge);

	EXPECT_EQ(result, fabs(test44_TempMatrix.getMatrixElement(3,0)-(0))<lowEdge);
	EXPECT_EQ(result, fabs(test44_TempMatrix.getMatrixElement(3,1)-(0))<lowEdge);
	EXPECT_EQ(result, fabs(test44_TempMatrix.getMatrixElement(3,2)-(0))<lowEdge);
	EXPECT_EQ(result, fabs(test44_TempMatrix.getMatrixElement(3,3)-(10))<lowEdge);

	StaticMatrix test33 = StaticMatrix(3,3);
	test33.setMatrixElement(0,0,1);
	test33.setMatrixElement(0,1,1);
	test33.setMatrixElement(0,2,1);

	test33.setMatrixElement(1,0,2);
	test33.setMatrixElement(1,1,-1);
	test33.setMatrixElement(1,2,-1);

	test33.setMatrixElement(2,0,2);
	test33.setMatrixElement(2,1,-4);
	test33.setMatrixElement(2,2,5);

	StaticMatrix test3_QMatrix = StaticMatrix(3,3);
	StaticMatrix test3_HouseholderTransMatrix = StaticMatrix(3,3);
	StaticMatrix test3_TempMatrix = StaticMatrix(3,3);

	qRDecompTest.reload(&test33,&test3_QMatrix,&test3_HouseholderTransMatrix,&test3_TempMatrix);
	qRDecompTest.calcQRMatrix();
	test33.regularZeroElement();

	lowEdge = test33.getLowEdge();

	EXPECT_EQ(result, fabs(test33.getMatrixElement(0,0)-(-3))<lowEdge);
	EXPECT_EQ(result, fabs(test33.getMatrixElement(0,1)-(3))<lowEdge);
	EXPECT_EQ(result, fabs(test33.getMatrixElement(0,2)-(-3))<lowEdge);

	EXPECT_EQ(result, fabs(test33.getMatrixElement(1,0)-(0))<lowEdge);
	EXPECT_EQ(result, fabs(test33.getMatrixElement(1,1)-(3))<lowEdge);
	EXPECT_EQ(result, fabs(test33.getMatrixElement(1,2)-(-3))<lowEdge);

	EXPECT_EQ(result, fabs(test33.getMatrixElement(2,0)-(0))<lowEdge);
	EXPECT_EQ(result, fabs(test33.getMatrixElement(2,1)-(0))<lowEdge);
	EXPECT_EQ(result, fabs(test33.getMatrixElement(2,2)-(3))<lowEdge);


	multi.reload(&test3_QMatrix,&test33, &test3_TempMatrix);
	multi.multiplyCalc();
	lowEdge = test3_TempMatrix.getLowEdge();

	EXPECT_EQ(result, fabs(test3_TempMatrix.getMatrixElement(0,0)-(1))<lowEdge);
	EXPECT_EQ(result, fabs(test3_TempMatrix.getMatrixElement(0,1)-(1))<lowEdge);
	EXPECT_EQ(result, fabs(test3_TempMatrix.getMatrixElement(0,2)-(1))<lowEdge);

	EXPECT_EQ(result, fabs(test3_TempMatrix.getMatrixElement(1,0)-(2))<lowEdge);
	EXPECT_EQ(result, fabs(test3_TempMatrix.getMatrixElement(1,1)-(-1))<lowEdge);
	EXPECT_EQ(result, fabs(test3_TempMatrix.getMatrixElement(1,2)-(-1))<lowEdge);

	EXPECT_EQ(result, fabs(test3_TempMatrix.getMatrixElement(2,0)-(2))<lowEdge);
	EXPECT_EQ(result, fabs(test3_TempMatrix.getMatrixElement(2,1)-(-4))<lowEdge);
	EXPECT_EQ(result, fabs(test3_TempMatrix.getMatrixElement(2,2)-(5))<lowEdge);
}

