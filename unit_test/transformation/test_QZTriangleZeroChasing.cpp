/*
 * test_QZTriangleZeroChasing.cpp
 *
 *  Created on: 2017年6月17日
 *      Author: looke
 */

#include "..\gtest_src\gtest\gtest.h"
#include "StaticVector.h"
#include "StaticMatrix.h"
#include "MatrixMultiplier.h"
#include "QZTriangleZeroChasing.h"
#include "math.h"


/*
 * 测试QZTriangleZeroChasing 0元降阶分解
 */
TEST(QZTriangleZeroChasingTest, postive)
{
	StaticMatrix inputMatrix_A = StaticMatrix(4,4);
	inputMatrix_A.setMatrixElement(0,0,1);
	inputMatrix_A.setMatrixElement(0,1,1);
	inputMatrix_A.setMatrixElement(0,2,1);
	inputMatrix_A.setMatrixElement(0,3,1);

	inputMatrix_A.setMatrixElement(1,0,1);
	inputMatrix_A.setMatrixElement(1,1,1);
	inputMatrix_A.setMatrixElement(1,2,1);
	inputMatrix_A.setMatrixElement(1,3,1);

	inputMatrix_A.setMatrixElement(2,0,0);
	inputMatrix_A.setMatrixElement(2,1,1);
	inputMatrix_A.setMatrixElement(2,2,1);
	inputMatrix_A.setMatrixElement(2,3,1);

	inputMatrix_A.setMatrixElement(3,0,0);
	inputMatrix_A.setMatrixElement(3,1,0);
	inputMatrix_A.setMatrixElement(3,2,1);
	inputMatrix_A.setMatrixElement(3,3,1);
	inputMatrix_A.printMatrix();

	StaticMatrix inputMatrix_B = StaticMatrix(4,4);
	inputMatrix_B.setMatrixElement(0,0,1);
	inputMatrix_B.setMatrixElement(0,1,2);
	inputMatrix_B.setMatrixElement(0,2,3);
	inputMatrix_B.setMatrixElement(0,3,4);

	inputMatrix_B.setMatrixElement(1,0,0);
	inputMatrix_B.setMatrixElement(1,1,0);
	inputMatrix_B.setMatrixElement(1,2,6);
	inputMatrix_B.setMatrixElement(1,3,7);

	inputMatrix_B.setMatrixElement(2,0,0);
	inputMatrix_B.setMatrixElement(2,1,0);
	inputMatrix_B.setMatrixElement(2,2,8);
	inputMatrix_B.setMatrixElement(2,3,9);

	inputMatrix_B.setMatrixElement(3,0,0);
	inputMatrix_B.setMatrixElement(3,1,0);
	inputMatrix_B.setMatrixElement(3,2,0);
	inputMatrix_B.setMatrixElement(3,3,10);

	StaticMatrix inputMatrix_Sub_A = StaticMatrix(4,4);
	StaticMatrix inputMatrix_Sub_B = StaticMatrix(4,4);
	StaticMatrix inputMatrix_Q_Total = StaticMatrix(4,4);
	StaticMatrix inputMatrix_Z_Total = StaticMatrix(4,4);
	StaticMatrix inputMatrix_QZ_Step = StaticMatrix(4,4);
	StaticMatrix inputMatrix_TempMatrix = StaticMatrix(4,4);

	QZTriangleZeroChasing qz0Chasing = QZTriangleZeroChasing(&inputMatrix_A, &inputMatrix_B, &inputMatrix_Sub_A, &inputMatrix_Sub_B, &inputMatrix_Q_Total,&inputMatrix_Z_Total,&inputMatrix_QZ_Step,&inputMatrix_TempMatrix);
	qz0Chasing.deflate();

	EXPECT_EQ(-999, 1);

}
