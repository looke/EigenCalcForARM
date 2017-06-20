/*
 * test_GeneralizedEigenSolverForReal.cpp
 *
 *  Created on: 2017年6月18日
 *      Author: looke
 */

#include "..\gtest_src\gtest\gtest.h"
#include "StaticVector.h"
#include "StaticMatrix.h"
#include "MatrixMultiplier.h"
#include "MatrixTransposer.h"
#include "QRDecomposition.h"
#include "GeneralizedEigenSolverForReal.h"
#include "math.h"
#include <iostream>
using namespace std;

/*
 * 测试GeneralizedEigenSolverForReal 求解特征值
 */
TEST(GeneralizedEigenSolverForReal_Test_Normal, postive)
{
	StaticMatrix test44_A = StaticMatrix(4,4);

	test44_A.setMatrixElement(0,0,1);
	test44_A.setMatrixElement(0,1,2);
	test44_A.setMatrixElement(0,2,3);
	test44_A.setMatrixElement(0,3,4);

	test44_A.setMatrixElement(1,0,5);
	test44_A.setMatrixElement(1,1,6);
	test44_A.setMatrixElement(1,2,7);
	test44_A.setMatrixElement(1,3,8);

	test44_A.setMatrixElement(2,0,0);
	test44_A.setMatrixElement(2,1,9);
	test44_A.setMatrixElement(2,2,10);
	test44_A.setMatrixElement(2,3,11);

	test44_A.setMatrixElement(3,0,0);
	test44_A.setMatrixElement(3,1,0);
	test44_A.setMatrixElement(3,2,12);
	test44_A.setMatrixElement(3,3,13);


	StaticMatrix test44_B = StaticMatrix(4,4);

	test44_B.setMatrixElement(0,0,1);
	test44_B.setMatrixElement(0,1,2);
	test44_B.setMatrixElement(0,2,3);
	test44_B.setMatrixElement(0,3,4);

	test44_B.setMatrixElement(1,0,0);
	test44_B.setMatrixElement(1,1,5);
	test44_B.setMatrixElement(1,2,6);
	test44_B.setMatrixElement(1,3,7);

	test44_B.setMatrixElement(2,0,0);
	test44_B.setMatrixElement(2,1,0);
	test44_B.setMatrixElement(2,2,8);
	test44_B.setMatrixElement(2,3,9);

	test44_B.setMatrixElement(3,0,0);
	test44_B.setMatrixElement(3,1,0);
	test44_B.setMatrixElement(3,2,0);
	test44_B.setMatrixElement(3,3,10);

	StaticVector test44_Vector = StaticVector(4);
	StaticMatrix test44_A_Deflate = StaticMatrix(4,4);
	StaticMatrix test44_B_Deflate = StaticMatrix(4,4);

	StaticMatrix test44_Q_Total = StaticMatrix(4,4);
	StaticMatrix test44_Z_Total = StaticMatrix(4,4);

	StaticMatrix test44_Q_Step = StaticMatrix(4,4);
	StaticMatrix test44_Z_Step = StaticMatrix(4,4);
	StaticMatrix test44_QZ_Step = StaticMatrix(4,4);

	StaticMatrix test44_Temp_Trans = StaticMatrix(4,4);
	StaticMatrix test44_Temp = StaticMatrix(4,4);

	GeneralizedEigenSolverForReal gEigenCalc = GeneralizedEigenSolverForReal(
			&test44_A,
			&test44_B,
			&test44_Vector,
			&test44_A_Deflate,
			&test44_B_Deflate,
			&test44_Q_Total,
			&test44_Z_Total,
			&test44_Q_Step,
			&test44_Z_Step,
			&test44_QZ_Step,
			&test44_Temp_Trans,
			&test44_Temp);

	//gEigenCalc.calcEigenValue();

	//EXPECT_EQ(-999, 1);
}
