/*
 * test_DoubleShiftQZIteration.cpp
 *
 *  Created on: 2017��6��15��
 *      Author: looke
 */

#include "..\gtest_src\gtest\gtest.h"
#include "StaticVector.h"
#include "StaticMatrix.h"
#include "MatrixMultiplier.h"
#include "MatrixTransposer.h"
#include "QRDecomposition.h"
#include "DoubleShiftQZIteration.h"
#include "math.h"
#include <iostream>
using namespace std;

/*
 * ����DoubleShiftQZIteration ��ʽ����
 */
TEST(DoubleShiftQZIterationIMITTest_Normal4x4, postive)
{
	double lowEdge;
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

	StaticMatrix test44_A_Original = StaticMatrix(4,4);
	test44_A_Original.copyMatrixElementNoCheck(&test44_A);

	StaticVector test44_Vector = StaticVector(3);
	StaticMatrix test44_Q_Total = StaticMatrix(4,4);
	StaticMatrix test44_Z_Total = StaticMatrix(4,4);
	StaticMatrix test44_QZ_Step = StaticMatrix(4,4);
	StaticMatrix test44_Temp_Trans = StaticMatrix(4,4);
	StaticMatrix test44_Temp = StaticMatrix(4,4);

	DoubleShiftQZIteration doubleQZ = DoubleShiftQZIteration(&test44_A, &test44_B, &test44_Vector, &test44_Q_Total, &test44_Z_Total, &test44_QZ_Step, &test44_Temp_Trans, &test44_Temp);
	doubleQZ.wilkinson_IM_QZIteration();
	lowEdge = test44_A.getLowEdge();
	EXPECT_GT(0.0000001, fabs(test44_A.getMatrixElement(0,0)/test44_B.getMatrixElement(0,0)-1));
}

/*
 * ����DoubleShiftQZIteration ��ʽ����
 */
TEST(DoubleShiftQZIterationIMITTest_Normal3x3, postive)
{
	double lowEdge;
	StaticMatrix test33_A = StaticMatrix(3,3);
	test33_A.setMatrixElement(0,0,274990.04191798158);
	test33_A.setMatrixElement(0,1,129700425.46824622);
	test33_A.setMatrixElement(0,2,27251931.289327145);

	test33_A.setMatrixElement(1,0,328007.87061981857);
	test33_A.setMatrixElement(1,1,142016533.00940323);
	test33_A.setMatrixElement(1,2,37156463.745292664);

	test33_A.setMatrixElement(2,0,0);
	test33_A.setMatrixElement(2,1,11957950.451266289);
	test33_A.setMatrixElement(2,2,28387678.2265172);

	StaticMatrix test33_B = StaticMatrix(3,3);
	test33_B.setMatrixElement(0,0,-0.000014789493880935526);
	test33_B.setMatrixElement(0,1,0.0013702168360741711);
	test33_B.setMatrixElement(0,2,0.0036075074094502123);

	test33_B.setMatrixElement(1,0,0);
	test33_B.setMatrixElement(1,1,-0.0080263383475187144);
	test33_B.setMatrixElement(1,2,-0.0047146282466187578);

	test33_B.setMatrixElement(2,0,0);
	test33_B.setMatrixElement(2,1,0);
	test33_B.setMatrixElement(2,2,-0.0018070479197437753);

	StaticMatrix test33_A_Original = StaticMatrix(3,3);
	test33_A_Original.copyMatrixElementNoCheck(&test33_A);

	StaticVector test33_Vector = StaticVector(3);
	StaticMatrix test33_Q_Total = StaticMatrix(3,3);
	StaticMatrix test33_Z_Total = StaticMatrix(3,3);
	StaticMatrix test33_QZ_Step = StaticMatrix(3,3);
	StaticMatrix test33_Temp_Trans = StaticMatrix(3,3);
	StaticMatrix test33_Temp = StaticMatrix(3,3);

	DoubleShiftQZIteration doubleQZ = DoubleShiftQZIteration(&test33_A, &test33_B, &test33_Vector, &test33_Q_Total, &test33_Z_Total, &test33_QZ_Step, &test33_Temp_Trans, &test33_Temp);
	//doubleQZ.wilkinson_IM_QZIteration();
	//lowEdge = test33_A.getLowEdge();
	//EXPECT_GT(0.0000001, fabs(test33_A.getMatrixElement(0,0)/test33_B.getMatrixElement(0,0)-1));
}
