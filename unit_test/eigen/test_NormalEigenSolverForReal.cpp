/*
 * test_NormalEigenSolverForReal.cpp
 *
 *  Created on: 2017��6��16��
 *      Author: looke
 */

#include "..\gtest_src\gtest\gtest.h"
#include "StaticVector.h"
#include "StaticMatrix.h"
#include "MatrixMultiplier.h"
#include "MatrixTransposer.h"
#include "QRDecomposition.h"
#include "NormalEigenSolverForReal.h"
#include "math.h"
#include <iostream>
using namespace std;


/*
 * ����NormalEigenSolverForReal �������ֵ
 */
TEST(NormalEigenSolverForReal_Test_Normal6x6, postive)
{
	double lowEdge;
	StaticMatrix test66 = StaticMatrix(6,6);
	test66.setMatrixElement(0,0,1);
	test66.setMatrixElement(0,1,2);
	test66.setMatrixElement(0,2,3);
	test66.setMatrixElement(0,3,4);
	test66.setMatrixElement(0,4,5);
	test66.setMatrixElement(0,5,6);

	test66.setMatrixElement(1,0,7);
	test66.setMatrixElement(1,1,8);
	test66.setMatrixElement(1,2,9);
	test66.setMatrixElement(1,3,10);
	test66.setMatrixElement(1,4,11);
	test66.setMatrixElement(1,5,12);

	test66.setMatrixElement(2,0,0);
	test66.setMatrixElement(2,1,13);
	test66.setMatrixElement(2,2,14);
	test66.setMatrixElement(2,3,15);
	test66.setMatrixElement(2,4,16);
	test66.setMatrixElement(2,5,17);

	test66.setMatrixElement(3,0,0);
	test66.setMatrixElement(3,1,0);
	test66.setMatrixElement(3,2,18);
	test66.setMatrixElement(3,3,19);
	test66.setMatrixElement(3,4,20);
	test66.setMatrixElement(3,5,21);

	test66.setMatrixElement(4,0,0);
	test66.setMatrixElement(4,1,0);
	test66.setMatrixElement(4,2,0);
	test66.setMatrixElement(4,3,22);
	test66.setMatrixElement(4,4,23);
	test66.setMatrixElement(4,5,24);

	test66.setMatrixElement(5,0,0);
	test66.setMatrixElement(5,1,0);
	test66.setMatrixElement(5,2,0);
	test66.setMatrixElement(5,3,0);
	test66.setMatrixElement(5,4,25);
	test66.setMatrixElement(5,5,26);

	StaticMatrix test66_Original = StaticMatrix(6,6);
	test66_Original.copyMatrixElementNoCheck(&test66);

	StaticVector test66_vector = StaticVector(6);
	StaticMatrix test66_QTMatrix = StaticMatrix(6,6);
	StaticMatrix test66_QMatrix = StaticMatrix(6,6);
	StaticMatrix test66_QQTMatrix = StaticMatrix(6,6);
	StaticMatrix test66_DeflatedMatrix = StaticMatrix(6,6);
	StaticMatrix test66_TempMatrix_Trans = StaticMatrix(6,6);
	StaticMatrix test66_TempMatrix = StaticMatrix(6,6);

	NormalEigenSolverForReal normalEigen = NormalEigenSolverForReal(&test66,&test66_vector,&test66_QTMatrix,&test66_QMatrix,&test66_QQTMatrix,&test66_DeflatedMatrix,&test66_TempMatrix_Trans,&test66_TempMatrix);
	normalEigen.calcEigenValue();

	lowEdge = test66.getLowEdge();
	EXPECT_GT(lowEdge,test66.getMatrixElement(1,0));
}
