/*
 * test_HessenbergTriangleFormular.cpp
 *
 *  Created on: 2017年6月13日
 *      Author: looke
 */

#include "..\gtest_src\gtest\gtest.h"
#include "StaticVector.h"
#include "StaticMatrix.h"
#include "MatrixMultiplier.h"
#include "QRDecomposition.h"
#include "HessenbergTriangleFormular.h"
#include "math.h"
#include <iostream>
using namespace std;

/*
 * 测试HessenbergTriangleFormular 普通4x4矩阵上Hessenberg-Triangle转换
 */
TEST(HessenbergTriangleFormular_HT_Test_Normal4x4, postive)
{
	double lowEdge, maxDiff;
	StaticMatrix AMatrix_44 = StaticMatrix(4,4);
	AMatrix_44.setMatrixElement(0,0,10);
	AMatrix_44.setMatrixElement(0,1,1);
	AMatrix_44.setMatrixElement(0,2,2);
	AMatrix_44.setMatrixElement(0,3,2);

	AMatrix_44.setMatrixElement(1,0,1);
	AMatrix_44.setMatrixElement(1,1,2);
	AMatrix_44.setMatrixElement(1,2,-1);
	AMatrix_44.setMatrixElement(1,3,-1);

	AMatrix_44.setMatrixElement(2,0,1);
	AMatrix_44.setMatrixElement(2,1,1);
	AMatrix_44.setMatrixElement(2,2,2);
	AMatrix_44.setMatrixElement(2,3,2);

	AMatrix_44.setMatrixElement(3,0,1);
	AMatrix_44.setMatrixElement(3,1,1);
	AMatrix_44.setMatrixElement(3,2,2);
	AMatrix_44.setMatrixElement(3,3,2);
	StaticMatrix AMatrix_44_Original = StaticMatrix(4,4);
	AMatrix_44_Original.copyMatrixElementNoCheck(&AMatrix_44);

	StaticMatrix BMatrix_44 = StaticMatrix(4,4);
	StaticMatrix BMatrix_44_Original = StaticMatrix(4,4);
	BMatrix_44_Original.copyMatrixElementNoCheck(&BMatrix_44);

	StaticMatrix QMatrix_44 = StaticMatrix(4,4);
	StaticMatrix ZMatrix_44 = StaticMatrix(4,4);

	StaticMatrix QZMatrix_Step_44 = StaticMatrix(4,4);
	StaticMatrix TempMatrix_Trans_44 = StaticMatrix(4,4);
	StaticMatrix TempMatrix_44 = StaticMatrix(4,4);

	HessenbergTriangleFormular HTFormular = HessenbergTriangleFormular(&AMatrix_44, &BMatrix_44, &QMatrix_44, &ZMatrix_44, &QZMatrix_Step_44, &TempMatrix_Trans_44, &TempMatrix_44);
	HTFormular.formularABMatrix();

	lowEdge = AMatrix_44.getLowEdge();
	EXPECT_LT(fabs(AMatrix_44.getMatrixElement(2,0)-(0)),lowEdge);
	EXPECT_LT(fabs(AMatrix_44.getMatrixElement(3,0)-(0)),lowEdge);
	EXPECT_LT(fabs(AMatrix_44.getMatrixElement(3,1)-(0)),lowEdge);

	lowEdge = BMatrix_44.getLowEdge();
	EXPECT_LT(fabs(BMatrix_44.getMatrixElement(1,0)-(0)),lowEdge);
	EXPECT_LT(fabs(BMatrix_44.getMatrixElement(2,0)-(0)),lowEdge);
	EXPECT_LT(fabs(BMatrix_44.getMatrixElement(3,0)-(0)),lowEdge);

	EXPECT_LT(fabs(BMatrix_44.getMatrixElement(2,1)-(0)),lowEdge);
	EXPECT_LT(fabs(BMatrix_44.getMatrixElement(3,1)-(0)),lowEdge);

	EXPECT_LT(fabs(BMatrix_44.getMatrixElement(3,2)-(0)),lowEdge);

	cout << "Q" << endl;
	QMatrix_44.printMatrix();

	cout << "Z" << endl;
	ZMatrix_44.printMatrix();

	cout << "A" << endl;
	AMatrix_44_Original.printMatrix();

	MatrixMultiplier m_multi = MatrixMultiplier(&QMatrix_44, &AMatrix_44_Original, &TempMatrix_44);
	m_multi.multiplyCalc();
	AMatrix_44_Original.copyMatrixElementNoCheck(&TempMatrix_44);
	m_multi.reload(&AMatrix_44_Original, &ZMatrix_44, &TempMatrix_44);
	m_multi.multiplyCalc();
	AMatrix_44_Original.copyMatrixElementNoCheck(&TempMatrix_44);
	cout << "Q * A * Z" << endl;
	AMatrix_44_Original.printMatrix();

	lowEdge = AMatrix_44_Original.getLowEdge();
	maxDiff = AMatrix_44_Original.calcMaxDifferentialNoCheck(&AMatrix_44);
	EXPECT_LT(maxDiff,lowEdge);
}
