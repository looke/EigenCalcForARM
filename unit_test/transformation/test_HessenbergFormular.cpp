/*
 * test_HessenbergFormular.cpp
 *
 *  Created on: 2017年6月13日
 *      Author: looke
 */
#include "..\gtest_src\gtest\gtest.h"
#include "StaticVector.h"
#include "StaticMatrix.h"
#include "MatrixMultiplier.h"
#include "MatrixTransposer.h"
#include "HessenbergFormular.h"
#include "math.h"


/*
 * 测试HessenbergFormular 普通4x4矩阵上Hessenberg转换
 */
TEST(HessenbergFormularUpperHessenTest_Normal4x4, postive)
{
	double lowEdge;
	//bool result = true;
	MatrixTransposer m_Transposer = MatrixTransposer();

	StaticMatrix test44 = StaticMatrix(4,4);
/*
	test44.setMatrixElement(0,0,54);
	test44.setMatrixElement(0,1,40);
	test44.setMatrixElement(0,2,10);
	test44.setMatrixElement(0,3,76);

	test44.setMatrixElement(1,0,47);
	test44.setMatrixElement(1,1,20);
	test44.setMatrixElement(1,2,94);
	test44.setMatrixElement(1,3,49);

	test44.setMatrixElement(2,0,26);
	test44.setMatrixElement(2,1,80);
	test44.setMatrixElement(2,2,94);
	test44.setMatrixElement(2,3,70);

	test44.setMatrixElement(3,0,3);
	test44.setMatrixElement(3,1,92);
	test44.setMatrixElement(3,2,83);
	test44.setMatrixElement(3,3,45);
*/
	test44.setMatrixElement(0,0,5);
	test44.setMatrixElement(0,1,-2);
	test44.setMatrixElement(0,2,2.8284271);
	test44.setMatrixElement(0,3,-4.2426407);

	test44.setMatrixElement(1,0,1);
	test44.setMatrixElement(1,1,0);
	test44.setMatrixElement(1,2,3.5355339);
	test44.setMatrixElement(1,3,-0.7071068);

	test44.setMatrixElement(2,0,0);
	test44.setMatrixElement(2,1,-1.41421356);
	test44.setMatrixElement(2,2,1);
	test44.setMatrixElement(2,3,0);

	test44.setMatrixElement(3,0,0);
	test44.setMatrixElement(3,1,1.41421356);
	test44.setMatrixElement(3,2,-4);
	test44.setMatrixElement(3,3,-1);

	StaticMatrix test44_Original = StaticMatrix(4,4);
	test44_Original.copyMatrixElementNoCheck(&test44);
	StaticMatrix test44_PreTrans = StaticMatrix(4,4);
	StaticMatrix test44_AfterTrans = StaticMatrix(4,4);
	StaticMatrix test44_Trans = StaticMatrix(4,4);
	StaticMatrix test44_Temp = StaticMatrix(4,4);

	HessenbergFormular hessenForm = HessenbergFormular(&test44, &test44_PreTrans, &test44_Trans, &test44_Temp);
	hessenForm.formularUpperHessnbergMatrix();
	lowEdge = test44.getLowEdge();
	//test44.printMatrix();

	EXPECT_LT(fabs(test44.getMatrixElement(2,0)-(0)),lowEdge);
	EXPECT_LT(fabs(test44.getMatrixElement(3,0)-(0)),lowEdge);
	EXPECT_LT(fabs(test44.getMatrixElement(3,1)-(0)),lowEdge);

	test44_AfterTrans.copyMatrixElementNoCheck(&test44_PreTrans);
	m_Transposer.transposeSquareMatrix(&test44_AfterTrans);
	MatrixMultiplier m_multi = MatrixMultiplier(&test44_AfterTrans,&test44,&test44_Temp);
	m_multi.multiplyCalc();
	test44.copyMatrixElementNoCheck(&test44_Temp);
	m_multi.reload(&test44,&test44_PreTrans,&test44_Temp);
	m_multi.multiplyCalc();
	test44.copyMatrixElementNoCheck(&test44_Temp);

	//test44.printMatrix();
	lowEdge = test44.getLowEdge();
	double maxDiff = test44.calcMaxDifferentialNoCheck(&test44_Original);
	EXPECT_LT(maxDiff, lowEdge);
}

/*
 * 测试HessenbergFormular 上Hessenberg4x4矩阵上Hessenberg转换
 */
TEST(HessenbergFormularUpperHessenTest_UpperHessen4x4, postive)
{
	double lowEdge;
	//bool result = true;
	MatrixTransposer m_Transposer = MatrixTransposer();

	StaticMatrix test44_up = StaticMatrix(4,4);

	test44_up.setMatrixElement(0,0,5);
	test44_up.setMatrixElement(0,1,-2);
	test44_up.setMatrixElement(0,2,2.8284271);
	test44_up.setMatrixElement(0,3,-4.2426407);

	test44_up.setMatrixElement(1,0,1);
	test44_up.setMatrixElement(1,1,0);
	test44_up.setMatrixElement(1,2,3.5355339);
	test44_up.setMatrixElement(1,3,-0.7071068);

	test44_up.setMatrixElement(2,0,0);
	test44_up.setMatrixElement(2,1,-1.41421356);
	test44_up.setMatrixElement(2,2,1);
	test44_up.setMatrixElement(2,3,0);

	test44_up.setMatrixElement(3,0,0);
	test44_up.setMatrixElement(3,1,0);
	test44_up.setMatrixElement(3,2,-4);
	test44_up.setMatrixElement(3,3,-1);

	StaticMatrix test44_Original = StaticMatrix(4,4);
	test44_Original.copyMatrixElementNoCheck(&test44_up);
	StaticMatrix test44_PreTrans = StaticMatrix(4,4);
	StaticMatrix test44_AfterTrans = StaticMatrix(4,4);
	StaticMatrix test44_Trans = StaticMatrix(4,4);
	StaticMatrix test44_Temp = StaticMatrix(4,4);

	HessenbergFormular hessenForm = HessenbergFormular(&test44_up, &test44_PreTrans, &test44_Trans, &test44_Temp);
	hessenForm.formularUpperHessnbergMatrix();
	lowEdge = test44_up.getLowEdge();
	//test44.printMatrix();

	EXPECT_EQ(5,test44_up.getMatrixElement(0,0));
	EXPECT_EQ(-2,test44_up.getMatrixElement(0,1));
	EXPECT_EQ(2.8284271,test44_up.getMatrixElement(0,2));
	EXPECT_EQ(-4.2426407,test44_up.getMatrixElement(0,3));

	EXPECT_EQ(1,test44_up.getMatrixElement(1,0));
	EXPECT_EQ(0,test44_up.getMatrixElement(1,1));
	EXPECT_EQ(3.5355339,test44_up.getMatrixElement(1,2));
	EXPECT_EQ(-0.7071068,test44_up.getMatrixElement(1,3));

	EXPECT_EQ(0,test44_up.getMatrixElement(2,0));
	EXPECT_EQ(-1.41421356,test44_up.getMatrixElement(2,1));
	EXPECT_EQ(1,test44_up.getMatrixElement(2,2));
	EXPECT_EQ(0,test44_up.getMatrixElement(2,3));

	EXPECT_EQ(0,test44_up.getMatrixElement(3,0));
	EXPECT_EQ(0,test44_up.getMatrixElement(3,1));
	EXPECT_EQ(-4,test44_up.getMatrixElement(3,2));
	EXPECT_EQ(-1,test44_up.getMatrixElement(3,3));

	EXPECT_EQ(1,test44_PreTrans.getMatrixElement(0,0));
	EXPECT_EQ(0,test44_PreTrans.getMatrixElement(0,1));
	EXPECT_EQ(0,test44_PreTrans.getMatrixElement(0,2));
	EXPECT_EQ(0,test44_PreTrans.getMatrixElement(0,3));

	EXPECT_EQ(0,test44_PreTrans.getMatrixElement(1,0));
	EXPECT_EQ(1,test44_PreTrans.getMatrixElement(1,1));
	EXPECT_EQ(0,test44_PreTrans.getMatrixElement(1,2));
	EXPECT_EQ(0,test44_PreTrans.getMatrixElement(1,3));

	EXPECT_EQ(0,test44_PreTrans.getMatrixElement(2,0));
	EXPECT_EQ(0,test44_PreTrans.getMatrixElement(2,1));
	EXPECT_EQ(1,test44_PreTrans.getMatrixElement(2,2));
	EXPECT_EQ(0,test44_PreTrans.getMatrixElement(2,3));

	EXPECT_EQ(0,test44_PreTrans.getMatrixElement(3,0));
	EXPECT_EQ(0,test44_PreTrans.getMatrixElement(3,1));
	EXPECT_EQ(0,test44_PreTrans.getMatrixElement(3,2));
	EXPECT_EQ(1,test44_PreTrans.getMatrixElement(3,3));





	test44_AfterTrans.copyMatrixElementNoCheck(&test44_PreTrans);
	m_Transposer.transposeSquareMatrix(&test44_AfterTrans);
	MatrixMultiplier m_multi = MatrixMultiplier(&test44_AfterTrans,&test44_up,&test44_Temp);
	m_multi.multiplyCalc();
	test44_up.copyMatrixElementNoCheck(&test44_Temp);
	m_multi.reload(&test44_up,&test44_PreTrans,&test44_Temp);
	m_multi.multiplyCalc();
	test44_up.copyMatrixElementNoCheck(&test44_Temp);

	//test44.printMatrix();
	lowEdge = test44_up.getLowEdge();
	double maxDiff = test44_up.calcMaxDifferentialNoCheck(&test44_Original);
	EXPECT_LT(maxDiff, lowEdge);
}
