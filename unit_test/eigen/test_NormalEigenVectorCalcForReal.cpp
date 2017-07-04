/*
 * test_NormalEigenVectorCalcForReal.cpp
 *
 *  Created on: 2017年6月20日
 *      Author: looke
 */

#include "..\gtest_src\gtest\gtest.h"
#include "StaticVector.h"
#include "StaticMatrix.h"
#include "MatrixMultiplier.h"
#include "MatrixTransposer.h"
#include "NormalEigenVectorCalcForReal.h"
#include "math.h"
#include <iostream>
using namespace std;

/*
 * 测试NormalEigenVectorCalcForReal 求解特征向量
 */
TEST(NormalEigenVectorCalcForReal_Test_Normal2x2, postive)
{
	StaticMatrix test22 = StaticMatrix(2,2);
	test22.setMatrixElement(0,0,1);
	test22.setMatrixElement(0,1,3);

	test22.setMatrixElement(1,0,0);
	test22.setMatrixElement(1,1,2);

	StaticMatrix test22_Sub = StaticMatrix(2,2);
	StaticMatrix test22_Sub_Inv = StaticMatrix(2,2);
	StaticMatrix test22_Sub_Temp = StaticMatrix(2,2);

	StaticVector eigenVector = StaticVector(2);
	NormalEigenVectorCalcForReal vectorCalc = NormalEigenVectorCalcForReal(&test22, &test22_Sub, &test22_Sub_Inv, &test22_Sub_Temp);
	vectorCalc.getEigenVector(1,&eigenVector);
	//eigenVector.printVector();
	eigenVector.normalizationVector();
	EXPECT_GT(0.0001, fabs(eigenVector.getElement(0) - 0.9487));
	EXPECT_GT(0.0001, fabs(eigenVector.getElement(1) - 0.3162));
	//EXPECT_EQ(3, eigenVector.getElement(0));
	//EXPECT_EQ(1, eigenVector.getElement(1));

	vectorCalc.getEigenVector(0,&eigenVector);
	//eigenVector.printVector();
	eigenVector.normalizationVector();
	EXPECT_EQ(1, eigenVector.getElement(0));
	EXPECT_EQ(0, eigenVector.getElement(1));
}

/*
 * 测试NormalEigenVectorCalcForReal 求解特征向量
 */
TEST(NormalEigenVectorCalcForReal_Test_Normal4x4, postive)
{
	StaticMatrix test44 = StaticMatrix(4,4);
	test44.setMatrixElement(0,0,1);
	test44.setMatrixElement(0,1,3);
	test44.setMatrixElement(0,2,5);
	test44.setMatrixElement(0,3,6);

	test44.setMatrixElement(1,0,0);
	test44.setMatrixElement(1,1,2);
	test44.setMatrixElement(1,2,4);
	test44.setMatrixElement(1,3,5);

	test44.setMatrixElement(2,0,0);
	test44.setMatrixElement(2,1,0);
	test44.setMatrixElement(2,2,11);
	test44.setMatrixElement(2,3,12);

	test44.setMatrixElement(3,0,0);
	test44.setMatrixElement(3,1,0);
	test44.setMatrixElement(3,2,0);
	test44.setMatrixElement(3,3,13);

	StaticMatrix test44_Sub = StaticMatrix(4,4);
	StaticMatrix test44_Sub_Inv = StaticMatrix(4,4);
	StaticMatrix test44_Sub_Temp = StaticMatrix(4,4);

	StaticVector eigenVector = StaticVector(4);
	NormalEigenVectorCalcForReal vectorCalc = NormalEigenVectorCalcForReal(&test44, &test44_Sub, &test44_Sub_Inv, &test44_Sub_Temp);
	vectorCalc.getEigenVector(3,&eigenVector);

	//eigenVector.printVector();

	//eigenVector.printVector();
	eigenVector.normalizationVector();
	EXPECT_GT(0.0001, fabs(eigenVector.getElement(0) - 0.4832));
	EXPECT_GT(0.0001, fabs(eigenVector.getElement(1) - 0.3482));
	EXPECT_GT(0.0001, fabs(eigenVector.getElement(2) - 0.7924));
	EXPECT_GT(0.0001, fabs(eigenVector.getElement(3) - 0.1321));
}

/*
 * 测试NormalEigenVectorCalcForReal 求解2重特征向量
 */
TEST(NormalEigenVectorCalcForReal_Test_Normal2x2_SameEigenValue, postive)
{
	StaticMatrix test22 = StaticMatrix(2,2);
	test22.setMatrixElement(0,0,1);
	test22.setMatrixElement(0,1,3);

	test22.setMatrixElement(1,0,0);
	test22.setMatrixElement(1,1,1);

	StaticMatrix test22_Sub = StaticMatrix(2,2);
	StaticMatrix test22_Sub_Inv = StaticMatrix(2,2);
	StaticMatrix test22_Sub_Temp = StaticMatrix(2,2);

	StaticVector eigenVector = StaticVector(2);
	NormalEigenVectorCalcForReal vectorCalc = NormalEigenVectorCalcForReal(&test22, &test22_Sub, &test22_Sub_Inv, &test22_Sub_Temp);
	vectorCalc.getEigenVector(1,&eigenVector);
	//eigenVector.printVector();
	eigenVector.normalizationVector();
	EXPECT_EQ(1, eigenVector.getElement(0));
	EXPECT_EQ(0, eigenVector.getElement(1));

	vectorCalc.getEigenVector(0,&eigenVector);
	//eigenVector.printVector();
	eigenVector.normalizationVector();
	EXPECT_EQ(1, eigenVector.getElement(0));
	EXPECT_EQ(0, eigenVector.getElement(1));
}

/*
 * 测试NormalEigenVectorCalcForReal 求解2重特征向量
 */
TEST(NormalEigenVectorCalcForReal_Test_Normal3x3_SameEigenValue, postive)
{
	StaticMatrix test33 = StaticMatrix(3,3);
	test33.setMatrixElement(0,0,2);
	test33.setMatrixElement(0,1,1);
	test33.setMatrixElement(0,2,1);

	test33.setMatrixElement(1,0,0);
	test33.setMatrixElement(1,1,1);
	test33.setMatrixElement(1,2,1);

	test33.setMatrixElement(2,0,0);
	test33.setMatrixElement(2,1,0);
	test33.setMatrixElement(2,2,1);

	StaticMatrix test33_Sub = StaticMatrix(3,3);
	StaticMatrix test33_Sub_Inv = StaticMatrix(3,3);
	StaticMatrix test33_Sub_Temp = StaticMatrix(3,3);

	StaticVector eigenVector = StaticVector(3);
	NormalEigenVectorCalcForReal vectorCalc = NormalEigenVectorCalcForReal(&test33, &test33_Sub, &test33_Sub_Inv, &test33_Sub_Temp);
	vectorCalc.getEigenVector(0,&eigenVector);
	//eigenVector.printVector();
	eigenVector.normalizationVector();
	EXPECT_EQ(1, eigenVector.getElement(0));
	EXPECT_EQ(0, eigenVector.getElement(1));
	EXPECT_EQ(0, eigenVector.getElement(2));

	vectorCalc.getEigenVector(1,&eigenVector);
	//eigenVector.printVector();
	eigenVector.normalizationVector();
	EXPECT_GT(0.0001, fabs(eigenVector.getElement(0) + 0.7071));
	EXPECT_GT(0.0001, fabs(eigenVector.getElement(1) +-0.7071));
	EXPECT_EQ(0, eigenVector.getElement(2));

	vectorCalc.getEigenVector(2,&eigenVector);
	//eigenVector.printVector();
	eigenVector.normalizationVector();
	EXPECT_GT(0.0001, fabs(eigenVector.getElement(0) + 0.7071));
	EXPECT_GT(0.0001, fabs(eigenVector.getElement(1) +-0.7071));
	EXPECT_EQ(0, eigenVector.getElement(2));

}

/*
 * 测试NormalEigenVectorCalcForReal 求带有复数特征值矩阵的实数解特征向量
 */
TEST(NormalEigenVectorCalcForReal_Test_Complex3x3, postive)
{
	StaticMatrix test33 = StaticMatrix(3,3);
	test33.setMatrixElement(0,0,1);
	test33.setMatrixElement(0,1,-1);
	test33.setMatrixElement(0,2,1);

	test33.setMatrixElement(1,0,1);
	test33.setMatrixElement(1,1,1);
	test33.setMatrixElement(1,2,1);

	test33.setMatrixElement(2,0,0);
	test33.setMatrixElement(2,1,0);
	test33.setMatrixElement(2,2,2);


	StaticMatrix test33_Sub = StaticMatrix(3,3);
	StaticMatrix test33_Sub_Inv = StaticMatrix(3,3);
	StaticMatrix test33_Sub_Temp = StaticMatrix(3,3);

	StaticVector eigenVector = StaticVector(3);
	NormalEigenVectorCalcForReal vectorCalc = NormalEigenVectorCalcForReal(&test33, &test33_Sub, &test33_Sub_Inv, &test33_Sub_Temp);
	vectorCalc.getEigenVector(2,&eigenVector);

	//eigenVector.printVector();
	eigenVector.normalizationVector();
	EXPECT_GT(0.0001, fabs(eigenVector.getElement(0) - 0));
	EXPECT_GT(0.0001, fabs(eigenVector.getElement(1) - 0.7071));
	EXPECT_GT(0.0001, fabs(eigenVector.getElement(2) - 0.7071));
}

/*
 * 测试NormalEigenVectorCalcForReal 求带有复数特征值矩阵的实数解特征向量
 */
TEST(NormalEigenVectorCalcForReal_Test_Complex4x4, postive)
{
	StaticMatrix test44 = StaticMatrix(4,4);
	test44.setMatrixElement(0,0,1);
	test44.setMatrixElement(0,1,1);
	test44.setMatrixElement(0,2,1);
	test44.setMatrixElement(0,3,1);

	test44.setMatrixElement(1,0,0);
	test44.setMatrixElement(1,1,1);
	test44.setMatrixElement(1,2,-1);
	test44.setMatrixElement(1,3,1);

	test44.setMatrixElement(2,0,0);
	test44.setMatrixElement(2,1,1);
	test44.setMatrixElement(2,2,1);
	test44.setMatrixElement(2,3,1);

	test44.setMatrixElement(3,0,0);
	test44.setMatrixElement(3,1,0);
	test44.setMatrixElement(3,2,0);
	test44.setMatrixElement(3,3,2);
	//test44.printMatrix();
	StaticMatrix test44_Sub = StaticMatrix(4,4);
	StaticMatrix test44_Sub_Inv = StaticMatrix(4,4);
	StaticMatrix test44_Sub_Temp = StaticMatrix(4,4);

	StaticVector eigenVector = StaticVector(4);
	NormalEigenVectorCalcForReal vectorCalc = NormalEigenVectorCalcForReal(&test44, &test44_Sub, &test44_Sub_Inv, &test44_Sub_Temp);
	vectorCalc.getEigenVector(3,&eigenVector);

	//eigenVector.printVector();
	eigenVector.normalizationVector();
	EXPECT_GT(0.0001, fabs(eigenVector.getElement(0) - 0.8165));
	EXPECT_GT(0.0001, fabs(eigenVector.getElement(1) - 0));
	EXPECT_GT(0.0001, fabs(eigenVector.getElement(2) - 0.4082));
	EXPECT_GT(0.0001, fabs(eigenVector.getElement(2) - 0.4082));
}
