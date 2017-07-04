/*
 * test_GeneralizedEigenVectorCalcForReal.cpp
 *
 *  Created on: 2017年6月22日
 *      Author: looke
 */
#include "..\gtest_src\gtest\gtest.h"
#include "StaticVector.h"
#include "StaticMatrix.h"
#include "MatrixMultiplier.h"
#include "MatrixTransposer.h"
#include "MatrixInverser.h"
#include "GeneralizedEigenVectorCalcForReal.h"
#include "GeneralizedEigenSolverForReal.h"

#include "math.h"
#include <iostream>
using namespace std;


/*
 * 测试GeneralizedEigenVectorCalcForReal 求解特征向量
 */
TEST(GeneralizedEigenVectorCalcForReal_Test_Normal2x2, postive)
{
	StaticMatrix test22_A = StaticMatrix(2,2);
	test22_A.setMatrixElement(0,0,1);
	test22_A.setMatrixElement(0,1,3);

	test22_A.setMatrixElement(1,0,0);
	test22_A.setMatrixElement(1,1,2);

	StaticMatrix test22_B = StaticMatrix(2,2);
	test22_B.setMatrixElement(0,0,1);
	test22_B.setMatrixElement(0,1,1);

	test22_B.setMatrixElement(1,0,0);
	test22_B.setMatrixElement(1,1,1);

	StaticMatrix test22_BinvA = StaticMatrix(2,2);
	StaticMatrix test22_Sub = StaticMatrix(2,2);
	StaticMatrix test22_Sub_Inv = StaticMatrix(2,2);
	StaticMatrix test22_Sub_Temp = StaticMatrix(2,2);

	StaticVector resultVector = StaticVector(2);

	GeneralizedEigenVectorCalcForReal generalizedVecotrCalc = GeneralizedEigenVectorCalcForReal(&test22_A,&test22_B,&test22_BinvA,&test22_Sub,&test22_Sub_Inv,&test22_Sub_Temp);
	generalizedVecotrCalc.getEigenVector(1,&resultVector);
	//resultVector.printVector();
}



/*
 * 测试GeneralizedEigenVectorCalcForReal 求解特征向量 B为奇异矩阵
 */
TEST(GeneralizedEigenVectorCalcForReal_Test_Singluar_B_4x4, postive)
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
	test44_B.setMatrixElement(0,0,0);
	test44_B.setMatrixElement(0,1,0);
	test44_B.setMatrixElement(0,2,2);
	test44_B.setMatrixElement(0,3,0);

	test44_B.setMatrixElement(1,0,0);
	test44_B.setMatrixElement(1,1,-1);
	test44_B.setMatrixElement(1,2,0);
	test44_B.setMatrixElement(1,3,0);

	test44_B.setMatrixElement(2,0,2);
	test44_B.setMatrixElement(2,1,0);
	test44_B.setMatrixElement(2,2,0);
	test44_B.setMatrixElement(2,3,0);

	test44_B.setMatrixElement(3,0,0);
	test44_B.setMatrixElement(3,1,0);
	test44_B.setMatrixElement(3,2,0);
	test44_B.setMatrixElement(3,3,0);

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

	gEigenCalc.calcEigenValue();

	StaticMatrix test44_BinvA = StaticMatrix(4,4);
	StaticMatrix test44_Sub = StaticMatrix(4,4);
	StaticMatrix test44_Sub_Inv = StaticMatrix(4,4);
	StaticMatrix test44_Sub_Temp = StaticMatrix(4,4);

	StaticVector resultVector = StaticVector(4);

	GeneralizedEigenVectorCalcForReal generalizedVecotrCalc = GeneralizedEigenVectorCalcForReal(&test44_A,&test44_B,&test44_BinvA,&test44_Sub,&test44_Sub_Inv,&test44_Sub_Temp);
	generalizedVecotrCalc.getEigenVector(2,&resultVector);
	//resultVector.printVector();

	StaticMatrix test41_BinvA_Vector = StaticMatrix(4,1);
	test41_BinvA_Vector.setMatrixElement(0,0,resultVector.getElement(0));
	test41_BinvA_Vector.setMatrixElement(1,0,resultVector.getElement(1));
	test41_BinvA_Vector.setMatrixElement(2,0,resultVector.getElement(2));
	test41_BinvA_Vector.setMatrixElement(3,0,resultVector.getElement(3));
	StaticMatrix test41_final_Vector = StaticMatrix(4,1);


	StaticMatrix test44_Z_Total_inv = StaticMatrix(4,4);
	test44_Z_Total_inv.resetMatrixToI();
	cout << "Z:" << endl;
	//test44_Z_Total.printMatrix();
	/*
	MatrixInverser m_Inverser = MatrixInverser(&test44_Z_Total, &test44_Z_Total_inv);
	m_Inverser.generateInverseMatrix();
	cout << "Z^-1:" << endl;
	test44_Z_Total_inv.printMatrix();
	*/
	MatrixMultiplier m_multi = MatrixMultiplier(&test44_Z_Total,&test41_BinvA_Vector,&test41_final_Vector);
	m_multi.multiplyCalc();
	resultVector.setElement(0,test41_final_Vector.getMatrixElement(0,0));
	resultVector.setElement(1,test41_final_Vector.getMatrixElement(1,0));
	resultVector.setElement(2,test41_final_Vector.getMatrixElement(2,0));
	resultVector.setElement(3,test41_final_Vector.getMatrixElement(3,0));
	resultVector.normalizationVector();

	EXPECT_GT(0.000001,fabs(resultVector.getElement(0) + 0.045656847));
	EXPECT_GT(0.000001,fabs(resultVector.getElement(1) + 0.00948163));
	EXPECT_GT(0.000001,fabs(resultVector.getElement(2) + 0.734004));
	EXPECT_GT(0.000001,fabs(resultVector.getElement(3) - 0.677542));
}
