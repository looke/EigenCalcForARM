/*
 * test_MatrixEquationSolver.cpp
 *
 *  Created on: 2017��7��1��
 *      Author: looke
 */

#include "..\gtest_src\gtest\gtest.h"
#include "StaticMatrix.h"
#include "MatrixEquationSolver.h"
#include "math.h"

/*
 * ����MatrixEquationSolver ��ⳣ�� 4Ԫ ���Է�����
 */
TEST(MatrixEquationSolver_Norm4x5_Test, postive)
{
	StaticMatrix leftMatrix45 = StaticMatrix(4,5);
	leftMatrix45.setMatrixElement(0,0,2);
	leftMatrix45.setMatrixElement(0,1,3);
	leftMatrix45.setMatrixElement(0,2,11);
	leftMatrix45.setMatrixElement(0,3,5);
	leftMatrix45.setMatrixElement(0,4,2);

	leftMatrix45.setMatrixElement(1,0,1);
	leftMatrix45.setMatrixElement(1,1,1);
	leftMatrix45.setMatrixElement(1,2,5);
	leftMatrix45.setMatrixElement(1,3,2);
	leftMatrix45.setMatrixElement(1,4,1);

	leftMatrix45.setMatrixElement(2,0,2);
	leftMatrix45.setMatrixElement(2,1,1);
	leftMatrix45.setMatrixElement(2,2,3);
	leftMatrix45.setMatrixElement(2,3,2);
	leftMatrix45.setMatrixElement(2,4,-3);

	leftMatrix45.setMatrixElement(3,0,1);
	leftMatrix45.setMatrixElement(3,1,1);
	leftMatrix45.setMatrixElement(3,2,3);
	leftMatrix45.setMatrixElement(3,3,4);
	leftMatrix45.setMatrixElement(3,4,-3);

	MatrixEquationSolver m_EqSolver = MatrixEquationSolver(&leftMatrix45);
	bool result = m_EqSolver.solveMatrixEquation();

	EXPECT_EQ(true, result);
	EXPECT_EQ(-2,leftMatrix45.getMatrixElement(0,4));
	EXPECT_EQ(0,leftMatrix45.getMatrixElement(1,4));
	EXPECT_EQ(1,leftMatrix45.getMatrixElement(2,4));
	EXPECT_EQ(-1,leftMatrix45.getMatrixElement(3,4));
}

/*
 * ����MatrixEquationSolver ��� 4Ԫ ���Է����� ���������
 */
TEST(MatrixEquationSolver_Norm4x5_InfinitRoot_Test, postive)
{
	StaticMatrix leftMatrix45 = StaticMatrix(4,5);
	leftMatrix45.setMatrixElement(0,0,2);
	leftMatrix45.setMatrixElement(0,1,3);
	leftMatrix45.setMatrixElement(0,2,11);
	leftMatrix45.setMatrixElement(0,3,5);
	leftMatrix45.setMatrixElement(0,4,2);

	leftMatrix45.setMatrixElement(1,0,1);
	leftMatrix45.setMatrixElement(1,1,1);
	leftMatrix45.setMatrixElement(1,2,5);
	leftMatrix45.setMatrixElement(1,3,2);
	leftMatrix45.setMatrixElement(1,4,1);

	leftMatrix45.setMatrixElement(2,0,2);
	leftMatrix45.setMatrixElement(2,1,1);
	leftMatrix45.setMatrixElement(2,2,3);
	leftMatrix45.setMatrixElement(2,3,2);
	leftMatrix45.setMatrixElement(2,4,-3);

	leftMatrix45.setMatrixElement(3,0,2);
	leftMatrix45.setMatrixElement(3,1,1);
	leftMatrix45.setMatrixElement(3,2,3);
	leftMatrix45.setMatrixElement(3,3,2);
	leftMatrix45.setMatrixElement(3,4,-3);

	MatrixEquationSolver m_EqSolver = MatrixEquationSolver(&leftMatrix45);
	bool result = m_EqSolver.solveMatrixEquation();

	EXPECT_TRUE(result);
}


/*
 * ����MatrixEquationSolver ��� 3Ԫ ���Է����� �޽����
 */
TEST(MatrixEquationSolver_Norm4x4_NoRoot_Test, postive)
{
	StaticMatrix leftMatrix44 = StaticMatrix(4,4);
	leftMatrix44.setMatrixElement(0,0,1);
	leftMatrix44.setMatrixElement(0,1,3);
	leftMatrix44.setMatrixElement(0,2,-2);
	leftMatrix44.setMatrixElement(0,3,4);

	leftMatrix44.setMatrixElement(1,0,3);
	leftMatrix44.setMatrixElement(1,1,2);
	leftMatrix44.setMatrixElement(1,2,-5);
	leftMatrix44.setMatrixElement(1,3,11);

	leftMatrix44.setMatrixElement(2,0,2);
	leftMatrix44.setMatrixElement(2,1,1);
	leftMatrix44.setMatrixElement(2,2,1);
	leftMatrix44.setMatrixElement(2,3,3);

	leftMatrix44.setMatrixElement(3,0,-2);
	leftMatrix44.setMatrixElement(3,1,1);
	leftMatrix44.setMatrixElement(3,2,3);
	leftMatrix44.setMatrixElement(3,3,-6);

	MatrixEquationSolver m_EqSolver = MatrixEquationSolver(&leftMatrix44);
	bool result = m_EqSolver.solveMatrixEquation();

	EXPECT_FALSE(result);
}
