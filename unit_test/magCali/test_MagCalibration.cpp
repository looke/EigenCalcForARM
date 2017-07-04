/*
 * test_MagCalibration.cpp
 *
 *  Created on: 2017��7��1��
 *      Author: looke
 */
#include "..\gtest_src\gtest\gtest.h"
#include "StaticVector.h"
#include "StaticMatrix.h"
#include "MatrixMultiplier.h"
#include "MatrixTransposer.h"
#include "QRDecomposition.h"
#include "GeneralizedEigenVectorCalcForReal.h"
#include "GeneralizedEigenSolverForReal.h"
#include "MagCalibration.h"
#include "math.h"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
using namespace std;

/*
 * ����MagCalibration ��ͨ10x10��������У�����
 */
TEST(MagCalibrationTest_Normal_10x10_100Sample, postive)
{
	StaticMatrix inMatrix1010 = StaticMatrix(10,10);

	//��ȡ�ļ����ݣ���ʼ��У���������
	string file = "D:\\looke\\CppWorkspace\\EigenCalcForARM\\magCaliVerify.txt";

	//string file = "D:\looke\CppWorkspace\EigenCalcForARM\fileReadTest.txt";
	ifstream infile;
	infile.open(file.data());
	EXPECT_TRUE(infile.is_open());


	if (!infile.is_open())
	{
		cout << "Open File Failed!" << endl;
		return;
	}
	char c;
	infile >> noskipws;
	string number="";
	number.clear();
	int rowNum = 0;
	int columnNum = 0;
	stringstream ss;
	bool isNegative = false;
	while (!infile.eof())
    {
        infile>>c;
        if(c<=57 && c>=48)
        {
        	number.push_back(c);
        }
        else if(10 == c)//���з��ţ���������+1
        {
        	if(number.length()>0)
        	{
        		ss << number;
        		double iValue;
        		ss >> iValue;
        		if(isNegative)
        		{
        			iValue = 0-iValue;
        		}
        		inMatrix1010.setMatrixElement(rowNum,columnNum,iValue);
        		columnNum++;
        		number.clear();
        		number = "";
        		ss.clear();
        		ss.str("");
        		isNegative = false;
        	}
        	rowNum++;
        	columnNum = 0;
        }
        else if(45 == c)//����
        {
        	isNegative = true;
        }
        else if(32 == c)//�ո����,number�ַ�����Ϊ�յĻ���ҪתΪint�󣬱����ھ����Ӧ��λ��
        {
        	if(number.length()>0)
        	{
        		ss << number;
        		double iValue;
        		ss >> iValue;
        		if(isNegative)
        		{
        			iValue = 0-iValue;
        		}
        		inMatrix1010.setMatrixElement(rowNum,columnNum,iValue);
        		columnNum++;
        		number.clear();
        		number = "";
        		ss.clear();
        		ss.str("");
        		isNegative = false;
        	}
        }

    }
    infile.close();             //�ر��ļ�������

    //MagCalibration magCali = MagCalibration(&inMatrix1010);
    //magCali.generateCaliInfo();
    //EXPECT_EQ(1,1);
}

/*
 * ����MagCalibration ��ͨ10x10��������У����� 20������
 */
TEST(MagCalibrationTest_Normal_10x10_20Sample, postive)
{
	StaticMatrix inMatrix1010 = StaticMatrix(10,10);

	//��ȡ�ļ����ݣ���ʼ��У���������
	string file = "D:\\looke\\CppWorkspace\\EigenCalcForARM\\magCaliVerify_20_A.txt";

	//string file = "D:\looke\CppWorkspace\EigenCalcForARM\fileReadTest.txt";
	ifstream infile;
	infile.open(file.data());
	EXPECT_TRUE(infile.is_open());


	if (!infile.is_open())
	{
		cout << "Open File Failed!" << endl;
		return;
	}
	char c;
	infile >> noskipws;
	string number="";
	number.clear();
	int rowNum = 0;
	int columnNum = 0;
	stringstream ss;
	bool isNegative = false;
	while (!infile.eof())
    {
        infile>>c;
        if(c<=57 && c>=48)
        {
        	number.push_back(c);
        }
        else if(10 == c)//���з��ţ���������+1
        {
        	if(number.length()>0)
        	{
        		ss << number;
        		double iValue;
        		ss >> iValue;
        		if(isNegative)
        		{
        			iValue = 0-iValue;
        		}
        		inMatrix1010.setMatrixElement(rowNum,columnNum,iValue);
        		columnNum++;
        		number.clear();
        		number = "";
        		ss.clear();
        		ss.str("");
        		isNegative = false;
        	}
        	rowNum++;
        	columnNum = 0;
        }
        else if(45 == c)//����
        {
        	isNegative = true;
        }
        else if(32 == c)//�ո����,number�ַ�����Ϊ�յĻ���ҪתΪint�󣬱����ھ����Ӧ��λ��
        {
        	if(number.length()>0)
        	{
        		ss << number;
        		double iValue;
        		ss >> iValue;
        		if(isNegative)
        		{
        			iValue = 0-iValue;
        		}
        		inMatrix1010.setMatrixElement(rowNum,columnNum,iValue);
        		columnNum++;
        		number.clear();
        		number = "";
        		ss.clear();
        		ss.str("");
        		isNegative = false;
        	}
        }

    }
    infile.close();             //�ر��ļ�������

    MagCalibration magCali = MagCalibration(&inMatrix1010);
    bool result = magCali.generateCaliInfo();

    EXPECT_TRUE(result);
    EXPECT_GT(0.0001, fabs(magCali.mag_X_HardIron_Cali_Para + 370.99684));
    EXPECT_GT(0.0001, fabs(magCali.mag_Y_HardIron_Cali_Para + 178.42237));
    EXPECT_GT(0.0001, fabs(magCali.mag_Z_HardIron_Cali_Para + 383.248701));

    StaticMatrix expectSoftIronCali = StaticMatrix(3,3);
    expectSoftIronCali.setMatrixElement(0,0,0.703350569345957);
    expectSoftIronCali.setMatrixElement(0,1,-0.00564015107245058);
    expectSoftIronCali.setMatrixElement(0,2,-0.0570317258816459);

    expectSoftIronCali.setMatrixElement(1,0,-0.00564015107245058);
    expectSoftIronCali.setMatrixElement(1,1,0.70395596781518);
    expectSoftIronCali.setMatrixElement(1,2,0.0821945828788549);

    expectSoftIronCali.setMatrixElement(2,0,-0.0570317258816459);
    expectSoftIronCali.setMatrixElement(2,1,0.0821945828788549);
    expectSoftIronCali.setMatrixElement(2,2,0.575822153801268);

    double maxDiff = expectSoftIronCali.calcMaxDifferentialNoCheck(&magCali.m_Mag_SoftIron_CaliMatrix);
    EXPECT_GT(0.0001, maxDiff);
}
