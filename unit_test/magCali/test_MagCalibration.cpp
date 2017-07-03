/*
 * test_MagCalibration.cpp
 *
 *  Created on: 2017年7月1日
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
 * 测试MagCalibration 普通10x10输入矩阵的校正结果
 */
TEST(MagCalibrationTest_Normal_10x10, postive)
{
	StaticMatrix inMatrix1010 = StaticMatrix(10,10);

	//读取文件内容，初始化校正输入矩阵
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
        else if(10 == c)//换行符号，矩阵行数+1
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
        else if(45 == c)//负号
        {
        	isNegative = true;
        }
        else if(32 == c)//空格符号,number字符串不为空的话需要转为int后，保存在矩阵对应的位置
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
    infile.close();             //关闭文件输入流

    MagCalibration magCali = MagCalibration(&inMatrix1010);
    magCali.generateCaliInfo();
    //EXPECT_EQ(1,1);
}
