/*
 * MatrixTransposer.cpp
 *
 *  Created on: 2017年2月25日
 *      Author: looke
 */

#include "..\include\matrixop\basic\MatrixTransposer.h"
using namespace std;

/*
 * 转置方阵
 */
void MatrixTransposer::transposeSquareMatrix(BasicMatrix* input_matrix) throw (length_error)
{
	int rowNumber = input_matrix->rowNum;
	int columnNumber = input_matrix->columnNum;

	if(rowNumber != columnNumber)
	{
		throw length_error(getMatrixNotSquareErrorMessage(rowNumber,columnNumber));
	}

	double temp;
	for(int i=0; i<rowNumber; i++)
	{
		for(int j=i+1; j<columnNumber; j++)
		{
			temp = input_matrix->getMatrixElement(i,j);
			input_matrix->setMatrixElement(i,j,input_matrix->getMatrixElement(j,i));
			input_matrix->setMatrixElement(j,i,temp);
		}
	}
};

//转置任意矩阵
void MatrixTransposer::transposeMatrix(BasicMatrix* p_opMatrix, BasicMatrix* p_resultMatrix) throw(length_error)
{
	if(p_opMatrix->rowNum != p_resultMatrix->columnNum || p_opMatrix->columnNum != p_resultMatrix->rowNum)
	{
		throw length_error(getMatrixNotSquareErrorMessage(p_opMatrix->rowNum,p_opMatrix->columnNum,p_resultMatrix->rowNum,p_resultMatrix->columnNum));
	}

	for(int i=0; i<p_opMatrix->rowNum; i++)
	{
		for(int j=0; j<p_opMatrix->columnNum ; j++)
		{
			p_resultMatrix->setMatrixElement(j,i,p_opMatrix->getMatrixElement(i,j));
		}
	}

};

/*
 * 生成非方阵异常信息
 */
string MatrixTransposer::getMatrixNotSquareErrorMessage(int rowNumber, int columnNumber)
{
	stringstream stream;

	stream<<rowNumber;
	string row = stream.str();
	stream.str("");
	stream<<columnNumber;
	string column = stream.str();

	string exInfo = "Matrix must be Square. Matrix Row: " +row + " Matrix Column: " + column + ".";
	return exInfo;
};

/*
 * 生成矩阵维度异常信息
 */
string MatrixTransposer::getMatrixNotSquareErrorMessage(int op_rowNumber, int op_columnNumber,int re_rowNumber, int re_columnNumber)
{
	stringstream stream;

	stream<<op_rowNumber;
	string op_row = stream.str();
	stream.str("");
	stream<<op_columnNumber;
	string op_column = stream.str();
	stream.str("");
	stream<<re_rowNumber;
	string re_row = stream.str();
	stream.str("");
	stream<<re_columnNumber;
	string re_column = stream.str();

	string exInfo = "Input Matrix and Result Matrix row and column must be match. Input Matrix Row: " +op_row + " Input Matrix Column: " + op_column + ". Result Matrix Row: " + re_row + " Result Matrix Column: " + re_column + ".";
	return exInfo;
};
