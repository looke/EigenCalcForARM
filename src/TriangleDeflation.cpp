/*
 * TriangleDeflation.cpp
 *
 *  Created on: 2017年6月22日
 *      Author: looke
 */
#include "TriangleDeflation.h"

TriangleDeflation::TriangleDeflation()
{
	this->deflateStart_new = 0;
	this->deflateEnd_new = 0;
};

/*
 * 查找并更新降阶点
 * 在全尺寸伴随Triangle矩阵中，从右下到左上开始，主对角线上查找0元
 * 找到0元后，返回所在索引
 */
bool TriangleDeflation::findNewDeflationPoint(BasicMatrix* p_TriangleMatrix, int deflationStart, int deflationEnd)
{
	bool result = false;
	double temp;
	int deflatedSize; //降阶点右下角矩阵维度
	this->deflateStart_new = deflationStart;
	this->deflateEnd_new = deflationEnd;
	int i=deflationEnd;


	return result;
};


int TriangleDeflation::getNewDeflationStart()
{
	return this->deflateStart_new;
};
int TriangleDeflation::getNewDeflationEnd()
{
	return this->deflateEnd_new;
};
