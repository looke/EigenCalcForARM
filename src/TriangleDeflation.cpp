/*
 * TriangleDeflation.cpp
 *
 *  Created on: 2017��6��22��
 *      Author: looke
 */
#include "TriangleDeflation.h"

TriangleDeflation::TriangleDeflation()
{
	this->deflateStart_new = 0;
	this->deflateEnd_new = 0;
};

/*
 * ���Ҳ����½��׵�
 * ��ȫ�ߴ����Triangle�����У������µ����Ͽ�ʼ�����Խ����ϲ���0Ԫ
 * �ҵ�0Ԫ�󣬷�����������
 */
bool TriangleDeflation::findNewDeflationPoint(BasicMatrix* p_TriangleMatrix, int deflationStart, int deflationEnd)
{
	bool result = false;
	double temp;
	int deflatedSize; //���׵����½Ǿ���ά��
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
