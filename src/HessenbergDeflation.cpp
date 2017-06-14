/*
 * HessenbergDeflation.cpp
 *
 *  Created on: 2017��6��1��
 *      Author: looke
 */

#include "HessenbergDeflation.h"

HessenbergDeflation::HessenbergDeflation()
{
	this->deflateStart_new = 0;
	this->deflateEnd_new = 0;
};

/*
 * ���Ҳ����½��׵�
 * ��ȫ�ߴ����Hessenberg�����У������µ����Ͽ�ʼ�������Խ����ϲ���0Ԫ
 * �ҵ�0Ԫ�󣬷���0Ԫ���½ǵ��Ӿ����Ƿ�ά�ȴ���2���Ƿ���Խ���double-shift������
 * ������ܵ��������½��׽�����ʶλ������Ѱ���µ�0Ԫ
 * ������Ե��������½��׿�ʼ��ʶλ������
 */
bool HessenbergDeflation::findNewDeflationPoint(BasicMatrix* p_HessenbergMatrix, int deflationStart, int deflationEnd)
{
	bool result = false;
	double temp;
	int deflatedSize; //���׵����½Ǿ���ά��
	this->deflateStart_new = deflationStart;
	this->deflateEnd_new = deflationEnd;
	int i=deflationEnd;
	while(i > 0)
	{
		//��ѯ�ζԽ���Ԫ��
		temp = p_HessenbergMatrix->getMatrixElement(i,i-1);
		if(temp == 0)
		{
			//�ζԽ���Ԫ��Ϊ0���������׵�
			deflatedSize = this->deflateEnd_new-i+1;

			if(deflatedSize>2)
			{
				//������׵����·��Ӿ���ά�ȴ���2������Խ��е�������
				this->deflateStart_new = i;
				//this->deflateEnd_new = deflationEnd;
				return true;
			}
			else
			{
				//������׵����·��Ӿ���ά��С�ڵ���2
				//���½��׿�ʼ���������ʶ
				this->deflateEnd_new = i - 1;
				//this->deflateStart_new = 0;
			}
			result = true;
		}
		i--;
	}

	return result;
};

int HessenbergDeflation::getNewDeflationStart()
{
	return this->deflateStart_new;
};
int HessenbergDeflation::getNewDeflationEnd()
{
	return this->deflateEnd_new;
};
