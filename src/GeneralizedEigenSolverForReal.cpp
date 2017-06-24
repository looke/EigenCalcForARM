/*
 * GeneralizedEigenSolverForReal.cpp
 *
 *  Created on: 2017��6��1��
 *      Author: looke
 */

#include "GeneralizedEigenSolverForReal.h"
#include "iostream"
using namespace std;

//GeneralizedEigenSolverForReal::GeneralizedEigenSolverForReal()
//{};

GeneralizedEigenSolverForReal::GeneralizedEigenSolverForReal(BasicMatrix* p_input_OpMatrix_A, BasicMatrix* p_input_OpMatrix_B, BasicVector* p_input_Vector, BasicMatrix* p_input_A_deflated, BasicMatrix* p_input_B_deflated, BasicMatrix* p_input_Q_Total, BasicMatrix* p_input_Z_Total, BasicMatrix* p_input_Q_Step, BasicMatrix* p_input_Z_Step, BasicMatrix* p_input_QZ_Step, BasicMatrix* p_input_TempMatrix_Trans,BasicMatrix* p_input_TempMatrix)
:m_ABInvCalc(),m_HessenbergDeflation(),m_Multiplier(p_input_OpMatrix_A,p_input_OpMatrix_A,p_input_TempMatrix),m_testMulti(p_input_OpMatrix_A,p_input_OpMatrix_A,p_input_TempMatrix),
 m_SingleShifeQZ(p_input_OpMatrix_A,p_input_OpMatrix_B,p_input_Q_Step,p_input_Z_Step,p_input_QZ_Step,p_input_TempMatrix_Trans,p_input_TempMatrix),
 m_DoubleShifeQZ(p_input_OpMatrix_A,p_input_OpMatrix_B,p_input_Vector,p_input_Q_Total,p_input_Z_Total,p_input_QZ_Step,p_input_TempMatrix_Trans,p_input_TempMatrix),
 m_HessenbergTriangleFormular(p_input_OpMatrix_A,p_input_OpMatrix_B,p_input_Q_Step,p_input_Z_Step,p_input_QZ_Step,p_input_TempMatrix_Trans,p_input_TempMatrix),
 m_QZTriangleZeroChasing(p_input_OpMatrix_A,p_input_OpMatrix_B,p_input_A_deflated,p_input_B_deflated,p_input_Q_Step, p_input_Z_Step,p_input_QZ_Step, p_input_TempMatrix)
{
	this->init(p_input_OpMatrix_A,
			p_input_OpMatrix_B,
			p_input_Vector,
			p_input_A_deflated,
			p_input_B_deflated,
			p_input_Q_Total,
			p_input_Z_Total,
			p_input_Q_Step,
			p_input_Z_Step,
			p_input_QZ_Step,
			p_input_TempMatrix_Trans,
			p_input_TempMatrix);
};

void GeneralizedEigenSolverForReal::init(BasicMatrix* p_input_OpMatrix_A, BasicMatrix* p_input_OpMatrix_B, BasicVector* p_input_Vector, BasicMatrix* p_input_A_deflated, BasicMatrix* p_input_B_deflated, BasicMatrix* p_input_Q_Total, BasicMatrix* p_input_Z_Total, BasicMatrix* p_input_Q_Step, BasicMatrix* p_input_Z_Step, BasicMatrix* p_input_QZ_Step, BasicMatrix* p_input_TempMatrix_Trans,BasicMatrix* p_input_TempMatrix)
{
	//ԭʼ��������
	this->p_OpMatrix_A = p_input_OpMatrix_A;
	this->p_OpMatrix_B = p_input_OpMatrix_B;

	this->p_OpTransVector = p_input_Vector;

	//ԭʼ�����ȫά������任����
	this->p_QMatrix_Total = p_input_Q_Total;
	this->p_ZMatrix_Total = p_input_Z_Total;

	//ԭʼ����� �����任���� Q\Z
	this->p_QMatrix_Step = p_input_Q_Step;
	this->p_ZMatrix_Step = p_input_Z_Step;

	//�ѽ��׵� ����Hessenberg����
	this->p_OpMatrix_Hessenberg_deflated = p_input_A_deflated;
	//�ѽ��׵� ����Triangle����
	this->p_OpMatrix_Triangle_deflated = p_input_B_deflated;

	//�м���̾���
	this->p_QZMatrix_Step = p_input_QZ_Step;
	//�м���̾���
	this->p_TempMatrix_Trans = p_input_TempMatrix_Trans;
	//�м���̾���
	this->p_TempMatrix = p_input_TempMatrix;

	//�����������ָʾ����ʼ��Ϊ0
	this->deflationStart = 0;
	//�����յ�����ָʾ����ʼ��Ϊn-1
	this->deflationEnd = p_OpMatrix_A->rowNum - 1;

	this->testForTemp_A_nxn = StaticMatrix(p_input_OpMatrix_A->rowNum,p_input_OpMatrix_A->columnNum);
	testForTemp_A_nxn.copyMatrixElementNoCheck(p_OpMatrix_A);

	this->testForTemp_B_nxn = StaticMatrix(p_input_OpMatrix_A->rowNum,p_input_OpMatrix_A->columnNum);
	testForTemp_B_nxn.copyMatrixElementNoCheck(p_OpMatrix_B);

	this->testTemp_nxn = StaticMatrix(p_input_OpMatrix_A->rowNum,p_input_OpMatrix_A->columnNum);
};

void GeneralizedEigenSolverForReal::reload(BasicMatrix* p_input_OpMatrix_A, BasicMatrix* p_input_OpMatrix_B, BasicVector* p_input_Vector, BasicMatrix* p_input_A_deflated, BasicMatrix* p_input_B_deflated, BasicMatrix* p_input_Q_Total, BasicMatrix* p_input_Z_Total, BasicMatrix* p_input_Q_Step, BasicMatrix* p_input_Z_Step, BasicMatrix* p_input_QZ_Step, BasicMatrix* p_input_TempMatrix_Trans,BasicMatrix* p_input_TempMatrix)
{
	this->init(p_input_OpMatrix_A,
			p_input_OpMatrix_B,
			p_input_Vector,
			p_input_A_deflated,
			p_input_B_deflated,
			p_input_Q_Total,
			p_input_Z_Total,
			p_input_Q_Step,
			p_input_Z_Step,
			p_input_QZ_Step,
			p_input_TempMatrix_Trans,
			p_input_TempMatrix);
};

/*
 * ����ԭʼHessenberg-Triangle��������
 */
void GeneralizedEigenSolverForReal::generateHTOpMatrix()
{
	this->m_HessenbergTriangleFormular.reload(this->p_OpMatrix_A,this->p_OpMatrix_B,p_QMatrix_Total,p_ZMatrix_Total,p_QZMatrix_Step,p_TempMatrix_Trans,p_TempMatrix);
	this->m_HessenbergTriangleFormular.formularABMatrix();
	//this->p_OpMatrix_Hessenberg->copyMatrixElementNoCheck(p_HessenbergTriangleFormular->getHessenbergMatrixA());
	//this->p_OpMatrix_Triangle->copyMatrixElementNoCheck(p_HessenbergTriangleFormular->getTriangleMatrixB());

	//QZ�任����-������
	//this->p_QMatrix_Iteration->copyMatrixElementNoCheck(p_HessenbergTriangleFormular->getMatrixQ_Total());
	//this->p_ZMatrix_Iteration->copyMatrixElementNoCheck(p_HessenbergTriangleFormular->getMatrixZ_Total());

	//��������ת������
	//updateQZMatrixTotal();

	cout << "--------generateHTOpMatrix full OpHessenbergMatrix----------" << endl;
	p_OpMatrix_A->printMatrix();
	cout << "--------generateHTOpMatrix full OpTriangleMatrix----------" << endl;
	p_OpMatrix_B->printMatrix();
};


/*
 * ���ݽ��׵� ���ɽ��׵�Hessenberg-Tirangle���� �����¹滮��Ӧ��ת������
 */
void GeneralizedEigenSolverForReal::generateDeflatedHTMatrixPair()
{
	int newSize = this->deflationEnd - this->deflationStart + 1;
	this->p_OpMatrix_Hessenberg_deflated->resizeMatrix(newSize,newSize);
	this->p_OpMatrix_Hessenberg_deflated->resetMatrixToI();
	this->p_OpMatrix_Triangle_deflated->resizeMatrix(newSize,newSize);
	this->p_OpMatrix_Triangle_deflated->resetMatrixToI();

	p_QMatrix_Step->resetMatrixToI();
	p_QMatrix_Step->resizeMatrix(newSize,newSize);

	p_ZMatrix_Step->resetMatrixToI();
	p_ZMatrix_Step->resizeMatrix(newSize,newSize);

	p_QZMatrix_Step->resetMatrixToI();
	p_QZMatrix_Step->resizeMatrix(newSize,newSize);
	/*
	p_QMatrix_Deflated_Iteration->resizeMatrix(newSize,newSize);
	p_QMatrix_Deflated_Iteration->resetMatrixToI();

	p_ZMatrix_Deflated_Iteration->resizeMatrix(newSize,newSize);
	p_ZMatrix_Deflated_Iteration->resetMatrixToI();
	*/

	double temp;
	//��ԭʼȫά��hessenberg-Triangle������ �������ѽ��׵�hessenberg-Triangle����
	for(int i=this->deflationStart,m=0; i<=this->deflationEnd; i++,m++)
	{
		for(int j=this->deflationStart,n=0; j<=this->deflationEnd; j++,n++)
		{
			temp = this->p_OpMatrix_A->getMatrixElement(i,j);
			this->p_OpMatrix_Hessenberg_deflated->setMatrixElement(m,n,temp);

			temp = this->p_OpMatrix_B->getMatrixElement(i,j);
			this->p_OpMatrix_Triangle_deflated->setMatrixElement(m,n,temp);
		}
	}

	cout << "--------generateDeflatedHTMatrixPair deflated OpHessenbergMatrix----------" << endl;
	p_OpMatrix_Hessenberg_deflated->printMatrix();
	cout << "--------generateDeflatedHTMatrixPair deflated OpTriangleMatrix----------" << endl;
	p_OpMatrix_Triangle_deflated->printMatrix();

};

/*
 * ���ѽ��׵ı任���� ������Ϊȫ�ߴ�任����
 */
void GeneralizedEigenSolverForReal::upgradeDeflatedQMatrix()
{
	int tailIndex = deflationEnd - deflationStart;
	//this->p_QMatrix_Iteration->resetMatrixToI();
	//this->p_ZMatrix_Iteration->resetMatrixToI();

	cout << "--------upgradeDeflatedQMatrix-------------" << endl;
	cout << "--------Deflation Start:" << this->deflationStart << "-----Deflation End:" << this->deflationEnd <<endl;
	cout << "--------Before Upgrade---Q_Deflated----------" << endl;
	p_QMatrix_Step->printMatrix();
	//this->p_QMatrix_Deflated_Iteration->printMatrix();
	//cout << "--------Before Upgrade---Z_Deflated----------" << endl;
	//this->p_ZMatrix_Deflated_Iteration->printMatrix();
	p_QMatrix_Step->resizeMatrix(this->p_OpMatrix_A->rowNum,this->p_OpMatrix_A->columnNum);
	p_QMatrix_Step->moveDiagonalSubMatrixDown(0,tailIndex,deflationStart);

	/*
	for(int i=this->deflationStart,m=0; i<=this->deflationEnd; i++,m++)
	{
		for(int j=this->deflationStart,n=0; j<=this->deflationEnd; j++,n++)
		{
			//����QT
			temp = this->p_QMatrix_Deflated_Iteration->getMatrixElement(m,n);
			this->p_QMatrix_Iteration->setMatrixElement(i,j,temp);

			//����Q
			temp = this->p_ZMatrix_Deflated_Iteration->getMatrixElement(m,n);
			this->p_ZMatrix_Iteration->setMatrixElement(i,j,temp);
		}
	}
	*/
	cout << "--------After Upgrade---Q_FullSize----------" << endl;
	p_QMatrix_Step->printMatrix();
	//this->p_QMatrix_Iteration->printMatrix();
	//cout << "--------After Upgrade---Z_FullSize----------" << endl;
	//this->p_ZMatrix_Iteration->printMatrix();
};

/*
 * ���ѽ��׵ı任���� ������Ϊȫ�ߴ�任����
 */
void GeneralizedEigenSolverForReal::upgradeDeflatedZMatrix()
{
	int tailIndex = deflationEnd - deflationStart;
	//this->p_QMatrix_Iteration->resetMatrixToI();
	//this->p_ZMatrix_Iteration->resetMatrixToI();

	cout << "--------upgradeDeflatedQZMatrix-------------" << endl;
	cout << "--------Deflation Start:" << this->deflationStart << "-----Deflation End:" << this->deflationEnd <<endl;
	//cout << "--------Before Upgrade---Q_Deflated----------" << endl;
	//this->p_QMatrix_Deflated_Iteration->printMatrix();
	cout << "--------Before Upgrade---Z_Deflated----------" << endl;
	//this->p_ZMatrix_Deflated_Iteration->printMatrix();
	p_ZMatrix_Step->printMatrix();

	p_ZMatrix_Step->resizeMatrix(this->p_OpMatrix_A->rowNum,this->p_OpMatrix_A->columnNum);
	p_ZMatrix_Step->moveDiagonalSubMatrixDown(0,tailIndex,deflationStart);

	/*
	for(int i=this->deflationStart,m=0; i<=this->deflationEnd; i++,m++)
	{
		for(int j=this->deflationStart,n=0; j<=this->deflationEnd; j++,n++)
		{
			//����QT
			temp = this->p_QMatrix_Deflated_Iteration->getMatrixElement(m,n);
			this->p_QMatrix_Iteration->setMatrixElement(i,j,temp);

			//����Q
			temp = this->p_ZMatrix_Deflated_Iteration->getMatrixElement(m,n);
			this->p_ZMatrix_Iteration->setMatrixElement(i,j,temp);
		}
	}
	*/
	//cout << "--------After Upgrade---Q_FullSize----------" << endl;
	//this->p_QMatrix_Iteration->printMatrix();
	cout << "--------After Upgrade---Z_FullSize----------" << endl;
	//this->p_ZMatrix_Iteration->printMatrix();
	p_ZMatrix_Step->printMatrix();
};

/*
 * ���������任���� �ϲ���Ϊ����任����
 */
void GeneralizedEigenSolverForReal::updateQMatrixTotal()
{
	//����Q_It * Q_Total
	this->m_Multiplier.reload(p_QMatrix_Step, p_QMatrix_Total, p_TempMatrix);
	this->m_Multiplier.multiplyCalc();
	p_QMatrix_Total->copyMatrixElementNoCheck(p_TempMatrix);

	cout << "--------updateQMatrixTotal----------" << endl;
	cout << "--------After Upgrade---Q_Total FullSize----------" << endl;
	this->p_QMatrix_Total->printMatrix();
	//cout << "--------After Upgrade---Z_Total FullSize----------" << endl;
	//this->p_ZMatrix_Total->printMatrix();
};

/*
 * ���������任���� �ϲ���Ϊ����任����
 */
void GeneralizedEigenSolverForReal::updateZMatrixTotal()
{
	//����Z_Total * Z_It
	this->m_Multiplier.reload(p_ZMatrix_Total, p_ZMatrix_Step, p_TempMatrix);
	this->m_Multiplier.multiplyCalc();
	p_ZMatrix_Total->copyMatrixElementNoCheck(p_TempMatrix);
	cout << "--------updateZMatrixTotal----------" << endl;
	//cout << "--------After Upgrade---Q_Total FullSize----------" << endl;
	//this->p_QMatrix_Total->printMatrix();
	cout << "--------After Upgrade---Z_Total FullSize----------" << endl;
	this->p_ZMatrix_Total->printMatrix();
};

//ʹ��Q�������H-T�����
void GeneralizedEigenSolverForReal::updateHTMatrixByQ()
{
	cout << "GeneralizedEigenSolverForReal--updateHTMatrixByQ----OP Hessenberg Matrix Before" << endl;
	this->p_OpMatrix_A->printMatrix();
	cout << "GeneralizedEigenSolverForReal--updateHTMatrixByQ----OP Triangle Matrix Before" << endl;
	this->p_OpMatrix_B->printMatrix();

	//���¾���Hessenberg A
	this->m_Multiplier.reload(p_QMatrix_Step, p_OpMatrix_A, p_TempMatrix);
	this->m_Multiplier.multiplyCalc();
	p_OpMatrix_A->copyMatrixElementNoCheck(p_TempMatrix);

	//���¾���Triangle B
	this->m_Multiplier.reload(p_QMatrix_Step, p_OpMatrix_B, p_TempMatrix);
	this->m_Multiplier.multiplyCalc();
	p_OpMatrix_B->copyMatrixElementNoCheck(p_TempMatrix);

	cout << "GeneralizedEigenSolverForReal--updateHTMatrixByQ----OP Hessenberg Matrix After" << endl;
	this->p_OpMatrix_A->printMatrix();
	cout << "GeneralizedEigenSolverForReal--updateHTMatrixByQ----OP Triangle Matrix After" << endl;
	this->p_OpMatrix_B->printMatrix();
};

//ʹ��Z�����ҳ�H-T�����
void GeneralizedEigenSolverForReal::updateHTMatrixByZ()
{
	cout << "GeneralizedEigenSolverForReal--updateHTMatrixByZ----OP Hessenberg Matrix Before" << endl;
	this->p_OpMatrix_A->printMatrix();
	cout << "GeneralizedEigenSolverForReal--updateHTMatrixByZ----OP Triangle Matrix Before" << endl;
	this->p_OpMatrix_B->printMatrix();

	//���¾���Hessenberg A
	this->m_Multiplier.reload(p_OpMatrix_A, p_ZMatrix_Step, p_TempMatrix);
	this->m_Multiplier.multiplyCalc();
	p_OpMatrix_A->copyMatrixElementNoCheck(p_TempMatrix);

	//���¾���Triangle B
	this->m_Multiplier.reload(p_OpMatrix_B, p_ZMatrix_Step,p_TempMatrix);
	this->m_Multiplier.multiplyCalc();
	p_OpMatrix_B->copyMatrixElementNoCheck(p_TempMatrix);

	cout << "GeneralizedEigenSolverForReal--updateHTMatrixByZ----OP Hessenberg Matrix After" << endl;
	this->p_OpMatrix_A->printMatrix();
	cout << "GeneralizedEigenSolverForReal--updateHTMatrixByZ----OP Triangle Matrix After" << endl;
	this->p_OpMatrix_B->printMatrix();
};


//��ʼ������ֵ������ؾ���
void GeneralizedEigenSolverForReal::initEigenCalcMatrix()
{
	this->p_OpMatrix_Hessenberg_deflated->resizeMatrix(this->p_OpMatrix_A->rowNum,this->p_OpMatrix_A->columnNum);
	p_OpMatrix_Hessenberg_deflated->copyMatrixElementNoCheck(p_OpMatrix_A);
	this->p_OpMatrix_Triangle_deflated->resizeMatrix(this->p_OpMatrix_B->rowNum,this->p_OpMatrix_B->columnNum);
	p_OpMatrix_Triangle_deflated->copyMatrixElementNoCheck(p_OpMatrix_B);

	p_QMatrix_Step->resetMatrixToI();
	p_QMatrix_Step->resizeMatrix(this->p_OpMatrix_A->rowNum,this->p_OpMatrix_A->columnNum);

	p_ZMatrix_Step->resetMatrixToI();
	p_ZMatrix_Step->resizeMatrix(this->p_OpMatrix_A->rowNum,this->p_OpMatrix_A->columnNum);

	p_QZMatrix_Step->resetMatrixToI();
	p_QZMatrix_Step->resizeMatrix(this->p_OpMatrix_A->rowNum,this->p_OpMatrix_A->columnNum);

	p_TempMatrix_Trans->resetMatrixToI();
	p_TempMatrix_Trans->resizeMatrix(this->p_OpMatrix_A->rowNum,this->p_OpMatrix_A->columnNum);

	p_TempMatrix->resetMatrixToI();
	p_TempMatrix->resizeMatrix(this->p_OpMatrix_A->rowNum,this->p_OpMatrix_A->columnNum);

	m_DoubleShifeQZ.reload(p_OpMatrix_Hessenberg_deflated, p_OpMatrix_Triangle_deflated, p_OpTransVector, p_QMatrix_Step, p_ZMatrix_Step,p_QZMatrix_Step, p_TempMatrix_Trans, p_TempMatrix);

	//p_QZMatrix_Deflated_Iteration->resizeMatrix(this->p_OpMatrix_Hessenberg->rowNum,this->p_OpMatrix_Hessenberg->columnNum);
	//p_QZMatrix_Deflated_Iteration->resetMatrixToI();

	//p_ZMatrix_Deflated_Iteration->resizeMatrix(this->p_OpMatrix_Hessenberg->rowNum,this->p_OpMatrix_Hessenberg->columnNum);
	//p_ZMatrix_Deflated_Iteration->resetMatrixToI();

	//ԭʼ����ĵ������任����
	//p_QMatrix_Iteration->resetMatrixToI();
	//p_ZMatrix_Iteration->resetMatrixToI();


};
//����resize��С������ؾ���
void GeneralizedEigenSolverForReal::resizeTransMatrix()
{
	int size = this->deflationEnd - this->deflationStart + 1;

	p_QMatrix_Step->resetMatrixToI();
	p_QMatrix_Step->resizeMatrix(size,size);

	p_ZMatrix_Step->resetMatrixToI();
	p_ZMatrix_Step->resizeMatrix(size,size);

	p_QZMatrix_Step->resetMatrixToI();
	p_QZMatrix_Step->resizeMatrix(size,size);

	p_TempMatrix_Trans->resetMatrixToI();
	p_TempMatrix_Trans->resizeMatrix(size,size);

	p_TempMatrix->resetMatrixToI();
	p_TempMatrix->resizeMatrix(size,size);
};
//�����ת����������Ϊȫά��
void GeneralizedEigenSolverForReal::upgradTransMatrix()
{
	p_TempMatrix_Trans->resizeMatrix(this->p_OpMatrix_A->rowNum,this->p_OpMatrix_A->columnNum);
	p_TempMatrix->resizeMatrix(this->p_OpMatrix_A->rowNum,this->p_OpMatrix_A->columnNum);
};

//��������ж��߼�
bool GeneralizedEigenSolverForReal::hasFinishedIteration()
{
	int distanceSize = this->deflationEnd - this->deflationStart + 1;
	if(distanceSize<=2)
	{
		return true;
	}
	return false;
};

//����2x2�Խǿ��Ƿ�Ϊ��������ֵ�ж�
bool GeneralizedEigenSolverForReal::isDiagonalBlockComplexEigen(BasicMatrix* p_Input_OpMatrix_A, BasicMatrix* p_Input_OpMatrix_B)
{
	this->m_ABInvCalc.generateABinvFirst2x2(p_Input_OpMatrix_A,p_Input_OpMatrix_B);
	double temp;
	double a,b,c,d;

	a = m_ABInvCalc.getABinv11();
	b = m_ABInvCalc.getABinv12();
	c = m_ABInvCalc.getABinv21();
	d = m_ABInvCalc.getABinv22();

	if(0 == c)
	{
		return false;
	}

	temp = 0.5*a + 0.5*d;
	temp = temp*temp;
	temp = temp-a*d + b*c;

	if(temp<0)
	{
		//�Ը��������ţ�����ֵΪ����
		return true;
	}

	return false;
};

//���Դ�ӡ��Q_Total * OP * Z_Total

void GeneralizedEigenSolverForReal::showQxOPxZ()
{
	//print QAZ
	m_testMulti.reload(this->p_QMatrix_Total, &testForTemp_A_nxn,&testTemp_nxn);
	m_testMulti.multiplyCalc();
	testForTemp_A_nxn.copyMatrixElementNoCheck(&testTemp_nxn);
	m_testMulti.reload(&testForTemp_A_nxn,this->p_ZMatrix_Total,&testTemp_nxn);
	m_testMulti.multiplyCalc();
	testForTemp_A_nxn.copyMatrixElementNoCheck(&testTemp_nxn);
	cout << "--------------Q_Total * Op_Matrix_A * Z_Total-----------------" << endl;
	testForTemp_A_nxn.printMatrix();

	//print QBZ
	m_testMulti.reload(this->p_QMatrix_Total, &testForTemp_B_nxn,&testTemp_nxn);
	m_testMulti.multiplyCalc();
	testForTemp_B_nxn.copyMatrixElementNoCheck(&testTemp_nxn);
	m_testMulti.reload(&testForTemp_B_nxn,this->p_ZMatrix_Total,&testTemp_nxn);
	m_testMulti.multiplyCalc();
	testForTemp_B_nxn.copyMatrixElementNoCheck(&testTemp_nxn);
	cout << "--------------Q_Total * Op_Matrix_B * Z_Total-----------------" << endl;
	testForTemp_B_nxn.printMatrix();
};

//��������ֵ
void GeneralizedEigenSolverForReal::calcEigenValue()
{
	this->deflationStart = 0;
	this->deflationEnd = this->p_OpMatrix_A->rowNum - 1;
	//ԭʼ���������任����
	p_QMatrix_Total->resetMatrixToI();
	p_ZMatrix_Total->resetMatrixToI();

	//����ԭʼȫά��Hessenberg-Triangle����
	generateHTOpMatrix();

	//��ӡ����ת������˻�����
	//showQxOPxZ();

	//������Triangle�Խ�����Ԫ�½�����
	m_QZTriangleZeroChasing.reload(p_OpMatrix_A,p_OpMatrix_B,p_OpMatrix_Hessenberg_deflated,p_OpMatrix_Triangle_deflated,p_QMatrix_Step, p_ZMatrix_Step,p_QZMatrix_Step, p_TempMatrix);
	m_QZTriangleZeroChasing.deflate();
	//p_QMatrix_Iteration->copyMatrixElementNoCheck(p_QZTriangleZeroChasing->getGivensMatrix_Q_Total());
	//p_ZMatrix_Iteration->copyMatrixElementNoCheck(p_QZTriangleZeroChasing->getGivensMatrix_Z_Total());

	//����ԭʼHessenberg-Trianlge ����ʹ��Q�������H-T�����
	//updateHTMatrixByQ();
	//��������ת������
	updateQMatrixTotal();

	//����ԭʼHessenberg-Trianlge ����ʹ��Z�����ҳ�H-T�����
	//updateHTMatrixByZ();
	//��������ת������
	updateZMatrixTotal();
	//��ӡ����ת������˻�����
	//showQxOPxZ();

	//��ʼ����ؾ���
	initEigenCalcMatrix();

	//��hessenberg����������ƶԽǻ�
	while(true)
	{
		//Hessenberg����������Σ�����Сֵ��Ϊ0Ԫ
		this->p_OpMatrix_A->regularZeroElement();
		this->p_OpMatrix_B->regularZeroElement();
		cout << "--------calcEigenValue After Regular---Fullsize OpHessenbergMatrix----------" << endl;
		p_OpMatrix_A->printMatrix();
		cout << "--------calcEigenValue After Regular---Fullsize OpTriangleMatrix----------" << endl;
		p_OpMatrix_B->printMatrix();

		bool hasNewDeflate = m_HessenbergDeflation.findNewDeflationPoint(p_OpMatrix_A, this->deflationStart,this->deflationEnd);
		this->deflationStart = m_HessenbergDeflation.getNewDeflationStart();
		this->deflationEnd = m_HessenbergDeflation.getNewDeflationEnd();

		if(this->hasFinishedIteration())
		{
			break;
		}

		if(hasNewDeflate)
		{
			//���ɽ���Hessenberg-Triangle����
			generateDeflatedHTMatrixPair();
			m_DoubleShifeQZ.reload(this->p_OpMatrix_Hessenberg_deflated,this->p_OpMatrix_Triangle_deflated,p_OpTransVector, p_QMatrix_Step, p_ZMatrix_Step,p_QZMatrix_Step, p_TempMatrix_Trans, p_TempMatrix);
		}
		//����resize���ת������
		resizeTransMatrix();
		m_DoubleShifeQZ.wilkinson_IM_QZIteration_Single();
		//p_QMatrix_Deflated_Iteration->copyMatrixElementNoCheck(p_DoubleShifeQZ->getQMatrix_Total());
		//p_ZMatrix_Deflated_Iteration->copyMatrixElementNoCheck(p_DoubleShifeQZ->getZMatrix_Total());

		//cout << "--------calcEigenValue deflated matrix After wilkinson--- Q_DF * OpHessenbergMatrix_Df * Z_Df----------" << endl;
		//this->p_Multiplier->reload(p_QMatrix_Deflated_Iteration,p_OpMatrix_Hessenberg_deflated);
		//this->p_Multiplier->multiplyCalc();
		//p_OpMatrix_Hessenberg_deflated->copyMatrixElementNoCheck(this->p_Multiplier->getMultiplyResult());
		//this->p_Multiplier->reload(p_OpMatrix_Hessenberg_deflated,p_ZMatrix_Deflated_Iteration);
		//this->p_Multiplier->multiplyCalc();
		//p_OpMatrix_Hessenberg_deflated->copyMatrixElementNoCheck(this->p_Multiplier->getMultiplyResult());
		//p_OpMatrix_Hessenberg_deflated->printMatrix();

		//cout << "--------calcEigenValue deflated matrix After wilkinson--- Q_DF * OpTriangleMatrix_Df * Z_Df----------" << endl;
		//this->p_Multiplier->reload(p_QMatrix_Deflated_Iteration,p_OpMatrix_Triangle_deflated);
		//this->p_Multiplier->multiplyCalc();
		//p_OpMatrix_Triangle_deflated->copyMatrixElementNoCheck(this->p_Multiplier->getMultiplyResult());
		//this->p_Multiplier->reload(p_OpMatrix_Triangle_deflated,p_ZMatrix_Deflated_Iteration);
		//this->p_Multiplier->multiplyCalc();
		//p_OpMatrix_Triangle_deflated->copyMatrixElementNoCheck(this->p_Multiplier->getMultiplyResult());
		//p_OpMatrix_Triangle_deflated->printMatrix();

		//������ת����������Ϊȫά��ת������
		upgradeDeflatedQMatrix();

		//�����ת����������Ϊȫά��
		upgradTransMatrix();

		//����ԭʼHessenberg-Trianlge ����ʹ��Q�������H-T�����
		updateHTMatrixByQ();
		//��������ת������
		updateQMatrixTotal();

		//������ת����������Ϊȫά��ת������
		upgradeDeflatedZMatrix();
		//����ԭʼHessenberg-Trianlge ����ʹ��Z�����ҳ�H-T�����
		updateHTMatrixByZ();
		//��������ת������
		updateZMatrixTotal();
		//��ӡ����ת������˻�����
		//showQxOPxZ();
	}
	//��ӡ����ת������˻�����
	//showQxOPxZ();
	int i = 0;
	double subDiagonalElement;
	while (i < this->p_OpMatrix_A->columnNum - 1)
	{
		subDiagonalElement = p_OpMatrix_A->getMatrixElement(i+1,i);
		if(0 == subDiagonalElement)
		{
			i++;
			continue;
		}
		p_OpMatrix_Hessenberg_deflated->resizeMatrix(2,2);
		this->p_OpMatrix_Hessenberg_deflated->setMatrixElement(0,0,p_OpMatrix_A->getMatrixElement(i,i));
		this->p_OpMatrix_Hessenberg_deflated->setMatrixElement(0,1,p_OpMatrix_A->getMatrixElement(i,i+1));
		this->p_OpMatrix_Hessenberg_deflated->setMatrixElement(1,0,p_OpMatrix_A->getMatrixElement(i+1,i));
		this->p_OpMatrix_Hessenberg_deflated->setMatrixElement(1,1,p_OpMatrix_A->getMatrixElement(i+1,i+1));

		p_OpMatrix_Triangle_deflated->resizeMatrix(2,2);
		this->p_OpMatrix_Triangle_deflated->setMatrixElement(0,0,p_OpMatrix_B->getMatrixElement(i,i));
		this->p_OpMatrix_Triangle_deflated->setMatrixElement(0,1,p_OpMatrix_B->getMatrixElement(i,i+1));
		this->p_OpMatrix_Triangle_deflated->setMatrixElement(1,0,p_OpMatrix_B->getMatrixElement(i+1,i));
		this->p_OpMatrix_Triangle_deflated->setMatrixElement(1,1,p_OpMatrix_B->getMatrixElement(i+1,i+1));

		//�ж϶Խ��ߵ�2x2������Ƿ���Ҫ�����Խǻ�
		if(!isDiagonalBlockComplexEigen(p_OpMatrix_Hessenberg_deflated,p_OpMatrix_Triangle_deflated))
		{
			lastStepIteration(i);
		}
		i++;
	}
	//��ӡ����ת������˻�����
	//showQxOPxZ();
};

//������Ϊ�Խǿ��Ժ����һ����������0,1�����ϵ�2x2�Խǿ���������ǻ�
void GeneralizedEigenSolverForReal::lastStepIteration(int startIndex)
{
	double temp;
	//�����趨deflate��ֹ��
	this->deflationStart = startIndex;
	this->deflationEnd = startIndex + 1;
	//���ɽ��׵�hessenberg-Triangle����
	generateDeflatedHTMatrixPair();
	//����resize���ת������
	resizeTransMatrix();

	m_SingleShifeQZ.reload(p_OpMatrix_Hessenberg_deflated,p_OpMatrix_Triangle_deflated,p_QMatrix_Step,p_ZMatrix_Step, p_QZMatrix_Step, p_TempMatrix_Trans,p_TempMatrix);

	while(true)
	{
		//����resize���ת������
		resizeTransMatrix();

		m_SingleShifeQZ.rayleigh_Quotient_IM_QZIteration_Step();
		//p_QMatrix_Deflated_Iteration->copyMatrixElementNoCheck(p_SingleShifeQZ->getQMatrix_Total());
		//p_ZMatrix_Deflated_Iteration->copyMatrixElementNoCheck(p_SingleShifeQZ->getZMatrix_Total());

		//������ת����������Ϊȫά��ת������
		upgradeDeflatedQMatrix();
		//�����ת����������Ϊȫά��
		upgradTransMatrix();

		//����ԭʼHessenberg-Trianlge ����ʹ��Q�������H-T�����
		updateHTMatrixByQ();
		//��������ת������
		updateQMatrixTotal();


		//������ת����������Ϊȫά��ת������
		upgradeDeflatedZMatrix();
		//����ԭʼHessenberg-Trianlge ����ʹ��Z�����ҳ�H-T�����
		updateHTMatrixByZ();
		//��������ת������
		updateZMatrixTotal();

		//��ӡ����ת������˻�����
		//showQxOPxZ();

		//Hessenberg����������Σ�����Сֵ��Ϊ0Ԫ
		p_OpMatrix_A->regularZeroElement();
		cout << "--------calcEigenValue lastStepIteration After Regular---Fullsize OpHessenbergMatrix----------" << endl;
		p_OpMatrix_A->printMatrix();

		p_OpMatrix_B->regularZeroElement();
		cout << "--------calcEigenValue lastStepIteration After Regular---Fullsize OpTriangleMatrix----------" << endl;
		p_OpMatrix_B->printMatrix();

		temp = this->p_OpMatrix_A->getMatrixElement(startIndex+1,startIndex);
		if(0 == temp)
		{
			break;
		}
	}
};

//��ȡH-T�����
BasicMatrix* GeneralizedEigenSolverForReal::getHessenbergMatrix()
{
	return this->p_OpMatrix_A;
};
BasicMatrix* GeneralizedEigenSolverForReal::getTriangleMatrix()
{
	return this->p_OpMatrix_B;
};

//��ȡQ\Z �ۺ�ת������
BasicMatrix* GeneralizedEigenSolverForReal::getQMatrix_Total()
{
	return this->p_QMatrix_Total;
};
BasicMatrix* GeneralizedEigenSolverForReal::getZMatrix_Total()
{
	return this->p_ZMatrix_Total;
};
