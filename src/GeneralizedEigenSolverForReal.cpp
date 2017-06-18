/*
 * GeneralizedEigenSolverForReal.cpp
 *
 *  Created on: 2017年6月1日
 *      Author: looke
 */

#include "GeneralizedEigenSolverForReal.h"
#include "iostream"
using namespace std;

//GeneralizedEigenSolverForReal::GeneralizedEigenSolverForReal()
//{};

GeneralizedEigenSolverForReal::GeneralizedEigenSolverForReal(BasicMatrix* p_input_OpMatrix_A, BasicMatrix* p_input_OpMatrix_B, BasicVector* p_input_Vector, BasicMatrix* p_input_A_deflated, BasicMatrix* p_input_B_deflated, BasicMatrix* p_input_Q_Total, BasicMatrix* p_input_Z_Total, BasicMatrix* p_input_Q_Step, BasicMatrix* p_input_Z_Step, BasicMatrix* p_input_QZ_Step, BasicMatrix* p_input_TempMatrix_Trans,BasicMatrix* p_input_TempMatrix)
:m_ABInvCalc(),m_HessenbergDeflation(),m_Multiplier(p_input_OpMatrix_A,p_input_OpMatrix_A,p_input_TempMatrix),
 m_SingleShifeQZ(p_input_OpMatrix_A,p_input_OpMatrix_B,p_input_Q_Step,p_input_Z_Step,p_input_QZ_Step,p_input_TempMatrix_Trans,p_input_TempMatrix),
 m_DoubleShifeQZ(p_input_OpMatrix_A,p_input_OpMatrix_B,p_input_Vector,p_input_Q_Step,p_input_Z_Step,p_input_QZ_Step,p_input_TempMatrix_Trans,p_input_TempMatrix),
 m_HessenbergTriangleFormular(p_input_OpMatrix_A,p_input_OpMatrix_B,p_input_Q_Total,p_input_Z_Total,p_input_QZ_Step,p_input_TempMatrix_Trans,p_input_TempMatrix),
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
	//原始操作矩阵
	this->p_OpMatrix_A = p_input_OpMatrix_A;
	this->p_OpMatrix_B = p_input_OpMatrix_B;

	this->p_OpTransVector = p_input_Vector;

	//原始矩阵的全维度总体变换矩阵
	this->p_QMatrix_Total = p_input_Q_Total;
	this->p_ZMatrix_Total = p_input_Z_Total;

	//原始矩阵的 迭代变换矩阵 Q\Z
	this->p_QMatrix_Step = p_input_Q_Step;
	this->p_ZMatrix_Step = p_input_Z_Step;

	//已降阶的 操作Hessenberg矩阵
	this->p_OpMatrix_Hessenberg_deflated = p_input_A_deflated;
	//已降阶的 操作Triangle矩阵
	this->p_OpMatrix_Triangle_deflated = p_input_B_deflated;

	//中间过程矩阵
	this->p_QZMatrix_Step = p_input_QZ_Step;
	//中间过程矩阵
	this->p_TempMatrix_Trans = p_input_TempMatrix_Trans;
	//中间过程矩阵
	this->p_TempMatrix = p_input_TempMatrix;

	//降阶起点索引指示，初始化为0
	this->deflationStart = 0;
	//降阶终点索引指示，初始化为n-1
	this->deflationEnd = p_OpMatrix_A->rowNum - 1;

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
 * 生成原始Hessenberg-Triangle操作矩阵
 */
void GeneralizedEigenSolverForReal::generateHTOpMatrix()
{
	this->m_HessenbergTriangleFormular.reload(this->p_OpMatrix_A,this->p_OpMatrix_B,p_QMatrix_Total,p_ZMatrix_Total,p_QZMatrix_Step,p_TempMatrix_Trans,p_TempMatrix);
	this->m_HessenbergTriangleFormular.formularABMatrix();
	//this->p_OpMatrix_Hessenberg->copyMatrixElementNoCheck(p_HessenbergTriangleFormular->getHessenbergMatrixA());
	//this->p_OpMatrix_Triangle->copyMatrixElementNoCheck(p_HessenbergTriangleFormular->getTriangleMatrixB());

	//QZ变换矩阵-单迭代
	//this->p_QMatrix_Iteration->copyMatrixElementNoCheck(p_HessenbergTriangleFormular->getMatrixQ_Total());
	//this->p_ZMatrix_Iteration->copyMatrixElementNoCheck(p_HessenbergTriangleFormular->getMatrixZ_Total());

	//更新总体转换矩阵
	//updateQZMatrixTotal();

	cout << "--------generateHTOpMatrix full OpHessenbergMatrix----------" << endl;
	p_OpMatrix_A->printMatrix();
	cout << "--------generateHTOpMatrix full OpTriangleMatrix----------" << endl;
	p_OpMatrix_B->printMatrix();
};


/*
 * 根据降阶点 生成降阶的Hessenberg-Tirangle矩阵 并重新规划相应的转换矩阵
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
	//从原始全维度hessenberg-Triangle矩阵中 拷贝出已降阶的hessenberg-Triangle矩阵
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
 * 将已降阶的变换矩阵 升级成为全尺寸变换矩阵
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
			//更新QT
			temp = this->p_QMatrix_Deflated_Iteration->getMatrixElement(m,n);
			this->p_QMatrix_Iteration->setMatrixElement(i,j,temp);

			//更新Q
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
 * 将已降阶的变换矩阵 升级成为全尺寸变换矩阵
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
			//更新QT
			temp = this->p_QMatrix_Deflated_Iteration->getMatrixElement(m,n);
			this->p_QMatrix_Iteration->setMatrixElement(i,j,temp);

			//更新Q
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
 * 将单迭代变换矩阵 合并成为总体变换矩阵
 */
void GeneralizedEigenSolverForReal::updateQMatrixTotal()
{
	//计算Q_It * Q_Total
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
 * 将单迭代变换矩阵 合并成为总体变换矩阵
 */
void GeneralizedEigenSolverForReal::updateZMatrixTotal()
{
	//计算Z_Total * Z_It
	this->m_Multiplier.reload(p_ZMatrix_Total, p_ZMatrix_Step, p_TempMatrix);
	this->m_Multiplier.multiplyCalc();
	p_ZMatrix_Total->copyMatrixElementNoCheck(p_TempMatrix);
	cout << "--------updateZMatrixTotal----------" << endl;
	//cout << "--------After Upgrade---Q_Total FullSize----------" << endl;
	//this->p_QMatrix_Total->printMatrix();
	cout << "--------After Upgrade---Z_Total FullSize----------" << endl;
	this->p_ZMatrix_Total->printMatrix();
};

//使用Q矩阵左乘H-T矩阵对
void GeneralizedEigenSolverForReal::updateHTMatrixByQ()
{
	cout << "GeneralizedEigenSolverForReal--updateHTMatrixByQ----OP Hessenberg Matrix Before" << endl;
	this->p_OpMatrix_A->printMatrix();
	cout << "GeneralizedEigenSolverForReal--updateHTMatrixByQ----OP Triangle Matrix Before" << endl;
	this->p_OpMatrix_B->printMatrix();

	//更新矩阵Hessenberg A
	this->m_Multiplier.reload(p_QMatrix_Step, p_OpMatrix_A, p_TempMatrix);
	this->m_Multiplier.multiplyCalc();
	p_OpMatrix_A->copyMatrixElementNoCheck(p_TempMatrix);

	//更新矩阵Triangle B
	this->m_Multiplier.reload(p_QMatrix_Step, p_OpMatrix_B, p_TempMatrix);
	this->m_Multiplier.multiplyCalc();
	p_OpMatrix_B->copyMatrixElementNoCheck(p_TempMatrix);

	cout << "GeneralizedEigenSolverForReal--updateHTMatrixByQ----OP Hessenberg Matrix After" << endl;
	this->p_OpMatrix_A->printMatrix();
	cout << "GeneralizedEigenSolverForReal--updateHTMatrixByQ----OP Triangle Matrix After" << endl;
	this->p_OpMatrix_B->printMatrix();
};

//使用Z矩阵右乘H-T矩阵对
void GeneralizedEigenSolverForReal::updateHTMatrixByZ()
{
	cout << "GeneralizedEigenSolverForReal--updateHTMatrixByZ----OP Hessenberg Matrix Before" << endl;
	this->p_OpMatrix_A->printMatrix();
	cout << "GeneralizedEigenSolverForReal--updateHTMatrixByZ----OP Triangle Matrix Before" << endl;
	this->p_OpMatrix_B->printMatrix();

	//更新矩阵Hessenberg A
	this->m_Multiplier.reload(p_OpMatrix_A, p_ZMatrix_Step, p_TempMatrix);
	this->m_Multiplier.multiplyCalc();
	p_OpMatrix_A->copyMatrixElementNoCheck(p_TempMatrix);

	//更新矩阵Triangle B
	this->m_Multiplier.reload(p_OpMatrix_B, p_ZMatrix_Step,p_TempMatrix);
	this->m_Multiplier.multiplyCalc();
	p_OpMatrix_B->copyMatrixElementNoCheck(p_TempMatrix);

	cout << "GeneralizedEigenSolverForReal--updateHTMatrixByZ----OP Hessenberg Matrix After" << endl;
	this->p_OpMatrix_A->printMatrix();
	cout << "GeneralizedEigenSolverForReal--updateHTMatrixByZ----OP Triangle Matrix After" << endl;
	this->p_OpMatrix_B->printMatrix();
};


//初始化特征值计算相关矩阵
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

	//原始矩阵的单迭代变换矩阵
	//p_QMatrix_Iteration->resetMatrixToI();
	//p_ZMatrix_Iteration->resetMatrixToI();


};
//重新resize缩小计算相关矩阵
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
//将相关转换矩阵升级为全维度
void GeneralizedEigenSolverForReal::upgradTransMatrix()
{
	p_TempMatrix_Trans->resizeMatrix(this->p_OpMatrix_A->rowNum,this->p_OpMatrix_A->columnNum);
	p_TempMatrix->resizeMatrix(this->p_OpMatrix_A->rowNum,this->p_OpMatrix_A->columnNum);
};

//解算结束判断逻辑
bool GeneralizedEigenSolverForReal::hasFinishedIteration()
{
	int distanceSize = this->deflationEnd - this->deflationStart + 1;
	if(distanceSize<=2)
	{
		return true;
	}
	return false;
};

//左上2x2对角块是否为复数特征值判断
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
		//对负数开根号，特征值为复数
		return true;
	}

	return false;
};

//测试打印，Q_Total * OP * Z_Total
/*
void GeneralizedEigenSolverForReal::showQxOPxZ()
{
	p_testForTemp_A_nxn->copyMatrixElementNoCheck(this->p_OpMatrix_A);
	p_testMulti->reload(this->p_QMatrix_Total, p_testForTemp_A_nxn);
	p_testMulti->multiplyCalc();
	p_testForTemp_A_nxn->copyMatrixElementNoCheck(p_testMulti->getMultiplyResult());
	p_testMulti->reload(p_testForTemp_A_nxn,this->p_ZMatrix_Total);
	p_testMulti->multiplyCalc();
	p_testForTemp_A_nxn->copyMatrixElementNoCheck(p_testMulti->getMultiplyResult());
	cout << "--------------Q_Total * Op_Matrix_A * Z_Total-----------------" << endl;
	p_testForTemp_A_nxn->printMatrix();

	p_testForTemp_B_nxn->copyMatrixElementNoCheck(this->p_OpMatrix_B);
	p_testMulti->reload(this->p_QMatrix_Total, p_testForTemp_B_nxn);
	p_testMulti->multiplyCalc();
	p_testForTemp_B_nxn->copyMatrixElementNoCheck(p_testMulti->getMultiplyResult());
	p_testMulti->reload(p_testForTemp_B_nxn,this->p_ZMatrix_Total);
	p_testMulti->multiplyCalc();
	p_testForTemp_B_nxn->copyMatrixElementNoCheck(p_testMulti->getMultiplyResult());
	cout << "--------------Q_Total * Op_Matrix_B * Z_Total-----------------" << endl;
	p_testForTemp_B_nxn->printMatrix();
};
*/
//计算特征值
void GeneralizedEigenSolverForReal::calcEigenValue()
{
	this->deflationStart = 0;
	this->deflationEnd = this->p_OpMatrix_A->rowNum - 1;
	//原始矩阵的总体变换矩阵
	p_QMatrix_Total->resetMatrixToI();
	p_ZMatrix_Total->resetMatrixToI();

	//生成原始全维度Hessenberg-Triangle矩阵
	generateHTOpMatrix();

	//打印总体转换矩阵乘积测试
	//showQxOPxZ();

	//待增加Triangle对角线零元下降过程
	m_QZTriangleZeroChasing.reload(p_OpMatrix_A,p_OpMatrix_B,p_OpMatrix_Hessenberg_deflated,p_OpMatrix_Triangle_deflated,p_QMatrix_Step, p_ZMatrix_Step,p_QZMatrix_Step, p_TempMatrix);
	m_QZTriangleZeroChasing.deflate();
	//p_QMatrix_Iteration->copyMatrixElementNoCheck(p_QZTriangleZeroChasing->getGivensMatrix_Q_Total());
	//p_ZMatrix_Iteration->copyMatrixElementNoCheck(p_QZTriangleZeroChasing->getGivensMatrix_Z_Total());

	//更新原始Hessenberg-Trianlge 矩阵使用Q矩阵左乘H-T矩阵对
	updateHTMatrixByQ();
	//更新总体转换矩阵
	updateQMatrixTotal();

	//更新原始Hessenberg-Trianlge 矩阵使用Z矩阵右乘H-T矩阵对
	updateHTMatrixByZ();
	//更新总体转换矩阵
	updateZMatrixTotal();
	//打印总体转换矩阵乘积测试
	//showQxOPxZ();

	//初始化相关矩阵
	initEigenCalcMatrix();

	//将hessenberg操作矩阵近似对角化
	while(true)
	{
		//Hessenberg矩阵进行整形，将极小值变为0元
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
			//生成降阶Hessenberg-Triangle矩阵
			generateDeflatedHTMatrixPair();
			m_DoubleShifeQZ.reload(this->p_OpMatrix_Hessenberg_deflated,this->p_OpMatrix_Triangle_deflated,p_OpTransVector, p_QMatrix_Step, p_ZMatrix_Step,p_QZMatrix_Step, p_TempMatrix_Trans, p_TempMatrix);
		}
		//重新resize相关转换矩阵
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

		//将降阶转换矩阵升级为全维度转换矩阵
		upgradeDeflatedQMatrix();

		//将相关转换矩阵升级为全维度
		upgradTransMatrix();

		//更新原始Hessenberg-Trianlge 矩阵使用Q矩阵左乘H-T矩阵对
		updateHTMatrixByQ();
		//更新总体转换矩阵
		updateQMatrixTotal();

		//将降阶转换矩阵升级为全维度转换矩阵
		upgradeDeflatedZMatrix();
		//更新原始Hessenberg-Trianlge 矩阵使用Z矩阵右乘H-T矩阵对
		updateHTMatrixByZ();
		//更新总体转换矩阵
		updateZMatrixTotal();
		//打印总体转换矩阵乘积测试
		//showQxOPxZ();
	}
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

		//判断对角线的2x2矩阵块是否需要继续对角化
		if(!isDiagonalBlockComplexEigen(p_OpMatrix_Hessenberg_deflated,p_OpMatrix_Triangle_deflated))
		{
			lastStepIteration(i);
		}
		i++;
	}
};

//初步化为对角块以后，最后一步迭代，将0,1索引上的2x2对角块进行上三角化
void GeneralizedEigenSolverForReal::lastStepIteration(int startIndex)
{
	double temp;
	//重新设定deflate起止点
	this->deflationStart = startIndex;
	this->deflationEnd = startIndex + 1;
	//生成降阶的hessenberg-Triangle矩阵
	generateDeflatedHTMatrixPair();
	//重新resize相关转换矩阵
	resizeTransMatrix();

	m_SingleShifeQZ.reload(p_OpMatrix_Hessenberg_deflated,p_OpMatrix_Triangle_deflated,p_QMatrix_Step,p_ZMatrix_Step, p_QZMatrix_Step, p_TempMatrix_Trans,p_TempMatrix);

	while(true)
	{
		//重新resize相关转换矩阵
		resizeTransMatrix();

		m_SingleShifeQZ.rayleigh_Quotient_IM_QZIteration_Step();
		//p_QMatrix_Deflated_Iteration->copyMatrixElementNoCheck(p_SingleShifeQZ->getQMatrix_Total());
		//p_ZMatrix_Deflated_Iteration->copyMatrixElementNoCheck(p_SingleShifeQZ->getZMatrix_Total());

		//将降阶转换矩阵升级为全维度转换矩阵
		upgradeDeflatedQMatrix();
		//将相关转换矩阵升级为全维度
		upgradTransMatrix();

		//更新原始Hessenberg-Trianlge 矩阵使用Q矩阵左乘H-T矩阵对
		updateHTMatrixByQ();
		//更新总体转换矩阵
		updateQMatrixTotal();


		//将降阶转换矩阵升级为全维度转换矩阵
		upgradeDeflatedZMatrix();
		//更新原始Hessenberg-Trianlge 矩阵使用Z矩阵右乘H-T矩阵对
		updateHTMatrixByZ();
		//更新总体转换矩阵
		updateZMatrixTotal();

		//打印总体转换矩阵乘积测试
		//showQxOPxZ();

		//Hessenberg矩阵进行整形，将极小值变为0元
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

//获取H-T矩阵对
BasicMatrix* GeneralizedEigenSolverForReal::getHessenbergMatrix()
{
	return this->p_OpMatrix_A;
};
BasicMatrix* GeneralizedEigenSolverForReal::getTriangleMatrix()
{
	return this->p_OpMatrix_B;
};

//获取Q\Z 综合转换矩阵
BasicMatrix* GeneralizedEigenSolverForReal::getQMatrix_Total()
{
	return this->p_QMatrix_Total;
};
BasicMatrix* GeneralizedEigenSolverForReal::getZMatrix_Total()
{
	return this->p_ZMatrix_Total;
};
