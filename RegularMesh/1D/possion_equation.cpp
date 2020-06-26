/**
 * @file   possion_equation.cpp
 * @author StudyBo <studybo@ubuntu18-04>
 * @date   Tue Jun 23 21:07:49 2020
 * 
 * @brief  Solve 1D possion equation, interval case.
 * 
 * 
 */

#include <iostream>
#include <cmath>

#include <AFEPack/AMGSolver.h>
#include <AFEPack/Geometry.h>
#include <AFEPack/TemplateElement.h>
#include <AFEPack/FEMSpace.h>
#include <AFEPack/Operator.h>
#include <AFEPack/Functional.h>
#include <AFEPack/EasyMesh.h>
#include <AFEPack/SparseMatrixTool.h>

#include <lac/sparse_matrix.h>
#include <lac/sparsity_pattern.h>
#include <lac/sparse_ilu.h>
#include <lac/vector.h>
#include <lac/full_matrix.h>
#include <lac/solver_cg.h>
#include <lac/sparse_mic.h>
#include <lac/sparse_decomposition.h>

#define PI (4.0 * atan(1.0))

double u (const double* p)
{
    return 0;
}

double f(const double* p)
{
    return 2;
}

/// 用来储存边界自由度编号的向量。
std::vector<unsigned int> BndDOF;
/// 从左到右，对 n 个区间单元中的第 i 个中的两个自由度进行编号， 0->1。
int eledof_1D(int i, int k)
{
    return (i + k);
}
/// n 个单元中第 i 个区间单元中两个节点的坐标
AFEPack::Point<1> elevtx_1D(int n, int i, int k)
{
    double x0 = 0.0;
    double x1 = 1.0;
    AFEPack::Point<1> pnt;
    double h = (x1 - x0) / n;
    pnt[0] = (i + k) * h;
    return pnt;
}
/// 存入边界自由度。
int Record_BndDOF(int n)
{
    BndDOF.push_back(0);
    BndDOF.push_back(n);
    return 0;
}

int main(int argc, char* argv[])
{
    /// 从 AFEPack 中读入 1D 模板单元格信息，基函数信息和坐标变换信息。
    TemplateGeometry<1> interval_template_geometry;
    interval_template_geometry.readData("interval.tmp_geo");
    CoordTransform<1, 1> interval_coord_transform;
    interval_coord_transform.readData("interval.crd_trs");
    TemplateDOF<1> interval_template_dof(interval_template_geometry);
    /// 一次元。
    interval_template_dof.readData("interval.1.tmp_dof");
    BasisFunctionAdmin<double, 1, 1> interval_basis_function(interval_template_dof);
    interval_basis_function.readData("interval.1.bas_fun");
    TemplateElement<double, 1, 1> template_element;
    ///模板单元初始化。
    template_element.reinit(interval_template_geometry,
	                    interval_template_dof,
	                    interval_coord_transform,
	                    interval_basis_function);
    double volume = template_element.volume();
    /// 取 2 阶代数精度。
    const QuadratureInfo<1>& quad_info = template_element.findQuadratureInfo(2);
    /// 积分点个数
    int n_quadrature_point = quad_info.n_quadraturePoint();
    /// 积分点
    std::vector<AFEPack::Point<1> > q_point = quad_info.quadraturePoint();
    /// 自由度个数
    int n_element_dof = template_element.n_dof();
    /// 基函数个数
    int n_bas = template_element.basisFunction().size();
    
    /// 产生一个具体单元顶点的缓存。一个区间的 2 个顶点。
    /// 这里其实是这 gv 表示全局的区间坐标，就是在物理计算
    /// 区域内一个单元的顶点坐标；而 lv 表示局部的区间坐标
    /// ，即参考单元的坐标。沿用了 AFEPack 的配置，是固定
    /// 的 [-1] - [1]。

    /// 全局点
    std::vector<AFEPack::Point<1> > gv(2);
    /// 观察一下模板中的自由度，基函数和基函数在具体积分点的取值情况，可从模板单元中读取。
    TemplateGeometry<1> &geo = template_element.geometry();
    const std::vector<AFEPack::Point<1> > &lv = geo.vertexArray();
    int n_vtx = geo.n_geometry(0);
    /// 设置实际的计算区间端点。
    double x0 = 0.0;
    double x1 = 1.0;
    /// 设置剖分段数。
    int n = 20;
    /// 自由度总数
    int dim = n + 1;
    /// 右端项初始化
    Vector<double> rhs(dim);
    /// 储存边界自由度。
    Record_BndDOF(n);
    /// 每行的非零元个数，每行最多3个非零元。
    std::vector<unsigned int> NoZeroPerRow(dim, 3);
    /// 区间端点所在行只有两个非零元。
    NoZeroPerRow[0] = 2;
    NoZeroPerRow[n] = 2;
    /// 建立稀疏矩阵模板。
    SparsityPattern sp_stiff_matrix(dim, NoZeroPerRow);
    /// 填充非零元素对应的行索引和列索引。
    int n_dof = template_element.n_dof();
    for (int i = 0; i < n; i++)
    {
	for (int dof1 = 0; dof1 < n_dof; dof1++)
	    for (int dof2 = 0; dof2 < n_dof; dof2++)
		sp_stiff_matrix.add(eledof_1D(i, dof1), eledof_1D(i, dof2));
    }
    /// 稀疏矩阵模板生成。
    sp_stiff_matrix.compress();
    /// 刚度矩阵初始化。
    SparseMatrix<double> stiff_mat(sp_stiff_matrix);
    /// 生成节点，单元，刚度矩阵和右端项。
    double h = (x1 - x0) / n;
    for (int i = 0; i < n; i++)
    {
	/// 实际单元顶点坐标
	for (int k = 0; k < n_vtx; k++)
	    gv[k] = elevtx_1D(n, i, k);

	/// 积分
	for (int l = 0; l < n_quadrature_point; l++)
	{
	    /// 积分点全局坐标
	    auto point = interval_coord_transform.local_to_global(q_point, lv, gv);
	    /// 积分点权重， Jacobi 变换系数。
	    double Jxy = quad_info.weight(l) * interval_coord_transform.local_to_global_jacobian(q_point[l], lv, gv) * volume;
	    
	    for (int base1 = 0; base1 < n_dof; base1++)
	    {
		for (int base2 = 0; base2 < n_dof; base2++)
		    stiff_mat.add(eledof_1D(i, base1), eledof_1D(i, base2),
				  Jxy * innerProduct(interval_basis_function[base1].gradient(point[l], gv), interval_basis_function[base2].gradient(point[l], gv)));
		/// 右端项
		rhs(eledof_1D(i, base1)) += Jxy * f(point[l]) * interval_basis_function[base1].value(point[l], gv);
	    }
	}
    }
    /// 边界条件处理
    for (int i = 0; i < BndDOF.size(); i++)
    {
	int bnd_dof = BndDOF[i];
	double x = x0 + h * bnd_dof;
	SparseMatrix<double>::iterator row_iterator = stiff_mat.begin(bnd_dof);
	double diag = row_iterator->value();
	AFEPack::Point<1> bnd_point;
	bnd_point[0] = x;
	double bnd_value = u(bnd_point);
	rhs(bnd_dof) = diag * bnd_value;
	for (++row_iterator; row_iterator != stiff_mat.end(bnd_dof); ++row_iterator)
	{
	    row_iterator->value() = 0.0;
	    int k = row_iterator->column();
	    SparseMatrix<double>::iterator col_iterator = stiff_mat.begin(k);
	    SparseMatrix<double>::iterator col_end = stiff_mat.end(k);
	    for (++col_iterator; col_iterator != col_end; ++col_iterator)
		if (col_iterator->column() == bnd_dof)
		    break;
	    if (col_iterator == col_end)
	    {
		std::cerr << "Error!" << std::endl;
		exit(-1);
	    }
	    rhs(k) -= col_iterator->value() * bnd_value;
	    col_iterator->value() = 0.0;
	}
    }
    /// 用代数多重网格 AMG 计算线性方程。
    AMGSolver solver(stiff_mat);
    /// 这里设置线性求解器的收敛判定为机器 epsilon 乘以矩阵的阶数，也就是自由度总数。这个参数基本上是理论可以达到的极限。
    Vector<double> solution(dim);
    double tol = std::numeric_limits<double>::epsilon() * dim;
    solver.solve(solution, rhs, tol, 10000);
    /// 输出到output.m，用Matlab或Octave运行，得到计算结果。
    std::ofstream fs;
    fs.open("output.m");
    fs << "x = " << x0 << ": 1/" << n << " : " << x1 << ";" << std::endl;
    fs << "u = [";
    for (int i = 0; i < dim; i++)
    {
	fs << " " << solution[i] << ",";
    }
    fs << "];" << std::endl;
    fs << "plot(x, u);" << std::endl;
    return 0;   
}

    
    
