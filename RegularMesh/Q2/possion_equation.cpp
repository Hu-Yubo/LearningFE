/**
 * @file   possion_equation.cpp
 * @author StudyBo <studybo@ubuntu18-04>
 * @date   Tue Jun 30 09:44:27 2020
 * 
 * @brief  Solve 2D possion equation, Q2 case.
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

/// 实际计算的矩形边界。
double xa = 0.0;
double xb = 1.0;
double ya = 0.0;
double yb = 1.0;
/// 用来存放边界自由度编号的向量。
std::vector<unsigned int> BndDof;

double u(const double* p)
{
    return sin(PI * p[0]) * sin(PI * p[1]);
}

double f(const double* p)
{
    return 2 * PI * PI * u(p);
}

/// 从 j 行 i 列，每一行 n 个网格的 Q2 剖分中映射 (i,j) 单元的第 k 个自由度编号。
int Q2_ele2dof(int n, int j, int i, int k)
{
    int dof = -1;
    switch(k)
    {
    case 0:
	dof = 2 * j * (2 * n + 1) + i * 2;
	break;
    case 1:
	dof = 2 * j * (2 * n + 1) + i * 2 + 2;
	break;
    case 2:
	dof = 2 * (j + 1)  * (2 * n + 1) + i * 2 + 2;
	break;
    case 3:
	dof = 2 * (j + 1) * (2 * n + 1) + i * 2;
	break;
    case 4:
	dof = 2 * j * (2 * n + 1) + i * 2 + 1;
	break;
    case 5:
	dof = (2 * j + 1) * (2 * n + 1) + i * 2 + 2;
	break;
    case 6:
	dof = 2 * (j + 1) * (2 * n + 1) + i * 2 + 1;
	break;
    case 7:
	dof = (2 * j + 1) * (2 * n + 1) + i * 2;
	break;
    case 8:
	dof = (2 * j + 1) * (2 * n + 1) + i * 2 + 1;
	break;	
    default:
	std::cerr << "Dof. no. error!" << std::endl;
	exit(-1);
    }
    return dof;
}

/// 每个编号为 dof 的自由度所对应的节点
AFEPack::Point<2> Dof_to_vtx(int n, int dof)
{
    AFEPack::Point<2> pnt;
    /// j 行 i 列
    int j = dof / (2 * n + 1);
    int i = dof % (2 * n + 1);
    pnt[0] = i * (xb - xa) / (n * 2);
    pnt[1] = j * (yb - ya) / (n * 2);
    return pnt;
}

/// 从 j 行 i 列，每一行 2n 个网格的 Q2 剖分中映射 (i,j) 单元的第 k 个顶点坐标。
AFEPack::Point<2> Q2_ele2vtx(int n, int j, int i, int k)
{
    AFEPack::Point<2> pnt;
    int dof = Q2_ele2dof(n, j, i, k);
    pnt = Dof_to_vtx(n, dof);
    return pnt;
}

/// 储存所有边界自由度编号。
int Store_BndDof(int n)
{
    for (int i = 0; i < 2 * n + 1; i++)
    {
	BndDof.push_back(i);
	BndDof.push_back((2 * n + 1) * (2 * n + 1) - 1 - i);
    }
    if (n >= 1)
    {
	for (int i = 1; i < 2 * n; i++)
	{
	    BndDof.push_back(i * (2 * n + 1));
	    BndDof.push_back(i * (2 * n + 1) + 2 * n);
	}
    }
    return 0;
}

int main(int argc, char* argv[])
{
    /// 从 AFEPack 中读入 Q2 模板单元格信息，基函数信息和坐标变换信息。
    TemplateGeometry<2> rectangle_template_geometry;
    rectangle_template_geometry.readData("rectangle.tmp_geo");
    CoordTransform<2, 2> rectangle_coord_transform;
    rectangle_coord_transform.readData("rectangle.crd_trs");
    TemplateDOF<2> rectangle_template_dof(rectangle_template_geometry);
    ///一次元
    rectangle_template_dof.readData("rectangle.2.tmp_dof");
    BasisFunctionAdmin<double, 2, 2> rectangle_basis_function(rectangle_template_dof);
    rectangle_basis_function.readData("rectangle.2.bas_fun");
    TemplateElement<double, 2, 2> template_element;
    /// 模板单元初始化
    template_element.reinit(rectangle_template_geometry,
	                    rectangle_template_dof,
	                    rectangle_coord_transform,
	                    rectangle_basis_function);
    double volume = template_element.volume();
    /// 取 4 次代数精度。
    const QuadratureInfo<2>& quad_info = template_element.findQuadratureInfo(4);
    /// 积分点个数
    int n_quadrature_point = quad_info.n_quadraturePoint();
    /// 积分点
    std::vector< AFEPack::Point<2> > q_point = quad_info.quadraturePoint();
    int n_dof = template_element.n_dof();
    /// 基函数个数
    int n_bas = template_element.basisFunction().size();

    /// 观察一下模板单元中的自由度，基函数和基函数在具体积分点的取值情况，可从模板单元中读取到。
    int n_vtx = 9;
    /// 局部点，参考单元的 9 个节点， 取点参照 rectangle.2.bas_fun。
    std::vector<AFEPack::Point<2> > lv(n_vtx);
    lv[0][0] = -1.0;  lv[0][1] = -1.0;
    lv[1][0] = 1.0;   lv[1][1] = -1.0;
    lv[2][0] = 1.0;   lv[2][1] = 1.0;
    lv[3][0] = -1.0;  lv[3][1] = 1.0;
    lv[4][0] = 0.0;   lv[4][1] = -1.0;
    lv[5][0] = 1.0;   lv[5][1] = 0.0;
    lv[6][0] = 0.0;   lv[6][1] = 1.0;
    lv[7][0] = -1.0;  lv[7][1] = 0.0;
    lv[8][0] = 0.0;   lv[8][1] = 0.0;
    /// 全局点
    std::vector<AFEPack::Point<2> > gv(n_vtx);
    /// 设置剖分段数。
    int n = 50;
    /// 自由度总数
    int dim = (2 * n + 1) * (2 * n + 1);
    /// 右端项
    Vector<double> rhs(dim);
    /// 记得先储存边界自由度！！！！！！！！
    Store_BndDof(n);
    /// 稀疏矩阵中每行对应的非零元个数，非零元个数最多为 25。
    std::vector<unsigned int> NoZeroPerRow(dim, 25);
    /// 听取了 LSJ 的建议，不需要确定具体的非零元个数，直接填充索引。
    /// 建立稀疏矩阵模板。
    SparsityPattern sp_stiff_matrix(dim, NoZeroPerRow);
    /// 填充非零元素对应的行索引和列索引，遍历顺序按照单元的顺序。
    for (int j = 0; j < n; j++)
    {
	for (int i = 0; i < n ; i++)
	{
	    for (int dof1 = 0; dof1 < n_dof; dof1++)
		for (int dof2 = 0; dof2 < n_dof; dof2++)
		{
		    sp_stiff_matrix.add(Q2_ele2dof(n, j, i, dof1),
			                Q2_ele2dof(n, j, i, dof2));
		}
	}
    }
    /// 稀疏矩阵模板生成。
    sp_stiff_matrix.compress();
    sp_stiff_matrix.print(std::cout);
    /// 刚度矩阵初始化。
    SparseMatrix<double> stiff_mat(sp_stiff_matrix);
    /// 生成节点，单元，刚度矩阵和右端项。
    for (int j = 0; j < n; j++)
	for (int i = 0; i < n; i++)
	{
	    /// 实际单元顶点坐标。。
	    for (int k = 0; k < n_vtx; k++)
		gv[k] = Q2_ele2vtx(n, j, i, k);
	    /// 积分
	    for (int l = 0; l < n_quadrature_point; l++)
	    {
		/// 积分点全局坐标。
		auto point = rectangle_coord_transform.local_to_global(q_point, lv, gv);
		/// 积分点权重，Jacobi变换系数。
		double Jxy = quad_info.weight(l) * rectangle_coord_transform.local_to_global_jacobian(q_point[l], lv, gv) * volume;

		for (int base1 = 0; base1 < n_dof; base1++)
		{
		    for (int base2 = 0; base2 < n_dof; base2++)
			stiff_mat.add(Q2_ele2dof(n, j, i, base1),
			              Q2_ele2dof(n, j, i, base2),
			              Jxy * innerProduct(rectangle_basis_function[base1].gradient(point[l], gv), rectangle_basis_function[base2].gradient(point[l], gv)));
		    /// 右端项
		    rhs(Q2_ele2dof(n, j, i, base1)) += Jxy * f(point[l]) * rectangle_basis_function[base1].value(point[l], gv);
		}
	    }
	}
    /// 处理边界条件
    for (int i = 0; i < BndDof.size(); i++)
    {
	int index = BndDof[i];
	SparseMatrix<double>::iterator row_iterator = stiff_mat.begin(index);
	SparseMatrix<double>::iterator row_end = stiff_mat.end(index);
	double diag = row_iterator->value();
	AFEPack::Point<2> bnd_point;
	bnd_point = Dof_to_vtx(n, index);
	double bnd_value = u(bnd_point);
	rhs(index) = diag * bnd_value;
	for (++row_iterator; row_iterator != row_end; ++row_iterator)
	{
	    row_iterator->value() = 0.0;
	    int k = row_iterator->column();
	    SparseMatrix<double>::iterator col_iterator = stiff_mat.begin(k);
	    SparseMatrix<double>::iterator col_end = stiff_mat.end(k);
	    for (++col_iterator; col_iterator != col_end; ++col_iterator)
		if (col_iterator->column() == index)
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
    /// 输出到 output.m
    std::ofstream fs;
    fs.open("output.m");
    fs << "x = " << xa << ":" << xb - xa << "/" << 2*n << ":" << xb << ";" << std::endl;
    fs << "y = " << ya << ":" << yb - ya << "/" << 2*n << ":" << yb << ";" << std::endl;
    fs << "[X, Y] = meshgrid(x, y);" << std::endl;
    fs << "u = [";
    for (int j = 0; j < 2 * n + 1; j++)
    {
	for (int i = 0; i < 2 * n + 1; i++)
	{
	    fs << solution[i + (2 * n + 1) * j] << " , ";
	}
	fs << ";" << std::endl;
    }
    fs << "];" << std::endl;
    fs << "surf(X, Y, u);" << std::endl;

    /// 计算 L2 误差。
    double error = 0;
    for (int dof = 0; dof < dim; dof++)
    {
	AFEPack::Point<2> pnt;
	pnt = Dof_to_vtx(n, dof);
	double d = (u(pnt) - solution[dof]);
	error += d * d;
    }
    error = std::sqrt(error);
    std::cerr << "\nL2 error = " << error << ", tol = " << tol << std::endl;
    return 0;
}
