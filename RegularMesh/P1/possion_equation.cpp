/**
 * @file   possion_equation.cpp
 * @author StudyBo <studybo@ubuntu18-04>
 * @date   Fri Jun 26 17:11:00 2020
 * 
 * @brief  Solve 2D possion equation, P1 case.
 * 
 * 
 */

#include <iostream>
#include <cmath>
#include <iomanip>

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

/* 如何划分单元格？                     
 * 首先将实际计算区域划分为规则的矩形网格，对于每个矩形网格，
 * 将其第0个和第2个顶点连线，划分为两个三角形，从左至右排序。
 * 2  ——  3 
 * |     /|    如图左所示，0-3-2, 0-1-3 构成两个 P1 单元片，
 * |    / |    且 0-3-2 编号在前。
 * |   /  |
 * |  /   |
 * | /    |
 * |/     |
 * 0  ——  1
 *
 */

/// 从 j 行 i 列，每一行 2n 个网格的 P1 剖分中映射 (i,j) 单元的第 k 个自由度编号。
int P1_ele2dof(int n, int j, int i, int k)
{
    int dof = -1;
    if ( i % 2 == 0)
    {
	switch(k)
	{
	case 0:
	    dof = j * (n + 1) + i / 2;
	    break;
	case 1:
	    dof = (j + 1) * (n + 1) + i / 2 + 1;
	    break;
	case 2:
	    dof = (j + 1) * (n + 1) + i / 2;
	    break;
	default:
	    std::cerr << "Dof. no. error!" << std::endl;
	    exit(-1);
	}
    }
    else
    {
	switch(k)
	{
	case 0:
	    dof = j * (n + 1) + i / 2;
	    break;
	case 1:
	    dof = j * (n + 1) + i / 2 + 1;
	    break;
	case 2:
	    dof = (j + 1) * (n + 1) + i / 2 + 1;
	    break;
	default:
	    std::cerr << "Dof. no. error!" << std::endl;
	    exit(-1);
	}
    }
    return dof;
}

/// 从 j 行 i 列，每一行 2n 个网格的 P1 剖分中映射 (i,j) 单元的第 k 个顶点坐标。
AFEPack::Point<2> P1_ele2vtx(int n, int j, int i, int k)
{
    AFEPack::Point<2> pnt;
    if (i % 2 == 0)
    {
	switch (k)
	{
	case 0:
	    pnt[0] = ((n - i / 2) * xa + i / 2 * xb) / n;
	    pnt[1] = ((n - j) * ya + j * yb) / n;
	    break;
	case 1:
	    pnt[0] = ((n - i / 2 - 1) * xa + (i / 2 + 1) * xb) / n;
	    pnt[1] = ((n - j - 1) * ya + (j + 1) * yb) / n;
	    break;
	case 2:
	    pnt[0] = ((n - i / 2) * xa + i / 2 * xb) / n;
	    pnt[1] = ((n - j - 1) * ya + (j + 1) * yb) / n;
	    break;
	default:
	    std::cerr << "Element vertex no. error!" << std::endl;
	    exit(-1);
	}
    }
    else
    {
	switch (k)
	{
	case 0:
	    pnt[0] = ((n - i / 2) * xa + i / 2 * xb) / n;
	    pnt[1] = ((n - j) * ya + j * yb) / n;
	    break;
	case 1:
	    pnt[0] = ((n - i / 2 - 1) * xa + (i / 2 + 1) * xb) / n;
	    pnt[1] = ((n - j) * ya + j * yb) / n;
	    break;
	case 2:
	    pnt[0] = ((n - i / 2 - 1) * xa + (i / 2 + 1) * xb) / n;
	    pnt[1] = ((n - j - 1) * ya + (j + 1) * yb) / n;
	    break;
	default:
	    std::cerr << "Element vertex no. error!" << std::endl;
	    exit(-1);
	}
    }
    return pnt;
}

/// 储存所有边界自由度编号。
int Store_BndDof(int n)
{
    for (int i = 0; i < n + 1; i++)
    {
	BndDof.push_back(i);
	BndDof.push_back((n + 1) * (n + 1) - 1 - i);
    }
    if (n > 1)
    {
	for (int i = 1; i < n; i++ )
	{
	    BndDof.push_back(i * (n + 1));
	    BndDof.push_back(i * (n + 1) + n);
	}
    }
    return 0;
}

/// 每个编号为 dof 的自由度所对应的节点
AFEPack::Point<2> Dof_to_vtx(int n, int dof)
{
    AFEPack::Point<2> pnt;
    /// j 行 i 列
    int j = dof / (n + 1);
    int i = dof % (n + 1);
    pnt[0] = i * (xb - xa) / n;
    pnt[1] = j * (yb - ya) / n;
    return pnt;
}

int main(int argc, char* argv[])
{
    /// 从 AFEPack 中读入 P1 模板单元格信息，基函数信息和坐标变换信息。
    TemplateGeometry<2> triangle_template_geometry;
    triangle_template_geometry.readData("triangle.tmp_geo");
    CoordTransform<2, 2> triangle_coord_transform;
    triangle_coord_transform.readData("triangle.crd_trs");
    TemplateDOF<2> triangle_template_dof(triangle_template_geometry);
    /// 一次元
    triangle_template_dof.readData("triangle.1.tmp_dof");
    BasisFunctionAdmin<double, 2, 2> triangle_basis_function(triangle_template_dof);
    triangle_basis_function.readData("triangle.1.bas_fun");
    TemplateElement<double, 2, 2> template_element;
    /// 模板单元初始化
    template_element.reinit(triangle_template_geometry,
			    triangle_template_dof,
			    triangle_coord_transform,
			    triangle_basis_function);
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
    
    /// 产生一个具体单元顶点的缓存。一个三角形的 3 个顶点。这里其实是这
    /// 三个顶点正好是 P1 单元的 3 个单元内自由度。gv 表示全局的三角形坐
    /// 标，就是在物理计算区域内一个网格的顶点坐标；而 lv 表示局部的三角
    /// 形坐标，即参考单元的坐标。沿用了 AFEPack 的配置，是固定的 [0,
    /// 0]-[0, 1]-[1, 0]。

    /// 全局点
    std::vector<AFEPack::Point<2> > gv(3);
    /// 观察一下模板单元中的自由度，基函数和基函数在具体积分点的取值情况，可从模板单元中读取到。
    TemplateGeometry<2> &geo = template_element.geometry();
    const std::vector<AFEPack::Point<2> > &lv = geo.vertexArray();
    int n_vtx = geo.n_geometry(0);
    /// 设置剖分段数。
    int n = 50;
    /// 自由度总数
    int dim = (n + 1) * (n + 1);
    /// 右端项
    Vector<double> rhs(dim);
    /// 记得先储存边界自由度！！！！！！！！
    Store_BndDof(n);
    /// 稀疏矩阵中每行对应的非零元个数，每行最多 7 个非零元。
    std::vector<unsigned int> NoZeroPerRow(dim, 7);
    /// 非角点的边界自由度所在行有 5 个非零元。
    for (int i = 0; i < BndDof.size(); i++)
	NoZeroPerRow[BndDof[i]] = 5;
    /// 角点自由度所在行有 3 或 4 个非零元。
    NoZeroPerRow[0] = 4;
    NoZeroPerRow[dim - 1] = 4;
    NoZeroPerRow[n] = 3;
    NoZeroPerRow[dim - n - 1] = 3;
    /// 建立稀疏矩阵模板。
    SparsityPattern sp_stiff_matrix(dim, NoZeroPerRow);
    /// 填充非零元素对应的行索引和列索引，遍历顺序按照单元的顺序。
    for (int j = 0; j < n; j++)
    {
	for (int i = 0; i < 2*n ; i++)
	{
	    for (int dof1 = 0; dof1 < n_dof; dof1++)
		for (int dof2 = 0; dof2 < n_dof; dof2++)
		{
		    sp_stiff_matrix.add(P1_ele2dof(n, j, i, dof1),
			                P1_ele2dof(n, j, i, dof2));
		}
	}
    }
    /// 稀疏矩阵模板生成。
    sp_stiff_matrix.compress();
    /// 刚度矩阵初始化。
    SparseMatrix<double> stiff_mat(sp_stiff_matrix);
    /// 生成节点，单元，刚度矩阵和右端项。
    // double h = (xb - xa) / n;
    for (int j = 0; j < n; j++)
	for (int i = 0; i < 2 * n; i++)
	{
	    /// 实际单元顶点坐标。。
	    for (int k = 0; k < n_vtx; k++)
		gv[k] = P1_ele2vtx(n ,j, i, k);
	    /// 积分
	    for (int l = 0; l < n_quadrature_point; l++)
	    {
		/// 积分点全局坐标。
		auto point = triangle_coord_transform.local_to_global(q_point, lv, gv);
		/// 积分点权重，Jacobi变换系数。
		double Jxy = quad_info.weight(l) * triangle_coord_transform.local_to_global_jacobian(q_point[l], lv, gv) * volume;

		for (int base1 = 0; base1 < n_dof; base1++)
		{
		    for (int base2 = 0; base2 < n_dof; base2++)
			stiff_mat.add(P1_ele2dof(n, j, i, base1),
			              P1_ele2dof(n, j, i, base2),
			              Jxy * innerProduct(triangle_basis_function[base1].gradient(point[l], gv), triangle_basis_function[base2].gradient(point[l], gv)));
		    /// 右端项
		    rhs(P1_ele2dof(n, j, i, base1)) += Jxy * f(point[l]) * triangle_basis_function[base1].value(point[l], gv);
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
    std::ofstream fs;

/// 输出到 output.m
    /*
    fs.open("output.m");
    fs << "x = " << xa << ":" << xb - xa << "/" << n << ":" << xb << ";" << std::endl;
    fs << "y = " << ya << ":" << yb - ya << "/" << n << ":" << yb << ";" << std::endl;
    fs << "[X, Y] = meshgrid(x, y);" << std::endl;
    fs << "u = [";
    for (int j = 0; j < n + 1; j++)
    {
	for (int i = 0; i < n + 1; i++)
	{
	    fs << solution[i + (n + 1) * j] << " , ";
	}
	fs << ";" << std::endl;
    }
    fs << "];" << std::endl;
    fs << "tri = [";
    for (int j = 0; j < n; j++)
	for (int i = 0; i < 2 * n; i++)
	{
	    for (int k = 0; k < 3; k++)
		fs << P1_ele2dof(n, j, i, k) + 1 << " ";
	    fs << ";" << std::endl;
	}
    fs << "];" << std::endl;
    fs << "vtkwrite(\"possion.vtk\", \"polydata\", \"triangle\", X, Y, u, tri);" << std::endl;
    */

    fs.open("possion.vtk");
    fs << "# vtk DataFile Version 2.0\n";
    fs << "VTK from Cpp\n";
    fs << "ASCII\n";
    fs << "DATASET STRUCTURED_GRID\n";
    fs << "DIMENSIONS " << n+1 << " " << n+1 << " " << 1;
    fs << "\nPOINTS " << dim << " float\n";
    for (int i = 0; i < dim; i++)
    {
	AFEPack::Point<2> P = Dof_to_vtx(n, i);
	fs << std::fixed << std::setprecision(2) << P[0] << " " << P[1] << " " << solution[i] << std::endl;
    }
    fs << std::endl;
    fs << "POINT_DATA " << dim << std::endl;
    fs << "SCALARS Value float\n";
    fs << "LOOKUP_TABLE default\n";
    for (int j = 0; j < dim; j++)
	fs << solution[j] << std::endl;

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
    

    
	    
