/**
 * @file   possion_equation_Q1.cpp
 * @author StudyBo <studybo@ubuntu18-04>
 * @date   Tue Jun 23 21:07:49 2020
 * 
 * @brief  Solve 2D possion equation, Q1 case.
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

double u(const double *p)
{
    return sin(PI * p[0]) * sin(PI * p[1]);
}

double f(const double *p)
{
    return 2 * PI * PI * u(p);
}

/// 从 j 行 i 列，每一行 n 个网格的 Q1 剖分中映射 (i,j) 单元的第 k 个自由度编号。
int Q1_ele2dof(int n, int j, int i, int k)
{
    int idx = -1;
    switch (k)
    {
    case 0:
	idx = j * (n + 1) + i;
	break;
    case 1:
	idx = j * (n + 1) + i + 1;
	break;
    case 2:
	idx = (j + 1) * (n + 1) + i + 1;
	break;
    case 3:
	idx = (j + 1) * (n + 1) + i;
	break;
    default:
	std::cerr << "Dof. no. error!" << std::endl;
	exit(-1);
    }
    return idx;
}

/// 每个自由度对应的顶点坐标
AFEPack::Point<2> Q1_ele2vtx(int n, int j, int i, int k)
{
    double x0 = 0.0;
    double x1 = 1.0;
    double y0 = 0.0;
    double y1 = 1.0;
    AFEPack::Point<2> pnt;
    switch(k)
    {
    case 0:
	pnt[0] = ((n - i) * x0 + i * x1) / n;
	pnt[1] = ((n - j) * y0 + j * y1) / n;
	break;
    case 1:
	pnt[0] = ((n - i - 1) * x0 + (i + 1) * x1) / n;
	pnt[1] = ((n - j) * y0 + j * y1) / n;
	break;
    case 2:
	pnt[0] = ((n - i - 1) * x0 + (i + 1) * x1) / n;
	pnt[1] = ((n - j - 1) * y0 + (j + 1) * y1) / n;
	break;
    case 3:
	pnt[0] = ((n - i) * x0 + i * x1) / n;
	pnt[1] = ((n - j - 1) * y0 + (j + 1) * y1) / n;
	break;
    default:
	std::cerr << "Element vertex no. error!" << std::endl;
	exit(-1);
    }
    return pnt;
}

int main(int argc, char* argv[])
{
    /// 从 AFEPack 中读入 Q1 模板单元格信息，基函数信息和坐标变换信息。
    TemplateGeometry<2> rectangle_template_geometry;
    rectangle_template_geometry.readData("rectangle.tmp_geo");
    CoordTransform<2, 2> rectangle_coord_transform;
    rectangle_coord_transform.readData("rectangle.crd_trs");
    TemplateDOF<2> rectangle_template_dof(rectangle_template_geometry);
    ///一次元
    rectangle_template_dof.readData("rectangle.1.tmp_dof");
    BasisFunctionAdmin<double, 2, 2> rectangle_basis_function(rectangle_template_dof);
    rectangle_basis_function.readData("rectangle.1.bas_fun");
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
    int n_element_dof = template_element.n_dof();
    /// 基函数个数
    int n_bas = template_element.basisFunction().size();

    /// 产生一个具体单元顶点的缓存。一个矩形的 4 个顶点。这里其实是这
    /// 四个顶点正好是 Q1 单元的 4 个单元内自由度。gv 表示全局的矩形坐
    /// 标，就是在物理计算区域内一个网格的顶点坐标；而 lv 表示局部的矩
    /// 形坐标，即参考单元的坐标。沿用了 AFEPack 的配置，是固定的 [-1,
    /// -1]-[-1, 1]-[1, 1]-[-1, 1]。
    
    /// 全局点
    std::vector<AFEPack::Point<2> > gv(4);
    /// 观察一下模板单元中的自由度，基函数和基函数在具体积分点的取值情况，可从模板单元中读取到。
    TemplateGeometry<2> &geo = template_element.geometry();
    const std::vector<AFEPack::Point<2> > &lv = geo.vertexArray();
    int n_vtx = geo.n_geometry(0);
    /// 设置实际的计算矩形边界。
    double x0 = 0.0;
    double x1 = 1.0;
    double y0 = 0.0;
    double y1 = 1.0;
    /// 设置剖分段数。
    int n = 20;
    /// 自由度总数
    int dim = (n + 1) * (n + 1);
    /// rhs是右端项
    Vector<double> rhs(dim);
    /// 每行对应的非零元个数 NoZeroPerRow，每行最多9个非零元。
    std::vector<unsigned int> NoZeroPerRow(dim, 9);
    /// 角点自由度所在行只有 4 个非零元。
    NoZeroPerRow[0] = 4;
    NoZeroPerRow[dim - 1] = 4;
    NoZeroPerRow[n] = 4;
    NoZeroPerRow[dim - n - 1] = 4;
    /// 非角点的边界自由度只有 6 个非零元。
    for (int i = 1; i < n; i++)
    {
	NoZeroPerRow[i] = 6;
	NoZeroPerRow[dim - 1 - i] = 6;
	NoZeroPerRow[i * (n + 1)] = 6;
	NoZeroPerRow[i * n + i + n] = 6;
    }
    /// 建立稀疏矩阵模板。
    SparsityPattern sp_stiff_matrix(dim, NoZeroPerRow);
    /// 填充非零元素对应的行索引和列索引，遍历顺序按照单元的顺序。
    for (int j = 0; j < n; j++)
    {
	for (int i = 0; i < n; i++)
	{
	    int n_dof = template_element.n_dof();
	    for (int dof1 = 0; dof1 < n_dof; dof1++)
		for (int dof2 = 0; dof2 < n_dof; dof2++)
		{
		    sp_stiff_matrix.add(Q1_ele2dof(n, j, i, dof1),
					Q1_ele2dof(n, j, i, dof2));
		}
	}
    }
    /// 稀疏矩阵模板生成。
    sp_stiff_matrix.compress();
    /// 刚度矩阵初始化。
    SparseMatrix<double> stiff_mat(sp_stiff_matrix);
    /// 生成节点，单元，刚度矩阵和右端项。
    double h = (x1 - x0) / n;
    for (int j = 0; j < n; j++)
	for (int i = 0; i < n; i++)
	{
	    /// 实际单元顶点坐标。
	    for (int k = 0; k < n_vtx; k++)
		gv[k] = Q1_ele2vtx(n ,j, i, k);

	    /// 积分
	    for (int l = 0; l < n_quadrature_point; l++)
	    {
		/// 积分点全局坐标。
		auto point = rectangle_coord_transform.local_to_global(q_point, lv, gv);
		/// 积分点权重，Jacobi变换系数。
		double Jxy = quad_info.weight(l) * rectangle_coord_transform.local_to_global_jacobian(q_point[l], lv, gv) * volume;

		for (int base1 = 0; base1 < template_element.n_dof(); base1++)
		{
		    for (int base2 = 0; base2 < template_element.n_dof(); base2++)
			stiff_mat.add(Q1_ele2dof(n, j, i, base1),
			              Q1_ele2dof(n, j, i, base2),
			              Jxy * innerProduct(rectangle_basis_function[base1].gradient(point[l], gv), rectangle_basis_function[base2].gradient(point[l], gv)));
		    /// 右端项
		    rhs(Q1_ele2dof(n, j, i, base1)) += Jxy * f(point[l]) * rectangle_basis_function[base1].value(point[l], gv);
		}
	    }
	}
    /// 边界条件处理。
    for (int index = 0; index < dim; index++)
    {
	/// 用非零元个数判断边界。
	if (NoZeroPerRow[index] == 4 || NoZeroPerRow[index] == 6)
	{
	    /// 首先计算该节点的实际坐标。
	    int x_num = index % (n + 1);
	    int y_num = index / (n + 1);
	    double x = x_num * h;
	    double y = y_num * h;
	    SparseMatrix<double>::iterator row_iterator = stiff_mat.begin(index);
	    SparseMatrix<double>::iterator row_end = stiff_mat.end(index);
	    double diag = row_iterator->value();
	    AFEPack::Point<2> bnd_point;
	    bnd_point[0] = x;
	    bnd_point[1] = y;
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
    }

    /// 用代数多重网格 AMG 计算线性方程。
    AMGSolver solver(stiff_mat);
    /// 这里设置线性求解器的收敛判定为机器 epsilon 乘以矩阵的阶数，也就是自由度总数。这个参数基本上是理论可以达到的极限。
    Vector<double> solution(dim);
    double tol = std::numeric_limits<double>::epsilon() * dim;
    solver.solve(solution, rhs, tol, 10000);
    std::ofstream fs;
    /// 输出到 output.m
    fs.open("output.m");
    fs << "x = 0:1/" << n << ":1;" << std::endl;
    fs << "y = 0:1/" << n << ":1;" << std::endl;
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
    fs << "surf(X, Y, u);" << std::endl;
    return 0;
}

    


    
