/**
 * @file   possion_equation.cpp
 * @author StudyBo <studybo@ubuntu18-04>
 * @date   Thu Jul  2 16:56:41 2020
 * 
 * @brief  Solve 3D possion equation, P1 tetrahedron case
 * 
 * 
 */

#include <iostream>
#include <cmath>
#include <algorithm>

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

#define PI 4.0 * atan(1.0)
/// 实际计算的长方体边界。
double xa = 0.0;
double xb = 1.0;
double ya = 0.0;
double yb = 1.0;
double za = 0.0;
double zb = 1.0;

/// 用来存放边界自由度编号的向量。
std::vector<unsigned int> BndDof;

double u(const double* p)
{
    return sin(PI * p[0]) * sin(PI * p[1]) * sin(PI * p[2]);
}

double f(const double* p)
{
    return 3 * PI * PI * u(p);
}

/* 如何划分单元格？                     
 * 首先将实际计算区域划分为规则的正六面体网格，对于每个正六面体，
 * 将其分为 5 个四面体：0-1-2-5，0-2-3-7，0-4-5-7，2-5-6-7，0-5-2-7。
 * (编号为顶点编号，并非自由度编号！！！)
 *
 *             7 ———————————— 6
 *            /|             /|
 *           / |            / |
 *          4 ———————————— 5  | 
 *          |  |           |  |
 *          |  3 ——————————|— 2
 *          | /            | /
 *          |/             |/
 *          0 ———————————— 1
 */

/// 从 j 行 i 列 k 层，每一行 5n 个 P1 剖分单元中映射 (i,j,k) 单元的第 m 个自由度编号。
/// 自由度编号按照从左至右，从前至后，从下至上的顺序。
int P1_ele3dof(int n, int k, int j, int i, int m)
{
    int dof = -1;
    int number = i % 5;
    i = i / 5;
    switch (number)
    {
    case 0:   /// 0-1-2-5 型单元
	switch (m)
	{
	case 0:
	    dof = k * (n + 1) * (n + 1) + j * (n + 1) + i;
	    break;
	case 1:
	    dof = k * (n + 1) * (n + 1) + j * (n + 1) + i + 1;
	    break;
	case 2:
	    dof = k * (n + 1) * (n + 1) + (j + 1) * (n + 1) + i + 1;
	    break;
	case 3:
	    dof = (k + 1) * (n + 1) * (n + 1) + j * (n + 1) + i + 1;
	    break;
	default:
	    std::cerr << "Dof. no. error!" << std::endl;
	    exit(-1);
	}
	break;
    case 1:   /// 0-2-3-7 型单元
	switch (m)
	{
	case 0:
	    dof = k * (n + 1) * (n + 1) + j * (n + 1) + i;
	    break;
	case 1:
	    dof = k * (n + 1) * (n + 1) + (j + 1) * (n + 1) + i + 1;
	    break;
	case 2:
	    dof = k * (n + 1) * (n + 1) + (j + 1) * (n + 1) + i;
	    break;
	case 3:
	    dof = (k + 1) * (n + 1) * (n + 1) + (j + 1) * (n + 1) + i;
	    break;
	default:
	    std::cerr << "Dof. no. error!" << std::endl;
	    exit(-1);
	}
	break;
    case 2:  /// 0-4-5-7 型单元
	switch (m)
	{
	case 0:
	    dof = k * (n + 1) * (n + 1) + j * (n + 1) + i;
	    break;
	case 1:
	    dof = (k + 1) * (n + 1) * (n + 1) + j * (n + 1) + i;
	    break;
	case 2:
	    dof = (k + 1) * (n + 1) * (n + 1) + j * (n + 1) + i + 1;
	    break;
	case 3:
	    dof = (k + 1) * (n + 1) * (n + 1) + (j + 1) * (n + 1) + i;
	    break;
	default:
	    std::cerr << "Dof. no. error!" << std::endl;
	    exit(-1);
	}
	break;
    case 3:  /// 2-5-6-7 型单元
	switch (m)
	{
	case 0:
	    dof = k * (n + 1) * (n + 1) + (j + 1) * (n + 1) + i + 1;
	    break;
	case 1:
	    dof = (k + 1) * (n + 1) * (n + 1) + j * (n + 1) + i + 1;
	    break;
	case 2:
	    dof = (k + 1) * (n + 1) * (n + 1) + (j + 1) * (n + 1) + i + 1;
	    break;
	case 3:
	    dof = (k + 1) * (n + 1) * (n + 1) + (j + 1) * (n + 1) + i;
	    break;
	default:
	    std::cerr << "Dof. no. error!" << std::endl;
	    exit(-1);
	}
	break;
    case 4:  /// 0-5-2-7 型单元
	switch (m)
	{
	case 0:
	    dof = k * (n + 1) * (n + 1) + j * (n + 1) + i;
	    break;
	case 1:
	    dof = (k + 1) * (n + 1) * (n + 1) + j * (n + 1) + i + 1;
	    break;
	case 2:
	    dof = k * (n + 1) * (n + 1) + (j + 1) * (n + 1) + i + 1;
	    break;
	case 3:
	    dof = (k + 1) * (n + 1) * (n + 1) + (j + 1) * (n + 1) + i;
	    break;
	default:
	    std::cerr << "Dof. no. error!" << std::endl;
	    exit(-1);
	}
	break;
    default:
	std::cerr << "Dof. no. error!" << std::endl;
	exit(-1);
    }
    return dof;
}

/// 每个编号为 dof 的自由度所对应的节点
AFEPack::Point<3> Dof_to_vtx(int n, int dof)
{
    AFEPack::Point<3> pnt;
    /// k 层 j 行 i 列
    int k = dof / (n * n + 2 * n + 1);
    dof = dof % (n * n + 2 * n + 1);
    int j = dof / (n + 1);
    int i = dof % (n + 1);
    pnt[0] = i * (xb - xa) / n;
    pnt[1] = j * (yb - ya) / n;
    pnt[2] = k * (zb - za) / n;
    return pnt;
}

/// 从 j 行 i 列 k 层，每一行 5n 个 P1 剖分单元中映射 (i,j,k) 单元的第 m 个顶点。
AFEPack::Point<3> P1_ele3vtx(int n, int k, int j, int i, int m)
{
    AFEPack::Point<3> pnt;
    int dof = P1_ele3dof(n, k, j, i, m);
    pnt = Dof_to_vtx(n, dof);
    return pnt;
}

/// 储存所有边界自由度编号。
int Store_BndDof(int n)
{
    for (int i = 0; i < (n * n + 2 * n + 1); i++)
    {
	BndDof.push_back(i);
	BndDof.push_back((n + 1) * (n + 1) * (n + 1) - i - 1);
    }
    if (n > 1)
    {
	for (int k = 1; k < n; k++)
	{
	    for (int i = 0; i < n + 1; i++)
	    {
		BndDof.push_back(i + k * (n + 1) * (n + 1));
		BndDof.push_back((k + 1) * (n + 1) * (n + 1) - 1 - i);
	    }
	    for (int i = 1; i < n; i++)
	    {
		BndDof.push_back(k * (n + 1) * (n + 1) + i * (n + 1));
		BndDof.push_back(k * (n + 1) * (n + 1) + i * (n + 1) + n);
	    }
	}
    }
    return 0;
}

int main(int argc, char* argv[])
{
    /// 从 AFEPack 中读入 3D_P1 模板单元格信息，基函数信息和坐标变换信息。
    TemplateGeometry<3> tetrahedron_template_geometry;
    tetrahedron_template_geometry.readData("tetrahedron.tmp_geo");
    CoordTransform<3, 3> tetrahedron_coord_transform;
    tetrahedron_coord_transform.readData("tetrahedron.crd_trs");
    TemplateDOF<3> tetrahedron_template_dof(tetrahedron_template_geometry);
    /// 一次元
    tetrahedron_template_dof.readData("tetrahedron.1.tmp_dof");
    BasisFunctionAdmin<double, 3, 3> tetrahedron_basis_function(tetrahedron_template_dof);
    tetrahedron_basis_function.readData("tetrahedron.1.bas_fun");
    TemplateElement<double, 3, 3> template_element;
    /// 模板单元初始化
    template_element.reinit(tetrahedron_template_geometry,
			    tetrahedron_template_dof,
			    tetrahedron_coord_transform,
			    tetrahedron_basis_function);
    double volume = template_element.volume();
    /// 取 4 次代数精度。
    const QuadratureInfo<3>& quad_info = template_element.findQuadratureInfo(4);
    /// 积分点个数
    int n_quadrature_point = quad_info.n_quadraturePoint();
    /// 积分点
    std::vector< AFEPack::Point<3> > q_point = quad_info.quadraturePoint();
    int n_dof = template_element.n_dof();
    /// 基函数个数
    int n_bas = template_element.basisFunction().size();
    /// 观察一下模板单元中的自由度，基函数和基函数在具体积分点的取值情况，可从模板单元中读取到。
    TemplateGeometry<3> &geo = template_element.geometry();
    const std::vector<AFEPack::Point<3> > &lv = geo.vertexArray();
    int n_vtx = lv.size();
    /// 全局点
    std::vector<AFEPack::Point<3> > gv(n_vtx);
    /// x,y,z 方向的剖分段数。
    int n = 20;
    /// 自由度总数
    int dim = (n + 1) * (n + 1) * (n + 1);
    /// 右端项
    Vector<double> rhs(dim);
    /// 记得先储存边界自由度！！！！！！！！
    Store_BndDof(n);
    /// 稀疏矩阵中非零元的个数，不知道最多多少个，但是不会超过 27 个。
    std::vector<unsigned int> NoZeroPerRow(dim, 27);
    /// 建立稀疏矩阵模板。
    SparsityPattern sp_stiff_matrix(dim, NoZeroPerRow);
    /// 填充非零元素对应的行索引和列索引，遍历顺序按照单元的顺序。
    for (int k = 0; k < n; k++)
	for (int j = 0; j < n; j++)
	{
	    for (int i = 0; i < 5 * n ; i++)
	    {
		for (int dof1 = 0; dof1 < n_dof; dof1++)
		    for (int dof2 = 0; dof2 < n_dof; dof2++)
		    {
			sp_stiff_matrix.add(P1_ele3dof(n, k, j, i, dof1),
					    P1_ele3dof(n, k, j, i, dof2));
		    }
	    }
	}
    
    /// 稀疏矩阵模板生成。
    sp_stiff_matrix.compress();
    /// 刚度矩阵初始化。
    SparseMatrix<double> stiff_mat(sp_stiff_matrix);
    /// 生成节点，单元，刚度矩阵和右端项。
    for (int k = 0; k < n; k++)
	for (int j = 0; j < n; j++)
	    for (int i = 0; i < 5 * n; i++)
	    {
		/// 实际单元顶点坐标。。
		for (int m = 0; m < n_vtx; m++)
		{
		    gv[m] = P1_ele3vtx(n, k, j, i, m);
		}
		/// 积分
		for (int l = 0; l < n_quadrature_point; l++)
		{
		    /// 积分点全局坐标。
		    auto point = tetrahedron_coord_transform.local_to_global(q_point, lv, gv);
		    /// 积分点权重，Jacobi变换系数。
		    double Jxy = quad_info.weight(l) * tetrahedron_coord_transform.local_to_global_jacobian(q_point[l], lv, gv) * volume;

		    for (int base1 = 0; base1 < n_dof; base1++)
		    {
			for (int base2 = 0; base2 < n_dof; base2++)
			    stiff_mat.add(P1_ele3dof(n, k, j, i, base1),
					  P1_ele3dof(n, k, j, i, base2),
					  Jxy * innerProduct(tetrahedron_basis_function[base1].gradient(point[l], gv), tetrahedron_basis_function[base2].gradient(point[l], gv)));
			/// 右端项
			rhs(P1_ele3dof(n, k, j, i, base1)) += Jxy * f(point[l]) * tetrahedron_basis_function[base1].value(point[l], gv);
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
	AFEPack::Point<3> bnd_point;
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

    /// 计算 L2 误差。
    double error = 0;
    for (int dof = 0; dof < dim; dof++)
    {
	AFEPack::Point<3> pnt;
	pnt = Dof_to_vtx(n, dof);
	double d = (u(pnt) - solution[dof]);
	error += d * d;
    }
    error = std::sqrt(error);
    std::cerr << "\nL2 error = " << error << ", tol = " << tol << std::endl;
    
    return 0;
}
	




	    
	    
	
	    
		
	    
	    
	
	    
	    
	    
	    
