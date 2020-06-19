/**
 * @file possion_equation_Q1.cpp
 * @author 
 * @date 
 * 
 * @brief 使用AFEpack里的函数，在结构化网格下实现二维一次元Possion方程计算
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

/// 从 j 行 i 列，每一列 n 个网格的 Q1 剖分中映射 (i,j) 单元的第 k 个自由度编号
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
    case 2:
	idx = (j + 1) * (n + 1) + i + 1;
    case 3:
	idx = (j + 1) * (n + 1) + i;
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
    rectangle_template_geometry.readData("rectangle.tmp.geo");
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
    int n_bas = rectangle_basis_function.size();

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

    Vector<double> rhs(dim);
}

    


    
