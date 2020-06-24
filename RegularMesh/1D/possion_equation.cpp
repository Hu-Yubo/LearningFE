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
    BndDOF.push_back(n + 1);
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
    
    
}

    
    
