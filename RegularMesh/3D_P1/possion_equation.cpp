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
    return 2 * PI * PI * u(p);
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
    int n = 2;
    Store_BndDof(n);
    for (int i = 0; i < BndDof.size(); i++)
    {
	std::cout << BndDof[i] << " ";
    }
    return 0;
}
	




	    
	    
	
	    
		
	    
	    
	
	    
	    
	    
	    
