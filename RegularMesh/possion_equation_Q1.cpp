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
    
