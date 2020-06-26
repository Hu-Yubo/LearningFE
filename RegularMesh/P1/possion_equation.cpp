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

double u(const double* p)
{
    return sin(PI * p[0] * sin(PI * p[1]));
}

double f(const double* p)
{
    return 2 * PI * PI * u(p);
}

/* 如何划分单元格？                     
 * 首先将实际计算区域划分为规则的矩形网格，对于每个矩形网格，
 * 将其第0个和第2个顶点连线，划分为两个三角形，从左至右排序。
 * 3 ———— 2 
 * |     /|    如图左所示，0-2-3, 0-1-2 构成两个 P1 单元片，
 * |    / |    且 0-2-3 编号在前。
 * |   /  |
 * |  /   |
 * | /    |
 * |/     |
 * 0 ———— 1
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
	    
