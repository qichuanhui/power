#include "pf.h"
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <cmath>
#include <cstring>
#include <algorithm>
#include <iostream>
#include <sys/time.h>
#include <pthread.h>
#include <cstdarg>

#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <Eigen/SparseLU>
#include <vector>
#include <thread>

// 增强稳定性的矩阵分解
EigenSolver* pf_fac_st(double *Bp_x, unsigned *Bp_i, unsigned *Bp_p, int node_num) {
    try {
        // 1. 创建稀疏矩阵 (CSC格式)
        Eigen::SparseMatrix<double> mat(node_num, node_num);
        std::vector<Eigen::Triplet<double>> triplets;
        triplets.reserve(Bp_p[node_num]);
        
        // 2. 填充矩阵并检查对称性
        for (int col = 0; col < node_num; ++col) {
            for (int idx = Bp_p[col]; idx < Bp_p[col + 1]; ++idx) {
                int row = Bp_i[idx];
                double val = Bp_x[idx];
                // 检查NaN/Inf
                if (std::isnan(val) || std::isinf(val)) {
                    std::cerr << "Error: Invalid matrix value at (" 
                              << row << "," << col << "): " << val << std::endl;
                    return nullptr;
                }
                triplets.emplace_back(row, col, val);
            }
        }
        
        // 3. 设置矩阵
        mat.setFromTriplets(triplets.begin(), triplets.end());
        mat.makeCompressed();
        
        // 4. 选择专用分解器
        EigenSolver* solver = nullptr;
        
            // 非对称矩阵使用带主元的LU
            Eigen::SparseLU<Eigen::SparseMatrix<double>>* lu = 
                new Eigen::SparseLU<Eigen::SparseMatrix<double>>();
            
            // 设置稳定性参数
            lu->setPivotThreshold(0.01);  // 主元阈值
            
            lu->compute(mat);
            
            if (lu->info() != Eigen::Success) {
                std::cerr << "LU decomposition failed: ";
                if (lu->info() == Eigen::NumericalIssue)
                    std::cerr << "Numerical issue detected";
                else
                    std::cerr << "Unknown error";
                std::cerr << std::endl;
                
                delete lu;
                return nullptr;
            }
            
            solver = reinterpret_cast<EigenSolver*>(lu);
        
        return solver;
    } catch (const std::exception& e) {
        std::cerr << "Exception in factorization: " << e.what() << std::endl;
        return nullptr;
    }
}

// 多线程分解结构体
struct Fac_Para_Mt {
    double *Bp_x;
    unsigned *Bp_i;
    unsigned *Bp_p;
    int node_num;
    EigenSolver* eigenSolver;
};

// 线程函数
void* pf_fac_mt(void* arg) {
    Fac_Para_Mt* para = static_cast<Fac_Para_Mt*>(arg);
    para->eigenSolver = pf_fac_st(para->Bp_x, para->Bp_i, para->Bp_p, para->node_num);
    return nullptr;
}

// 顺序执行两个矩阵分解
void pf_fac_seq(double* Bp_x, unsigned* Bp_i, unsigned* Bp_p,
               double* Bpp_x, unsigned* Bpp_i, unsigned* Bpp_p, int node_num,
               EigenSolver** EigenSolver_p, EigenSolver** EigenSolver_pp) {
    *EigenSolver_p = pf_fac_st(Bp_x, Bp_i, Bp_p, node_num);
    *EigenSolver_pp = pf_fac_st(Bpp_x, Bpp_i, Bpp_p, node_num);
}

// 并行执行两个矩阵分解
void pf_fac_par(double* Bp_x, unsigned* Bp_i, unsigned* Bp_p,
               double* Bpp_x, unsigned* Bpp_i, unsigned* Bpp_p, int node_num,
               EigenSolver** EigenSolver_p, EigenSolver** EigenSolver_pp) {
    Fac_Para_Mt fac_para_mt[2];
    
    // 配置第一个矩阵参数
    fac_para_mt[0] = {Bp_x, Bp_i, Bp_p, node_num, nullptr};
    
    // 配置第二个矩阵参数
    fac_para_mt[1] = {Bpp_x, Bpp_i, Bpp_p, node_num, nullptr};

    // 创建并运行线程
    std::thread threads[2];
    threads[0] = std::thread(pf_fac_mt, &fac_para_mt[0]);
    threads[1] = std::thread(pf_fac_mt, &fac_para_mt[1]);
    
    // 等待线程完成
    threads[0].join();
    threads[1].join();

    // 返回结果
    *EigenSolver_p = fac_para_mt[0].eigenSolver;
    *EigenSolver_pp = fac_para_mt[1].eigenSolver;
}