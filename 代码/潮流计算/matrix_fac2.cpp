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
#include <vector>
#include <thread>

// 包含CKTSO库头文件
#include "cktso.h"


// 增强稳定性的矩阵分解
ICktSo pf_fac_st(double *Bp_x, unsigned *Bp_i, unsigned *Bp_p, int node_num) {
    // 1. 创建CKTSO求解器实例
    ICktSo inst = nullptr;
    int* iparm = nullptr;
    const long long* oparm = nullptr;
    int ret = CKTSO_CreateSolver(&inst, &iparm, &oparm);
    if (ret != 0 || inst == nullptr) {
        std::cerr << "Failed to create CKTSO solver instance. Error: " << ret << std::endl;
        return nullptr;
    }
    
    // 设置输入参数（使用默认值或根据需要调整）
    // 例如：iparm[1] = 1000; // 设置主元容差
    
    // 2. 分析矩阵
    // 将unsigned数组转换为int数组（CKTSO要求int类型）
    std::vector<int> ap_vec(Bp_p, Bp_p + node_num + 1);
    std::vector<int> ai_vec(Bp_i, Bp_i + ap_vec[node_num]);
    
    ret = CKTSO_Analyze(inst, false, node_num, 
                         ap_vec.data(), ai_vec.data(), 
                         Bp_x, 
                         16); 
    
    if (ret != 0) {
        std::cerr << "CKTSO_Analyze failed. Error: " << ret << std::endl;
        CKTSO_DestroySolver(inst);
        return nullptr;
    }
    
    // 3. 分解矩阵
    ret = CKTSO_Factorize(inst, Bp_x, false); // 使用完整分解（非快速模式）
    if (ret != 0) {
        std::cerr << "CKTSO_Factorize failed. Error: " << ret << std::endl;
        CKTSO_DestroySolver(inst);
        return nullptr;
    }
    
    // 4. 返回求解器实例
    return inst;
}

// 多线程分解结构体
struct Fac_Para_Mt {
    double *Bp_x;
    unsigned *Bp_i;
    unsigned *Bp_p;
    int node_num;
    ICktSo solver; // 修改为CKTSOSolver类型
};

// 线程函数
void* pf_fac_mt(void* arg) {
    Fac_Para_Mt* para = static_cast<Fac_Para_Mt*>(arg);
    para->solver = pf_fac_st(para->Bp_x, para->Bp_i, para->Bp_p, para->node_num);
    return nullptr;
}

// 顺序执行两个矩阵分解
void pf_fac_seq(double* Bp_x, unsigned* Bp_i, unsigned* Bp_p,
               double* Bpp_x, unsigned* Bpp_i, unsigned* Bpp_p, int node_num,
               ICktSo* solver_p, ICktSo* solver_pp) {
    *solver_p = pf_fac_st(Bp_x, Bp_i, Bp_p, node_num);
    *solver_pp = pf_fac_st(Bpp_x, Bpp_i, Bpp_p, node_num);
}

// 并行执行两个矩阵分解
void pf_fac_par(double* Bp_x, unsigned* Bp_i, unsigned* Bp_p,
               double* Bpp_x, unsigned* Bpp_i, unsigned* Bpp_p, int node_num,
               ICktSo* solver_p, ICktSo* solver_pp) {
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
    *solver_p = fac_para_mt[0].solver;
    *solver_pp = fac_para_mt[1].solver;
}

// 求解器销毁函数
void pf_destroy_solver(ICktSo solver) {
    if (solver) {
        CKTSO_DestroySolver(solver);
    }
}