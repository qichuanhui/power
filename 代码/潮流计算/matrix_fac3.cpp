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
#include <pthread.h>
// 包含CKTSO库头文件
#include "cktso.h"


// 增强稳定性的矩阵分解
ICktSo pf_fac_st(double *Bp_x, unsigned *Bp_i, unsigned *Bp_p, int node_num) {
    struct timeval start, finish;
    gettimeofday(&start, 0);
    // 1. 创建CKTSO求解器实例
    ICktSo inst = nullptr;
    int* iparm = nullptr;
    const long long* oparm = nullptr;
    int ret = CKTSO_CreateSolver(&inst, &iparm, &oparm);
    if (ret != 0 || inst == nullptr) {
        std::cerr << "Failed to create CKTSO solver instance. Error: " << ret << std::endl;
        return nullptr;
    }
    gettimeofday(&finish, 0);
    double cost_time;
    cost_time = 1000000 * (finish.tv_sec - start.tv_sec) + finish.tv_usec - start.tv_usec;
    printf("createsolver cost the time is : %lf s.\n", cost_time / 1000000.0);
    // 设置输入参数（使用默认值或根据需要调整）
    // 例如：iparm[1] = 1000; // 设置主元容差
    
    // 2. 分析矩阵
    // 将unsigned数组转换为int数组（CKTSO要求int类型）
    gettimeofday(&start, 0);
    const int* col_ptr = reinterpret_cast<const int*>(Bp_p);
    const int* row_idx = reinterpret_cast<const int*>(Bp_i);
    ret = CKTSO_Analyze(inst, false, node_num, 
                         col_ptr, row_idx, 
                         Bp_x, 
                         0); 
    
    if (ret != 0) {
        std::cerr << "CKTSO_Analyze failed. Error: " << ret << std::endl;
        CKTSO_DestroySolver(inst);
        return nullptr;
    }
     gettimeofday(&finish, 0);
     cost_time = 1000000 * (finish.tv_sec - start.tv_sec) + finish.tv_usec - start.tv_usec;
    printf("analyze cost the time is : %lf s.\n", cost_time / 1000000.0);
    // 3. 分解矩阵
    gettimeofday(&start, 0);
    ret = CKTSO_Factorize(inst, Bp_x, true); // 使用完整分解（非快速模式）
    if (ret != 0) {
        std::cerr << "CKTSO_Factorize failed. Error: " << ret << std::endl;
        CKTSO_DestroySolver(inst);
        return nullptr;
    }
    gettimeofday(&finish, 0);
     cost_time = 1000000 * (finish.tv_sec - start.tv_sec) + finish.tv_usec - start.tv_usec;
    printf("factorize cost the time is : %lf s.\n", cost_time / 1000000.0);
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
 
    struct Fac_Para_Mt fac_para_mt[2];
	fac_para_mt[0].Bp_x=Bp_x;
	fac_para_mt[0].Bp_i=Bp_i;
	fac_para_mt[0].Bp_p=Bp_p;
	fac_para_mt[0].node_num=node_num;
	fac_para_mt[1].Bp_x=Bpp_x;
	fac_para_mt[1].Bp_i=Bpp_i;
	fac_para_mt[1].Bp_p=Bpp_p;
	fac_para_mt[1].node_num=node_num;

	pthread_t fac_p_arr[2];
	pthread_create(&fac_p_arr[0],NULL,pf_fac_mt,&fac_para_mt[0]);
	pthread_create(&fac_p_arr[1],NULL,pf_fac_mt,&fac_para_mt[1]);
	pthread_join(fac_p_arr[0],NULL);
	pthread_join(fac_p_arr[1],NULL);

	*solver_p=fac_para_mt[0].solver;
	*solver_pp=fac_para_mt[1].solver;
}

// 求解器销毁函数
void pf_destroy_solver(ICktSo solver) {
    if (solver) {
        CKTSO_DestroySolver(solver);
    }
}