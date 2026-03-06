#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <sys/time.h>
#include "pf2.h"
#include <math.h>
#include <unistd.h>
#include <vector>
#include <cstring> // for memset
// XRT 头文件
#include <xrt/xrt_bo.h>
#include <xrt/xrt_device.h>
#include <xrt/xrt_kernel.h>
#include <experimental/xrt_xclbin.h>
#include <experimental/xrt_graph.h>

#define NUM_KERNEL 4
// HBM Banks requirements
#define MAX_DDR_BANKCOUNT 4
const int bank[MAX_DDR_BANKCOUNT] = {0, 1, 2, 3}; // 使用 bank 索引

double* host_node_P[NUM_KERNEL];
double* host_node_Q[NUM_KERNEL];
double* host_node_Vm[NUM_KERNEL];
double* host_node_Va[NUM_KERNEL];
int * host_node_type[NUM_KERNEL];
int * host_edge_off[NUM_KERNEL];
int* host_edge_src[NUM_KERNEL];
double* host_edge_G[NUM_KERNEL];
double* host_edge_B[NUM_KERNEL];
double* host_deltaP[NUM_KERNEL];
double* host_deltaQ[NUM_KERNEL];
double* host_max_delta;

int host_node_begin[NUM_KERNEL], host_node_end[NUM_KERNEL];

bool init_flag;
bool Vm_flag;
bool Va_flag;
bool deltaP_flag;
bool deltaQ_flag;
bool compute_flag;

#define GROUP_SIZE 1
#define MAX_ITER 10
int debug = 0;
double total_time[MAX_ITER * 10];
int total_time_count = 0;


void program_exit(std::string casename, double precision, int node_num, int edge_num, int iter,
                 double Ybus_time, double fac_time, double rhs_time, double back_for_time, double convert_time,
                 struct Ret_Ybus ret_Ybus, struct Graph_PF_IN graph_in, double *total_time, int total_time_count);

int main(int argc, char *argv[]) {
    struct timeval start, finish;
    double cost_time;
    double Ybus_time = 0, fac_time = 0, rhs_time = 0, back_for_time = 0;
    double convert_time = 0;

    if (argc != 4) {
        std::cout << "Usage: " << argv[0] << " <xclbin> <casename> <precision>" << std::endl;
        return EXIT_FAILURE;
    }
    std::string binaryFile = argv[1];
    std::string casename = argv[2];
    double precision = atof(argv[3]);
    printf("dataset:%s\tprecision:%f\n", casename.c_str(), precision);

    // ========== XRT 初始化 ==========
    auto device = xrt::device(1); // 使用第二个设备
    auto uuid = device.load_xclbin(binaryFile);
    auto krnl_pf = xrt::kernel(device, uuid, "krnl_pf");

    // ========== 加载数据 ==========
    struct Graph_PF_IN graph_in = pf_load(casename);
    gettimeofday(&start, 0);
    struct Ret_Ybus ret_Ybus = pf_Ybus_par(graph_in);
    gettimeofday(&finish, 0);
    cost_time = 1000000 * (finish.tv_sec - start.tv_sec) + finish.tv_usec - start.tv_usec;
    printf("node_num=%d\tedge_num=%d\n", ret_Ybus.graph.node_num, ret_Ybus.graph.edge_num);
    printf("initial and create matrix: cost the time is : %lf s.\n", cost_time / (double)1000000);
    Ybus_time = cost_time /(double) 1000000;
    total_time[total_time_count++] = cost_time /(double)1000000*1000;

    gettimeofday(&start, 0);
    ICktSo nicslu_p = nullptr,nicslu_pp = nullptr;
    pf_fac_par(ret_Ybus.Bp_x, ret_Ybus.Bp_i, ret_Ybus.Bp_p,
               ret_Ybus.Bpp_x, ret_Ybus.Bpp_i, ret_Ybus.Bpp_p,
               ret_Ybus.graph.node_num, &nicslu_p, &nicslu_pp);
    gettimeofday(&finish, 0);
    cost_time = 1000000 * (finish.tv_sec - start.tv_sec) + finish.tv_usec - start.tv_usec;
    printf("matrix factorization: cost the time is : %lf s.\n", cost_time / 1000000.0);
    fac_time = cost_time / 1000000.0;
    total_time[total_time_count++] = cost_time / 1000.0;

    struct Node_PF *node = ret_Ybus.graph.node;
    struct Edge_PF *edge = ret_Ybus.graph.edge;
    unsigned *off = ret_Ybus.graph.off;
    int node_num = ret_Ybus.graph.node_num;
    int edge_num = ret_Ybus.graph.edge_num;

    // 计算节点分区
    {
        int block_size = node_num / (8 * NUM_KERNEL) * 8;
        host_node_begin[0] = 0;
        host_node_end[0] = block_size;
        host_node_begin[1] = block_size;
        host_node_end[1] = 2 * block_size;
        host_node_begin[2] = 2 * block_size;
        host_node_end[2] = 3 * block_size;
        host_node_begin[3] = 3 * block_size;
        host_node_end[3] = node_num;
    }

    // ========== 主机内存分配 ==========
    gettimeofday(&start, 0);
    double *deltaP = (double *)malloc(node_num * sizeof(double));
    double *deltaQ = (double *)malloc(node_num * sizeof(double));
    double *deltaP_debug = (double *)malloc(node_num * sizeof(double));
    double *deltaQ_debug = (double *)malloc(node_num * sizeof(double));

    int posix_ret;
    for (int i = 0; i < NUM_KERNEL; i++) {
        int l_node_num = host_node_end[i] - host_node_begin[i];
        size_t node_PQ_len = ((l_node_num * 8) % 64) ? ((l_node_num * 8) / 64 * 64 + 64) : (l_node_num * 8);
        posix_ret = posix_memalign((void **)&host_node_P[i], 4096, node_PQ_len);
        memset(&host_node_P[i][l_node_num], 0, node_PQ_len - l_node_num * 8);
        posix_ret = posix_memalign((void **)&host_node_Q[i], 4096, node_PQ_len);
        memset(&host_node_Q[i][l_node_num], 0, node_PQ_len - l_node_num * 8);

        size_t node_V_len = ((l_node_num * 8) % 64) ? ((l_node_num * 8) / 64 * 64 + 64) : (l_node_num * 8);
        posix_ret = posix_memalign((void **)&host_node_Vm[i], 4096, node_V_len);
        memset(&host_node_Vm[i][l_node_num], 0x0c, node_V_len - l_node_num * 8);
        posix_ret = posix_memalign((void **)&host_node_Va[i], 4096, node_V_len);
        memset(&host_node_Va[i][l_node_num], 0, node_V_len - l_node_num * 8);

        size_t node_type_len = ((l_node_num * 4) % 64) ? ((l_node_num * 4) / 64 * 64 + 64) : (l_node_num * 4);
        posix_ret = posix_memalign((void **)&host_node_type[i], 4096, node_type_len);
        memset(&host_node_type[i][l_node_num], 0, node_type_len - l_node_num * 4);

        size_t edge_off_len = (((l_node_num + 1) * 4) % 64) ? (((l_node_num + 1) * 4) / 64 * 64 + 64) : ((l_node_num + 1) * 4);
        posix_ret = posix_memalign((void **)&host_edge_off[i], 4096, edge_off_len);
        memset(&host_edge_off[i][l_node_num + 1], 0, edge_off_len - (l_node_num + 1) * 4);

        int l_edge_num = off[host_node_end[i]] - off[host_node_begin[i]];
        size_t edge_src_len = ((l_edge_num * 4) % 64) ? ((l_edge_num * 4) / 64 * 64 + 64) : (l_edge_num * 4);
        posix_ret = posix_memalign((void **)&host_edge_src[i], 4096, edge_src_len);
        memset(&host_edge_src[i][l_edge_num], 0, edge_src_len - l_edge_num * 4);

        size_t edge_GB_len = ((l_edge_num * 8) % 64) ? ((l_edge_num * 8) / 64 * 64 + 64) : (l_edge_num * 8);
        posix_ret = posix_memalign((void **)&host_edge_G[i], 4096, edge_GB_len);
        memset(&host_edge_G[i][l_edge_num], 0, edge_GB_len - l_edge_num * 8);
        posix_ret = posix_memalign((void **)&host_edge_B[i], 4096, edge_GB_len);
        memset(&host_edge_B[i][l_edge_num], 0, edge_GB_len - l_edge_num * 8);

        size_t node_deltaPQ_len = ((l_node_num * 8) % 64) ? ((l_node_num * 8) / 64 * 64 + 64) : (l_node_num * 8);
        posix_ret = posix_memalign((void **)&host_deltaP[i], 4096, node_deltaPQ_len);
        memset(&host_deltaP[i][l_node_num], 0, node_deltaPQ_len - l_node_num * 8);
        posix_ret = posix_memalign((void **)&host_deltaQ[i], 4096, node_deltaPQ_len);
        memset(&host_deltaQ[i][l_node_num], 0, node_deltaPQ_len - l_node_num * 8);
    }
    posix_ret = posix_memalign((void **)&host_max_delta, 4096, 64);
    memset(host_max_delta, 0, 64);

    // 初始化主机数据
    for (int i = 0; i < NUM_KERNEL; i++) {
        int l_node_begin = host_node_begin[i];
        int l_node_num = host_node_end[i] - host_node_begin[i];
        for (int j = 0; j < l_node_num; j++) {
            host_node_P[i][j] = node[l_node_begin + j].P;
            host_node_Q[i][j] = node[l_node_begin + j].Q;
            host_node_Vm[i][j] = node[l_node_begin + j].Vm;
            host_node_Va[i][j] = node[l_node_begin + j].Va;
            host_node_type[i][j] = node[l_node_begin + j].type;
            host_edge_off[i][j] = off[l_node_begin + j];
        }
        host_edge_off[i][l_node_num] = off[l_node_begin + l_node_num];
    }

    for (int i = 0; i < NUM_KERNEL; i++) {
        int l_edge_begin = off[host_node_begin[i]];
        int l_edge_num = off[host_node_end[i]] - off[host_node_begin[i]];
        for (int j = 0; j < l_edge_num; j++) {
            host_edge_src[i][j] = edge[l_edge_begin + j].src;
            host_edge_G[i][j] = edge[l_edge_begin + j].G;
            host_edge_B[i][j] = edge[l_edge_begin + j].B;
        }
    }
    gettimeofday(&finish, 0);
    cost_time = 1000000 * (finish.tv_sec - start.tv_sec) + finish.tv_usec - start.tv_usec;
    printf("host memory initialize: cost the time is : %lf s.\n", cost_time / 1000000.0);

    // ========== XRT 缓冲区创建 ==========
    gettimeofday(&start, 0);
    std::vector<xrt::bo> buffer_node_P(NUM_KERNEL);
    std::vector<xrt::bo> buffer_node_Q(NUM_KERNEL);
    std::vector<xrt::bo> buffer_node_Vm(NUM_KERNEL);
    std::vector<xrt::bo> buffer_node_Va(NUM_KERNEL);
    std::vector<xrt::bo> buffer_node_type(NUM_KERNEL);
    std::vector<xrt::bo> buffer_edge_off(NUM_KERNEL);
    std::vector<xrt::bo> buffer_edge_src(NUM_KERNEL);
    std::vector<xrt::bo> buffer_edge_G(NUM_KERNEL);
    std::vector<xrt::bo> buffer_edge_B(NUM_KERNEL);
    std::vector<xrt::bo> buffer_deltaP(NUM_KERNEL);
    std::vector<xrt::bo> buffer_deltaQ(NUM_KERNEL);
    xrt::bo buffer_max_delta;

    for (int i = 0; i < NUM_KERNEL; i++) {
        int l_node_num = host_node_end[i] - host_node_begin[i];
        size_t node_PQ_len = ((l_node_num * 8) % 64) ? ((l_node_num * 8) / 64 * 64 + 64) : (l_node_num * 8);
        
        // 创建缓冲区并指定内存bank
        buffer_node_P[i] = xrt::bo(device, node_PQ_len, krnl_pf.group_id(bank[i]));
        buffer_node_Q[i] = xrt::bo(device, node_PQ_len, krnl_pf.group_id(bank[i]));
        
        size_t node_V_len = ((l_node_num * 8) % 64) ? ((l_node_num * 8) / 64 * 64 + 64) : (l_node_num * 8);
        buffer_node_Vm[i] = xrt::bo(device, node_V_len, krnl_pf.group_id(bank[i]));
        buffer_node_Va[i] = xrt::bo(device, node_V_len, krnl_pf.group_id(bank[i]));
        
        size_t node_type_len = ((l_node_num * 4) % 64) ? ((l_node_num * 4) / 64 * 64 + 64) : (l_node_num * 4);
        buffer_node_type[i] = xrt::bo(device, node_type_len, krnl_pf.group_id(bank[i]));
        
        size_t edge_off_len = (((l_node_num + 1) * 4) % 64) ? (((l_node_num + 1) * 4) / 64 * 64 + 64) : ((l_node_num + 1) * 4);
        buffer_edge_off[i] = xrt::bo(device, edge_off_len, krnl_pf.group_id(bank[i]));
        
        int l_edge_num = off[host_node_end[i]] - off[host_node_begin[i]];
        size_t edge_src_len = ((l_edge_num * 4) % 64) ? ((l_edge_num * 4) / 64 * 64 + 64) : (l_edge_num * 4);
        buffer_edge_src[i] = xrt::bo(device, edge_src_len, krnl_pf.group_id(bank[i]));
        
        size_t edge_GB_len = ((l_edge_num * 8) % 64) ? ((l_edge_num * 8) / 64 * 64 + 64) : (l_edge_num * 8);
        buffer_edge_G[i] = xrt::bo(device, edge_GB_len, krnl_pf.group_id(bank[i]));
        buffer_edge_B[i] = xrt::bo(device, edge_GB_len, krnl_pf.group_id(bank[i]));
        
        size_t node_deltaPQ_len = ((l_node_num * 8) % 64) ? ((l_node_num * 8) / 64 * 64 + 64) : (l_node_num * 8);
        buffer_deltaP[i] = xrt::bo(device, node_deltaPQ_len, krnl_pf.group_id(bank[i]));
        buffer_deltaQ[i] = xrt::bo(device, node_deltaPQ_len, krnl_pf.group_id(bank[i]));
    }
    
    buffer_max_delta = xrt::bo(device, 64, krnl_pf.group_id(bank[0]));
    gettimeofday(&finish, 0);
    cost_time = 1000000 * (finish.tv_sec - start.tv_sec) + finish.tv_usec - start.tv_usec;
    printf("device memory allocate: cost the time is : %lf s.\n", cost_time / 1000000.0);

    // ========== 主计算循环 ==========
int iter;
double maxDeltaP = 12345678;
double maxDeltaQ = 12345678;

for (iter = 0; iter < MAX_ITER; ++iter) {
    maxDeltaP = 0;
    maxDeltaQ = 0;
    
    // ========== 更新主机缓冲区中的Vm和Va ==========
    for (int i = 0; i < NUM_KERNEL; i++) {
        int l_node_begin = host_node_begin[i];
        int l_node_num = host_node_end[i] - host_node_begin[i];
        for (int j = 0; j < l_node_num; j++) {
            host_node_Vm[i][j] = node[l_node_begin + j].Vm;
            host_node_Va[i][j] = node[l_node_begin + j].Va;
        }
    }

    // ========== 数据传输到设备 ==========
    gettimeofday(&start, 0);
    
    // 1. 首次迭代需要传输所有数据
    if (iter == 0) {
        for (int i = 0; i < NUM_KERNEL; i++) {
            // 映射缓冲区并复制数据
            auto p_ptr = buffer_node_P[i].map<double*>();
            std::memcpy(p_ptr, host_node_P[i], buffer_node_P[i].size());
            
            auto q_ptr = buffer_node_Q[i].map<double*>();
            std::memcpy(q_ptr, host_node_Q[i], buffer_node_Q[i].size());
            
            auto type_ptr = buffer_node_type[i].map<int*>();
            std::memcpy(type_ptr, host_node_type[i], buffer_node_type[i].size());
            
            auto off_ptr = buffer_edge_off[i].map<int*>();
            std::memcpy(off_ptr, host_edge_off[i], buffer_edge_off[i].size());
            
            auto src_ptr = buffer_edge_src[i].map<int*>();
            std::memcpy(src_ptr, host_edge_src[i], buffer_edge_src[i].size());
            
            auto g_ptr = buffer_edge_G[i].map<double*>();
            std::memcpy(g_ptr, host_edge_G[i], buffer_edge_G[i].size());
            
            auto b_ptr = buffer_edge_B[i].map<double*>();
            std::memcpy(b_ptr, host_edge_B[i], buffer_edge_B[i].size());
            
            // 同步到设备
            buffer_node_P[i].sync(XCL_BO_SYNC_BO_TO_DEVICE);
            buffer_node_Q[i].sync(XCL_BO_SYNC_BO_TO_DEVICE);
            buffer_node_type[i].sync(XCL_BO_SYNC_BO_TO_DEVICE);
            buffer_edge_off[i].sync(XCL_BO_SYNC_BO_TO_DEVICE);
            buffer_edge_src[i].sync(XCL_BO_SYNC_BO_TO_DEVICE);
            buffer_edge_G[i].sync(XCL_BO_SYNC_BO_TO_DEVICE);
            buffer_edge_B[i].sync(XCL_BO_SYNC_BO_TO_DEVICE);
        }
    }
    
    // 2. 每次迭代都需要传输更新的Vm和Va
    for (int i = 0; i < NUM_KERNEL; i++) {
        auto vm_ptr = buffer_node_Vm[i].map<double*>();
        std::memcpy(vm_ptr, host_node_Vm[i], buffer_node_Vm[i].size());
        buffer_node_Vm[i].sync(XCL_BO_SYNC_BO_TO_DEVICE);
        
        auto va_ptr = buffer_node_Va[i].map<double*>();
        std::memcpy(va_ptr, host_node_Va[i], buffer_node_Va[i].size());
        buffer_node_Va[i].sync(XCL_BO_SYNC_BO_TO_DEVICE);
    }
    
    gettimeofday(&finish, 0);
    cost_time = 1000000 * (finish.tv_sec - start.tv_sec) + finish.tv_usec - start.tv_usec;
    printf("iter %d: write device: cost the time is : %lf s.\n", iter, cost_time / 1000000.0);
    convert_time += cost_time / 1000000.0;
    total_time[total_time_count++] = cost_time / 1000.0;
    
    // ========== 内核执行 ==========
    gettimeofday(&start, 0);
    
    // 设置内核标志（保持原逻辑）
    if (iter == 0) {
        init_flag = true;
        Vm_flag = true;
        Va_flag = true;
        deltaP_flag = true;
        deltaQ_flag = false;
        compute_flag = true;
    } else {
        init_flag = false;
        Vm_flag = true;
        Va_flag = false;
        deltaP_flag = true;
        deltaQ_flag = false;
        compute_flag = true;
    }
    
    // 创建内核执行对象并传递所有参数
auto run = krnl_pf(
                buffer_node_P[0], buffer_node_P[1], buffer_node_P[2], buffer_node_P[3],
                buffer_node_Q[0], buffer_node_Q[1], buffer_node_Q[2], buffer_node_Q[3],
                buffer_node_Vm[0], buffer_node_Vm[1], buffer_node_Vm[2], buffer_node_Vm[3],
                buffer_node_Va[0], buffer_node_Va[1], buffer_node_Va[2], buffer_node_Va[3],
                buffer_node_type[0], buffer_node_type[1], buffer_node_type[2], buffer_node_type[3],
                buffer_edge_off[0], buffer_edge_off[1], buffer_edge_off[2], buffer_edge_off[3],
                buffer_edge_src[0], buffer_edge_src[1], buffer_edge_src[2], buffer_edge_src[3],
                buffer_edge_G[0], buffer_edge_G[1], buffer_edge_G[2], buffer_edge_G[3],
                buffer_edge_B[0], buffer_edge_B[1], buffer_edge_B[2], buffer_edge_B[3],
                buffer_deltaP[0], buffer_deltaP[1], buffer_deltaP[2], buffer_deltaP[3],
                buffer_deltaQ[0], buffer_deltaQ[1], buffer_deltaQ[2], buffer_deltaQ[3],
                buffer_max_delta,
                host_node_begin[0], host_node_begin[1], host_node_begin[2], host_node_begin[3],
                host_node_end[0], host_node_end[1], host_node_end[2], host_node_end[3],
                node_num,
                init_flag, Vm_flag, Va_flag, deltaP_flag, deltaQ_flag, compute_flag
            );
    
    // 等待内核完成
    run.wait();
    
    gettimeofday(&finish, 0);
    cost_time = 1000000 * (finish.tv_sec - start.tv_sec) + finish.tv_usec - start.tv_usec;
    printf("iter %d: compute rhs: cost the time is : %lf s.\n", iter, cost_time / 1000000.0);
    rhs_time += cost_time / 1000000.0;
    total_time[total_time_count++] = cost_time / 1000.0;
    
    // ========== 读取结果 ==========
    gettimeofday(&start, 0);
    
    // 同步结果到主机
    for (int i = 0; i < NUM_KERNEL; i++) {
        buffer_deltaP[i].sync(XCL_BO_SYNC_BO_FROM_DEVICE);
        buffer_deltaQ[i].sync(XCL_BO_SYNC_BO_FROM_DEVICE);
    }
    buffer_max_delta.sync(XCL_BO_SYNC_BO_FROM_DEVICE);
    
    // 映射结果缓冲区并复制数据
    for (int i = 0; i < NUM_KERNEL; i++) {
        auto deltaP_ptr = buffer_deltaP[i].map<double*>();
        std::memcpy(host_deltaP[i], deltaP_ptr, buffer_deltaP[i].size());
        
        auto deltaQ_ptr = buffer_deltaQ[i].map<double*>();
        std::memcpy(host_deltaQ[i], deltaQ_ptr, buffer_deltaQ[i].size());
    }
    auto max_delta_ptr = buffer_max_delta.map<double*>();
    std::memcpy(host_max_delta, max_delta_ptr, buffer_max_delta.size());
    
    // 更新主机数据
    for (int i = 0; i < NUM_KERNEL; i++) {
        int l_node_begin = host_node_begin[i];
        int l_node_end = host_node_end[i];
        int l_node_num = l_node_end - l_node_begin;
        for (int j = 0; j < l_node_num; j++) {
            deltaP[l_node_begin + j] = host_deltaP[i][j];
            deltaQ[l_node_begin + j] = host_deltaQ[i][j];
        }
    }
    maxDeltaP = host_max_delta[0];
    maxDeltaQ = host_max_delta[1];
    
    gettimeofday(&finish, 0);
    cost_time = 1000000 * (finish.tv_sec - start.tv_sec) + finish.tv_usec - start.tv_usec;
    printf("iter %d: read device: cost the time is : %lf s.\n", iter, cost_time / 1000000.0);
    convert_time += cost_time / 1000000.0;
    total_time[total_time_count++] = cost_time / 1000.0;

    std::cout << "iter " << iter << ": maxDeltaP, " << maxDeltaP << ", maxDeltaQ, " << maxDeltaQ << std::endl;
    // 检查收敛
    if (maxDeltaP < precision && maxDeltaQ < precision) {
        break;
    }

    // ========== 求解并更新Va ==========
    gettimeofday(&start, 0);
    CKTSO_Solve(nicslu_p, deltaP, deltaP, false, false);
    for (int i = 0; i < node_num; ++i) {
        node[i].Va -= deltaP[i];
        // 添加角度规范化
    //    if (node[i].Va >= M_PI) node[i].Va -= 2 * M_PI;
    //    if (node[i].Va < -M_PI) node[i].Va += 2 * M_PI;
    }
    gettimeofday(&finish, 0);
    cost_time = 1000000 * (finish.tv_sec - start.tv_sec) + finish.tv_usec - start.tv_usec;
    printf("back_forward and update Va: cost the time is : %lf s.\n", cost_time / 1000000.0);
    back_for_time += cost_time / 1000000.0;
    total_time[total_time_count++] = cost_time / 1000.0;

    // ========== 更新Va后立即传输到设备 ==========
    gettimeofday(&start, 0);
    for (int i = 0; i < NUM_KERNEL; i++) {
        int l_node_begin = host_node_begin[i];
        int l_node_num = host_node_end[i] - host_node_begin[i];
        for (int j = 0; j < l_node_num; j++) {
            host_node_Va[i][j] = node[l_node_begin + j].Va;
        }
        auto va_ptr = buffer_node_Va[i].map<double*>();
        std::memcpy(va_ptr, host_node_Va[i], buffer_node_Va[i].size());
        buffer_node_Va[i].sync(XCL_BO_SYNC_BO_TO_DEVICE);
    }
    gettimeofday(&finish, 0);
    cost_time = 1000000 * (finish.tv_sec - start.tv_sec) + finish.tv_usec - start.tv_usec;
    printf("iter %d: write device (Va): cost the time is : %lf s.\n", iter, cost_time / 1000000.0);
    convert_time += cost_time / 1000000.0;
    total_time[total_time_count++] = cost_time / 1000.0;

    // ========== 第二次内核执行 ==========
    gettimeofday(&start, 0);
    
    // 设置内核标志（保持原逻辑）
    init_flag = false;
    Vm_flag = false;
    Va_flag = true;
    deltaP_flag = false;
    deltaQ_flag = true;
    compute_flag = true;
    
  // 创建内核执行对象
            auto run2 = krnl_pf(
                buffer_node_P[0], buffer_node_P[1], buffer_node_P[2], buffer_node_P[3],
                buffer_node_Q[0], buffer_node_Q[1], buffer_node_Q[2], buffer_node_Q[3],
                buffer_node_Vm[0], buffer_node_Vm[1], buffer_node_Vm[2], buffer_node_Vm[3],
                buffer_node_Va[0], buffer_node_Va[1], buffer_node_Va[2], buffer_node_Va[3],
                buffer_node_type[0], buffer_node_type[1], buffer_node_type[2], buffer_node_type[3],
                buffer_edge_off[0], buffer_edge_off[1], buffer_edge_off[2], buffer_edge_off[3],
                buffer_edge_src[0], buffer_edge_src[1], buffer_edge_src[2], buffer_edge_src[3],
                buffer_edge_G[0], buffer_edge_G[1], buffer_edge_G[2], buffer_edge_G[3],
                buffer_edge_B[0], buffer_edge_B[1], buffer_edge_B[2], buffer_edge_B[3],
                buffer_deltaP[0], buffer_deltaP[1], buffer_deltaP[2], buffer_deltaP[3],
                buffer_deltaQ[0], buffer_deltaQ[1], buffer_deltaQ[2], buffer_deltaQ[3],
                buffer_max_delta,
                host_node_begin[0], host_node_begin[1], host_node_begin[2], host_node_begin[3],
                host_node_end[0], host_node_end[1], host_node_end[2], host_node_end[3],
                node_num,
                init_flag, Vm_flag, Va_flag, deltaP_flag, deltaQ_flag, compute_flag
            );
    
    // 等待内核完成
    run2.wait();
    
    gettimeofday(&finish, 0);
    cost_time = 1000000 * (finish.tv_sec - start.tv_sec) + finish.tv_usec - start.tv_usec;
    printf("iter %d: compute rhs (Q): cost the time is : %lf s.\n", iter, cost_time / 1000000.0);
    rhs_time += cost_time / 1000000.0;
    total_time[total_time_count++] = cost_time / 1000.0;
    
    // ========== 读取结果 ==========
    gettimeofday(&start, 0);
    
    // 同步结果到主机
    for (int i = 0; i < NUM_KERNEL; i++) {
        buffer_deltaQ[i].sync(XCL_BO_SYNC_BO_FROM_DEVICE);
    }
    buffer_max_delta.sync(XCL_BO_SYNC_BO_FROM_DEVICE);
    
    // 映射结果缓冲区并复制数据
    for (int i = 0; i < NUM_KERNEL; i++) {
        auto deltaQ_ptr = buffer_deltaQ[i].map<double*>();
        std::memcpy(host_deltaQ[i], deltaQ_ptr, buffer_deltaQ[i].size());
    }
    auto max_delta_ptr2 = buffer_max_delta.map<double*>();
    std::memcpy(host_max_delta, max_delta_ptr2, buffer_max_delta.size());
    
    // 更新主机数据
    for (int i = 0; i < NUM_KERNEL; i++) {
        int l_node_begin = host_node_begin[i];
        int l_node_end = host_node_end[i];
        int l_node_num = l_node_end - l_node_begin;
        for (int j = 0; j < l_node_num; j++) {
            deltaQ[l_node_begin + j] = host_deltaQ[i][j];
        }
    }
    maxDeltaP = host_max_delta[0];
    maxDeltaQ = host_max_delta[1];
    gettimeofday(&finish, 0);
    cost_time = 1000000 * (finish.tv_sec - start.tv_sec) + finish.tv_usec - start.tv_usec;
    printf("iter %d: read device (Q): cost the time is : %lf s.\n", iter, cost_time / 1000000.0);
    convert_time += cost_time / 1000000.0;
    total_time[total_time_count++] = cost_time / 1000.0;
    std::cout << "iter " << iter << ": maxDeltaP, " << maxDeltaP << ", maxDeltaQ, " << maxDeltaQ << std::endl;
    // 检查收敛
    if (maxDeltaP < precision && maxDeltaQ < precision) {
        break;
    }
    // ========== 求解并更新Vm ==========
    gettimeofday(&start, 0);
    CKTSO_Solve(nicslu_pp, deltaQ , deltaQ, false, false);
    for (int i = 0; i < node_num; ++i) {
        node[i].Vm -= deltaQ[i];
    }
    gettimeofday(&finish, 0);
    cost_time = 1000000 * (finish.tv_sec - start.tv_sec) + finish.tv_usec - start.tv_usec;
    printf("back_forward and update Vm: cost the time is : %lf s.\n", cost_time / 1000000.0);
    back_for_time += cost_time / 1000000.0;
    total_time[total_time_count++] = cost_time / 1000.0;
}
    // ========== 结果输出和清理 ==========
    printf("----------------information---------------\n");
    printf("dataset                   : %s\n", casename.c_str());
    printf("precision                 : %lf\n", precision);
    printf("nodes                     : %d\n", node_num);
    printf("edges                     : %d\n", edge_num);
    printf("iterations                : %d\n", iter);
    printf("------------------------------------------\n");
    printf("----------------excute time---------------\n");
    printf("initial and create matrix : %lfs\n", Ybus_time/3);
    printf("matrix factorization      : %lfs\n", fac_time/3);
    printf("compute rhs               : %lfs\n", rhs_time/3);
    printf("back and forward          : %lfs\n", back_for_time/3);
    printf("convert                   : %lfs\n", convert_time/3);
    printf("total                     : %lfs\n", (Ybus_time + fac_time + rhs_time + back_for_time + convert_time)/3);
    printf("------------------------------------------\n");

    power_show_csv(ret_Ybus.graph, graph_in);
    power_show2(ret_Ybus.graph, precision, iter, (Ybus_time + fac_time + rhs_time + back_for_time + convert_time) * 1000, total_time, total_time_count);

    // ========== 资源清理 ==========
    CKTSO_DestroySolver(nicslu_p);
    CKTSO_DestroySolver(nicslu_pp);
    free(node);
    free(edge);
    free(off);
    free(graph_in.node);
    free(graph_in.edge);
    free(graph_in.off);
    
    for (int i = 0; i < NUM_KERNEL; i++) {
        free(host_node_P[i]);
        free(host_node_Q[i]);
        free(host_node_Vm[i]);
        free(host_node_Va[i]);
        free(host_node_type[i]);
        free(host_edge_off[i]);
        free(host_edge_src[i]);
        free(host_edge_G[i]);
        free(host_edge_B[i]);
        free(host_deltaP[i]);
        free(host_deltaQ[i]);
    }
    free(host_max_delta);
    free(deltaP);
    free(deltaQ);
    free(deltaP_debug);
    free(deltaQ_debug);

    return EXIT_SUCCESS;
}

void program_exit(std::string casename, double precision, int node_num, int edge_num, int iter,
                 double Ybus_time, double fac_time, double rhs_time, double back_for_time, double convert_time,
                 struct Ret_Ybus ret_Ybus, struct Graph_PF_IN graph_in, double *total_time, int total_time_count) {
    // 实现与main函数中相同的输出逻辑
    printf("----------------information---------------\n");
    printf("dataset                   : %s\n", casename.c_str());
    printf("precision                 : %lf\n", precision);
    printf("nodes                     : %d\n", node_num);
    printf("edges                     : %d\n", edge_num);
    printf("iterations                : %d\n", iter);
    printf("------------------------------------------\n");
    printf("----------------excute time---------------\n");
    printf("initial and create matrix : %lfs\n", Ybus_time/3);
    printf("matrix factorization      : %lfs\n", fac_time/3);
    printf("compute rhs               : %lfs\n", rhs_time/3);
    printf("back and forward          : %lfs\n", back_for_time/3);
    printf("convert                   : %lfs\n", convert_time/3);
    printf("total                     : %lfs\n", (Ybus_time + fac_time + rhs_time + back_for_time + convert_time)/3);
    printf("------------------------------------------\n");

    power_show_csv(ret_Ybus.graph, graph_in);
    power_show2(ret_Ybus.graph, precision, iter, (Ybus_time + fac_time + rhs_time + back_for_time + convert_time) * 1000, total_time, total_time_count);
}