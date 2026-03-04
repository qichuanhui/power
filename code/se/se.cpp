#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <sys/time.h>
#include "se.h"
#include <math.h>
#include <xrt/xrt_kernel.h>
#include <xrt/xrt_device.h>
#include <xrt/xrt_bo.h>
#include <unistd.h>
#include <vector>
#include <cstring>

#define NUM_KERNEL 4
// HBM Banks requirements
#define MAX_DDR_BANKCOUNT 4

struct Node_SE_CONST { // 32byte
    double P, Q, Ri_V, M_Vm;
};

struct Edge_SE_CONST { // 128byte
    int src, Kcount;
    double K;
    double G, B, hB, BIJ;
    double M_P_TLPF, M_Q_TLPF, Ri_eP, Ri_eQ;
    double M_P_TLPF_reverse, M_Q_TLPF_reverse, Ri_eP_reverse, Ri_eQ_reverse;
    double filling1, filling2;
};

struct Node_SE_V* host_node_V[NUM_KERNEL];
struct Node_SE_CONST* host_node_const[NUM_KERNEL];
int* host_edge_off[NUM_KERNEL];
struct Edge_SE_CONST* host_edge_const[NUM_KERNEL];
struct SE_H_r_PQ* host_H_r_PQ[NUM_KERNEL];
double* host_max_delta;

int host_node_begin[NUM_KERNEL], host_node_end[NUM_KERNEL];

#define GROUP_SIZE 1
#define MAX_ITER 20
int debug = 0;
double total_time[MAX_ITER * 10];
int total_time_count = 0;

int main(int argc, char* argv[]) {
    struct timeval start, finish;
    double cost_time;
    double H_time = 0, fac_time = 0, rhs_time = 0, back_for_time = 0;
    double convert_time = 0;

    if (argc != 4) {
        std::cout << "Usage: " << argv[0] << " <xclbin> <casename> <precision>" << std::endl;
        return EXIT_FAILURE;
    }
    std::string binaryFile = argv[1];
    std::string casename = argv[2];
    double precision = atof(argv[3]);

    // XRT 初始化
    unsigned int device_index = 1;
    auto device = xrt::device(device_index);
    auto uuid = device.load_xclbin(binaryFile);

    // 创建内核
    auto krnl_se = xrt::kernel(device, uuid, "krnl_se");

    struct Graph_SE_IN graph_in = se_load(casename);
    gettimeofday(&start, 0); // 计时开始
    struct Ret_H ret_H = se_H_seq(graph_in);
    gettimeofday(&finish, 0); // 计时结束
    cost_time = 1000000 * (finish.tv_sec - start.tv_sec) + finish.tv_usec - start.tv_usec;
    printf("initial and create matrix: cost the time is: %lf s.\n", cost_time / (double)1000000);
    H_time = cost_time / (double)1000000;
    total_time[total_time_count++] = cost_time / (double)1000000 * 1000;

    gettimeofday(&start, 0); // 计时开始
    ICktSo nicslu_p = nullptr, nicslu_pp = nullptr;
    se_fac_seq(ret_H.Bp_x, ret_H.Bp_i, ret_H.Bp_p,
        ret_H.Bpp_x, ret_H.Bpp_i, ret_H.Bpp_p,
        ret_H.graph.node_num, &nicslu_p, &nicslu_pp);
    gettimeofday(&finish, 0); // 计时结束
    cost_time = 1000000 * (finish.tv_sec - start.tv_sec) + finish.tv_usec - start.tv_usec;
    printf("matrix factorization: cost the time is: %lf s.\n", cost_time / (double)1000000);
    fac_time = cost_time / (double)1000000;
    total_time[total_time_count++] = cost_time / (double)1000000 * 1000;

    struct Node_SE* node = ret_H.graph.node;
    struct Edge_SE* edge = ret_H.graph.edge;
    unsigned* off = ret_H.graph.off;
    int node_num = ret_H.graph.node_num;
    int edge_num = ret_H.graph.edge_num;

    int slackbus = ret_H.slackbus;

    {
        int i = 0;

        host_node_begin[0] = 0;
        while (off[++i] < edge_num / 4);
        i = i / 8 * 8;
        while (off[i] - off[host_node_begin[0]] > 1900) i -= 8;
        host_node_end[0] = i;

        host_node_begin[1] = i;
        while (off[++i] < 2 * edge_num / 4);
        i = i / 8 * 8;
        while (off[i] - off[host_node_begin[1]] > 1900) i -= 8;
        host_node_end[1] = i;

        host_node_begin[2] = i;
        while (off[++i] < 3 * edge_num / 4);
        i = i / 8 * 8;
        while (off[i] - off[host_node_begin[2]] > 1900) i -= 8;
        host_node_end[2] = i;

        host_node_begin[3] = i;
        host_node_end[3] = node_num;
    }

    gettimeofday(&start, 0); // 计时开始
    double* H_r_P, * H_r_Q;
    H_r_P = (double*)malloc(node_num * sizeof(double));
    H_r_Q = (double*)malloc(node_num * sizeof(double));

    double* debug_H_r_P, * debug_H_r_Q;
    debug_H_r_P = (double*)malloc(node_num * sizeof(double));
    debug_H_r_Q = (double*)malloc(node_num * sizeof(double));

    int posix_ret;
    for (int i = 0; i < NUM_KERNEL; i++) {
        int l_node_num = host_node_end[i] - host_node_begin[i];

        const size_t node_V_len = l_node_num * sizeof(struct Node_SE_V);
        posix_ret = posix_memalign((void**)&host_node_V[i], 4096, node_V_len);

        const size_t node_const_len = l_node_num * sizeof(struct Node_SE_CONST);
        posix_ret = posix_memalign((void**)&host_node_const[i], 4096, node_const_len);

        const size_t edge_off_len = (l_node_num + 1) * sizeof(int);
        posix_ret = posix_memalign((void**)&host_edge_off[i], 4096, edge_off_len);
        memset(&host_edge_off[i][l_node_num + 1], 0, edge_off_len - (l_node_num + 1) * 4);

        int l_edge_num = off[host_node_end[i]] - off[host_node_begin[i]];

        const size_t edge_const_len = l_edge_num * sizeof(struct Edge_SE_CONST);
        posix_ret = posix_memalign((void**)&host_edge_const[i], 4096, edge_const_len);

        const size_t node_deltaPQ_len = l_node_num * sizeof(struct SE_H_r_PQ);
        posix_ret = posix_memalign((void**)&host_H_r_PQ[i], 4096, node_deltaPQ_len);
    }
    posix_ret = posix_memalign((void**)&host_max_delta, 4096, 64);
    memset(host_max_delta, 0, 64);

    for (int i = 0; i < NUM_KERNEL; i++) {
        int l_node_begin = host_node_begin[i];
        int l_node_num = host_node_end[i] - host_node_begin[i];
        for (int j = 0; j < l_node_num; j++) {
            host_node_V[i][j].Vm = node[l_node_begin + j].Vm;
            host_node_V[i][j].Va = node[l_node_begin + j].Va;
            host_edge_off[i][j] = off[l_node_begin + j];

            host_node_const[i][j].P = node[l_node_begin + j].P;
            host_node_const[i][j].Q = node[l_node_begin + j].Q;
            host_node_const[i][j].M_Vm = node[l_node_begin + j].M_Vm;
            host_node_const[i][j].Ri_V = node[l_node_begin + j].Ri_V;
        }
        host_edge_off[i][l_node_num] = off[l_node_begin + l_node_num];
    }
    gettimeofday(&finish, 0); // 计时结束
    cost_time = 1000000 * (finish.tv_sec - start.tv_sec) + finish.tv_usec - start.tv_usec;

    for (int i = 0; i < NUM_KERNEL; i++) {
        int l_edge_begin = off[host_node_begin[i]];
        int l_edge_num = off[host_node_end[i]] - off[host_node_begin[i]];
        for (int j = 0; j < l_edge_num; j++) {
            host_edge_const[i][j].src = edge[l_edge_begin + j].src;
            host_edge_const[i][j].Kcount = edge[l_edge_begin + j].Kcount;
            host_edge_const[i][j].K = edge[l_edge_begin + j].K;

            host_edge_const[i][j].G = edge[l_edge_begin + j].G;
            host_edge_const[i][j].B = edge[l_edge_begin + j].B;
            host_edge_const[i][j].hB = edge[l_edge_begin + j].hB;
            host_edge_const[i][j].BIJ = edge[l_edge_begin + j].BIJ;

            host_edge_const[i][j].M_P_TLPF = edge[l_edge_begin + j].M_P_TLPF;
            host_edge_const[i][j].M_Q_TLPF = edge[l_edge_begin + j].M_Q_TLPF;
            host_edge_const[i][j].Ri_eP = edge[l_edge_begin + j].Ri_eP;
            host_edge_const[i][j].Ri_eQ = edge[l_edge_begin + j].Ri_eQ;
            host_edge_const[i][j].M_P_TLPF_reverse = edge[l_edge_begin + j].M_P_TLPF_reverse;
            host_edge_const[i][j].M_Q_TLPF_reverse = edge[l_edge_begin + j].M_Q_TLPF_reverse;
            host_edge_const[i][j].Ri_eP_reverse = edge[l_edge_begin + j].Ri_eP_reverse;
            host_edge_const[i][j].Ri_eQ_reverse = edge[l_edge_begin + j].Ri_eQ_reverse;

            host_edge_const[i][j].filling1 = 0;
            host_edge_const[i][j].filling2 = 0;
        }
    }

    // XRT 缓冲区对象
    std::vector<xrt::bo> buffer_node_V;
    std::vector<xrt::bo> buffer_node_const;
    std::vector<xrt::bo> buffer_edge_off;
    std::vector<xrt::bo> buffer_edge_const;
    std::vector<xrt::bo> buffer_H_r_PQ;
    xrt::bo buffer_max_delta;

    // 创建缓冲区
    for (int i = 0; i < NUM_KERNEL; i++) {
        int l_node_num = host_node_end[i] - host_node_begin[i];
        const size_t node_V_len = l_node_num * sizeof(struct Node_SE_V);
        buffer_node_V.push_back(xrt::bo(device, node_V_len, krnl_se.group_id(i)));

        const size_t node_const_len = l_node_num * sizeof(struct Node_SE_CONST);
        buffer_node_const.push_back(xrt::bo(device, node_const_len, krnl_se.group_id(i)));

        const size_t edge_off_len = (l_node_num + 1) * sizeof(int);
        buffer_edge_off.push_back(xrt::bo(device, edge_off_len, krnl_se.group_id(i)));

        int l_edge_num = off[host_node_end[i]] - off[host_node_begin[i]];
        const size_t edge_const_len = l_edge_num * sizeof(struct Edge_SE_CONST);
        buffer_edge_const.push_back(xrt::bo(device, edge_const_len, krnl_se.group_id(i)));

        const size_t node_deltaPQ_len = l_node_num * sizeof(struct SE_H_r_PQ);
        buffer_H_r_PQ.push_back(xrt::bo(device, node_deltaPQ_len, krnl_se.group_id(i)));
    }
    buffer_max_delta = xrt::bo(device, 64, krnl_se.group_id(0));

    // 数据映射和传输
    auto transfer_data = [&](bool first_iter) {
        for (int i = 0; i < NUM_KERNEL; i++) {
            int l_node_num = host_node_end[i] - host_node_begin[i];
            const size_t node_V_len = l_node_num * sizeof(struct Node_SE_V);
            auto ptr_node_V = buffer_node_V[i].map<struct Node_SE_V*>();
            memcpy(ptr_node_V, host_node_V[i], node_V_len);
            buffer_node_V[i].sync(XCL_BO_SYNC_BO_TO_DEVICE);

            if (first_iter) {
                const size_t node_const_len = l_node_num * sizeof(struct Node_SE_CONST);
                auto ptr_node_const = buffer_node_const[i].map<struct Node_SE_CONST*>();
                memcpy(ptr_node_const, host_node_const[i], node_const_len);
                buffer_node_const[i].sync(XCL_BO_SYNC_BO_TO_DEVICE);

                const size_t edge_off_len = (l_node_num + 1) * sizeof(int);
                auto ptr_edge_off = buffer_edge_off[i].map<int*>();
                memcpy(ptr_edge_off, host_edge_off[i], edge_off_len);
                buffer_edge_off[i].sync(XCL_BO_SYNC_BO_TO_DEVICE);

                int l_edge_num = off[host_node_end[i]] - off[host_node_begin[i]];
                const size_t edge_const_len = l_edge_num * sizeof(struct Edge_SE_CONST);
                auto ptr_edge_const = buffer_edge_const[i].map<struct Edge_SE_CONST*>();
                memcpy(ptr_edge_const, host_edge_const[i], edge_const_len);
                buffer_edge_const[i].sync(XCL_BO_SYNC_BO_TO_DEVICE);
            }
        }
    };

    // 运行内核
    auto run_kernel = [&](bool init_flag) {
        auto run = krnl_se(
            buffer_node_V[0],     // arg0
            buffer_node_V[1],     // arg1
            buffer_node_V[2],     // arg2
            buffer_node_V[3],     // arg3
            buffer_node_const[0], // arg4
            buffer_node_const[1], // arg5
            buffer_node_const[2], // arg6
            buffer_node_const[3], // arg7
            buffer_edge_off[0],   // arg8
            buffer_edge_off[1],   // arg9
            buffer_edge_off[2],   // arg10
            buffer_edge_off[3],   // arg11
            buffer_edge_const[0], // arg12
            buffer_edge_const[1], // arg13
            buffer_edge_const[2], // arg14
            buffer_edge_const[3], // arg15
            buffer_H_r_PQ[0],     // arg16
            buffer_H_r_PQ[1],     // arg17
            buffer_H_r_PQ[2],     // arg18
            buffer_H_r_PQ[3],     // arg19
            host_node_begin[0],   // arg20
            host_node_begin[1],   // arg21
            host_node_begin[2],   // arg22
            host_node_begin[3],   // arg23
            host_node_end[0],     // arg24
            host_node_end[1],     // arg25
            host_node_end[2],     // arg26
            host_node_end[3],     // arg27
            node_num,             // arg28
            init_flag             // arg29
        );
        run.wait();
    };

    // 读取结果
    auto read_results = [&]() {
        for (int i = 0; i < NUM_KERNEL; i++) {
            buffer_H_r_PQ[i].sync(XCL_BO_SYNC_BO_FROM_DEVICE);
            auto ptr_H_r_PQ = buffer_H_r_PQ[i].map<struct SE_H_r_PQ*>();

            int l_node_begin = host_node_begin[i];
            int l_node_end = host_node_end[i];
            int l_node_num = l_node_end - l_node_begin;

            for (int j = 0; j < l_node_num; j++) {
                H_r_P[l_node_begin + j] = ptr_H_r_PQ[j].H_r_P;
                H_r_Q[l_node_begin + j] = ptr_H_r_PQ[j].H_r_Q;
            }
        }
        H_r_P[slackbus] = 0;
    };

    double max_dVm = 12345678;
    double max_dVa = 12345678;
    int iter;
    for (iter = 0; iter < MAX_ITER && (max_dVm > precision || max_dVa > precision); ++iter) {
        max_dVm = 0;
        max_dVa = 0;

        if (true) {
            // 更新主机数据
            for (int i = 0; i < NUM_KERNEL; i++) {
                int l_node_begin = host_node_begin[i];
                int l_node_end = host_node_end[i];
                int l_node_num = l_node_end - l_node_begin;
                for (int j = 0; j < l_node_num; j++) {
                    host_node_V[i][j].Vm = node[l_node_begin + j].Vm;
                    host_node_V[i][j].Va = node[l_node_begin + j].Va;
                }
            }

            gettimeofday(&start, 0);
            transfer_data(iter == 0);
            gettimeofday(&finish, 0);
            cost_time = 1000000 * (finish.tv_sec - start.tv_sec) + finish.tv_usec - start.tv_usec;
            printf("iter %d: write device: cost time: %lf s\n", iter / 2, cost_time / (double)1000000);
            convert_time += cost_time / (double)1000000;
            total_time[total_time_count++] = cost_time / (double)1000000 * 1000;

            gettimeofday(&start, 0);
            run_kernel(iter == 0);
            gettimeofday(&finish, 0);
            cost_time = 1000000 * (finish.tv_sec - start.tv_sec) + finish.tv_usec - start.tv_usec;
            printf("iter %d: compute rhs: cost time: %lf s\n", iter / 2, cost_time / (double)1000000);
            rhs_time += cost_time / (double)1000000;
            total_time[total_time_count++] = cost_time / (double)1000000 * 1000;

            gettimeofday(&start, 0);
            read_results();
            gettimeofday(&finish, 0);
            cost_time = 1000000 * (finish.tv_sec - start.tv_sec) + finish.tv_usec - start.tv_usec;
            printf("iter %d: read device: cost time: %lf s\n", iter / 2, cost_time / (double)1000000);
            convert_time += cost_time / (double)1000000;
            total_time[total_time_count++] = cost_time / (double)1000000 * 1000;

            if (debug == 1) {
                for (int i = 0; i < node_num; ++i) {
                    double deltaP = 0, deltaQ = 0;
                    double tmp_H_r_P = 0, tmp_H_r_Q = 0;
                    for (int j = off[i]; j < off[i + 1]; ++j) {
                        struct Edge_SE e = edge[j];
                        int p = e.src;

                        double tap_ratio = fabs(e.K / e.Kcount);
                        double tap_ratio_square = (e.K / e.Kcount) * (e.K / e.Kcount);

                        if (e.K == 0) {
                            tmp_H_r_P += e.BIJ * (e.M_P_TLPF - (node[i].Vm * node[i].Vm * e.G - node[i].Vm * node[p].Vm * (e.G * cos(node[i].Va - node[p].Va) + (-e.B) * sin(node[i].Va - node[p].Va)))) * e.Ri_eP;
                            tmp_H_r_Q += (e.B - e.hB) * (e.M_Q_TLPF - (-node[i].Vm * node[i].Vm * (-e.B + 0.5 * e.hB) - node[i].Vm * node[p].Vm * (e.G * sin(node[i].Va - node[p].Va) - (-e.B) * cos(node[i].Va - node[p].Va)))) * e.Ri_eQ;
                        }
                        else if (e.K > 0) {
                            tmp_H_r_P += e.BIJ * (e.M_P_TLPF - (node[i].Vm * node[i].Vm * (e.G / tap_ratio_square) - node[i].Vm * node[p].Vm * ((e.G / tap_ratio) * cos(node[i].Va - node[p].Va) + (-e.B / tap_ratio) * sin(node[i].Va - node[p].Va)))) * e.Ri_eP;
                            tmp_H_r_Q += (e.B / tap_ratio - e.hB) * (e.M_Q_TLPF - (-node[i].Vm * node[i].Vm * (-e.B + 0.5 * e.hB) / tap_ratio_square - node[i].Vm * node[p].Vm * ((e.G / tap_ratio) * sin(node[i].Va - node[p].Va) - (-e.B / tap_ratio) * cos(node[i].Va - node[p].Va)))) * e.Ri_eQ;
                        }
                        else {
                            tmp_H_r_P += e.BIJ * (e.M_P_TLPF - (node[i].Vm * node[i].Vm * e.G - node[i].Vm * node[p].Vm * ((e.G / tap_ratio) * cos(node[i].Va - node[p].Va) + (-e.B / tap_ratio) * sin(node[i].Va - node[p].Va)))) * e.Ri_eP;
                            tmp_H_r_Q += (e.B / tap_ratio - e.hB) * (e.M_Q_TLPF - (-node[i].Vm * node[i].Vm * (-e.B + 0.5 * e.hB) - node[i].Vm * node[p].Vm * ((e.G / tap_ratio) * sin(node[i].Va - node[p].Va) - (-e.B / tap_ratio) * cos(node[i].Va - node[p].Va)))) * e.Ri_eQ;
                        }

                        if ((-e.K) == 0) {
                            tmp_H_r_P += (-1) * e.BIJ * (e.M_P_TLPF_reverse - (node[p].Vm * node[p].Vm * e.G - node[p].Vm * node[i].Vm * (e.G * cos(node[p].Va - node[i].Va) + (-e.B) * sin(node[p].Va - node[i].Va)))) * e.Ri_eP_reverse;
                            tmp_H_r_Q += (-1) * e.B * (e.M_Q_TLPF_reverse - (-node[p].Vm * node[p].Vm * (-e.B + 0.5 * e.hB) - node[p].Vm * node[i].Vm * (e.G * sin(node[p].Va - node[i].Va) - (-e.B) * cos(node[p].Va - node[i].Va)))) * e.Ri_eQ_reverse;
                        }
                        else if ((-e.K) > 0) {
                            tmp_H_r_P += (-1) * e.BIJ * (e.M_P_TLPF_reverse - (node[p].Vm * node[p].Vm * (e.G / tap_ratio_square) - node[p].Vm * node[i].Vm * ((e.G / tap_ratio) * cos(node[p].Va - node[i].Va) + (-e.B / tap_ratio) * sin(node[p].Va - node[i].Va)))) * e.Ri_eP_reverse;
                            tmp_H_r_Q += (-1) * (e.B / tap_ratio) * (e.M_Q_TLPF_reverse - (-node[p].Vm * node[p].Vm * (-e.B + 0.5 * e.hB) / tap_ratio_square - node[p].Vm * node[i].Vm * ((e.G / tap_ratio) * sin(node[p].Va - node[i].Va) - (-e.B / tap_ratio) * cos(node[p].Va - node[i].Va)))) * e.Ri_eQ_reverse;
                        }
                        else {
                            tmp_H_r_P += (-1) * e.BIJ * (e.M_P_TLPF_reverse - (node[p].Vm * node[p].Vm * e.G - node[p].Vm * node[i].Vm * ((e.G / tap_ratio) * cos(node[p].Va - node[i].Va) + (-e.B / tap_ratio) * sin(node[p].Va - node[i].Va)))) * e.Ri_eP_reverse;
                            tmp_H_r_Q += (-1) * (e.B / tap_ratio) * (e.M_Q_TLPF_reverse - (-node[p].Vm * node[p].Vm * (-e.B + 0.5 * e.hB) - node[p].Vm * node[i].Vm * ((e.G / tap_ratio) * sin(node[p].Va - node[i].Va) - (-e.B / tap_ratio) * cos(node[p].Va - node[i].Va)))) * e.Ri_eQ_reverse;
                        }
                    }
                    deltaP = node[i].P - (deltaP + node[i].Vm * node[i].Vm * node[i].sumG);
                    tmp_H_r_P += (-1) * node[i].sumBi * deltaP * node[i].Ri_vP;
                    if (node[i].type != 3) {
                        debug_H_r_P[i] = tmp_H_r_P;
                    }
                    else {
                        debug_H_r_P[i] = 0;
                    }

                    deltaQ = node[i].Q - (deltaQ - node[i].Vm * node[i].Vm * node[i].sumB);
                    double deltaVm = node[i].M_Vm - node[i].Vm;
                    tmp_H_r_Q += (-1) * (2 * node[i].sumB + node[i].sumBedge) * deltaQ * node[i].Ri_vQ + deltaVm * node[i].Ri_V;
                    debug_H_r_Q[i] = tmp_H_r_Q;
                }
            }

            if (debug == 1) {
                int error_count = 0;
                for (int j = 0; j < node_num; j++) {
                    if ((fabs(debug_H_r_P[j] - H_r_P[j]) > 0.001 || fabs(debug_H_r_Q[j] - H_r_Q[j]) > 0.001)
                        && error_count < 20) {
                        printf("debug_H_r_P[%d]=%lf,debug_H_r_Q[%d]=%lf,H_r_P[%d]=%lf,H_r_Q[%d]=%lf\n",
                            j, debug_H_r_P[j], j, debug_H_r_Q[j],
                            j, H_r_P[j], j, H_r_Q[j]);
                        error_count++;
                    }
                    H_r_P[j] = debug_H_r_P[j];
                    H_r_Q[j] = debug_H_r_Q[j];
                }
            }
        }

        gettimeofday(&start, 0); // 计时开始
        if (iter % 2 == 0) {
            CKTSO_Solve(nicslu_p, H_r_P, H_r_P, false, false);
            max_dVa = 0;
            for (int i = 0; i < node_num; i++) {
                if (fabs(H_r_P[i]) > max_dVa) {
                    max_dVa = fabs(H_r_P[i]);
                }
                node[i].Va += H_r_P[i];
            }
        }
        else {
            CKTSO_Solve(nicslu_pp, H_r_Q, H_r_Q, false, false);
            max_dVm = 0;
            for (int i = 0; i < node_num; i++) {
                if (fabs(H_r_Q[i]) > max_dVm) {
                    max_dVm = fabs(H_r_Q[i]);
                }
                node[i].Vm += H_r_Q[i];
            }
        }

        gettimeofday(&finish, 0); // 计时结束
        cost_time = 1000000 * (finish.tv_sec - start.tv_sec) + finish.tv_usec - start.tv_usec;
        if ((iter % 2) == 0)
            printf("back_forward and update Va: cost the time is: %lf s.\n", cost_time / (double)1000000);
        else
            printf("back_forward and update Vm: cost the time is: %lf s.\n", cost_time / (double)1000000);
        back_for_time += cost_time / (double)1000000;
        total_time[total_time_count++] = cost_time / (double)1000000 * 1000;

        printf("iter %d: max_dVa=%lf, max_dVm=%lf\n", iter / 2, max_dVa, max_dVm);
    }

    printf("----------------information---------------\n");
    printf("dataset                   : %s\n", casename.c_str());
    printf("precision                 : %lf\n", precision);
    printf("nodes                     : %d\n", node_num);
    printf("edges                     : %d\n", edge_num);
    printf("iterations                : %d\n", (iter) / 2);
    printf("------------------------------------------\n");
    printf("----------------excute time---------------\n");
    printf("initial and create matrix : %lfs\n", H_time*100);
    printf("matrix factorization      : %lfs\n", fac_time*100);
    printf("compute rhs               : %lfs\n", rhs_time*100);
    printf("back and forward          : %lfs\n", back_for_time*100);
    printf("convert                   : %lfs\n", convert_time*100);
    printf("total                     : %lfs\n", (H_time + fac_time + rhs_time + back_for_time + convert_time)*100);
    printf("------------------------------------------\n");

    power_show(ret_H.graph, graph_in);
    power_show2(ret_H.graph, precision, (iter + 1) / 2,
                (H_time + fac_time + rhs_time + back_for_time + convert_time) * 1000,
                total_time, total_time_count);

} // <-- 这里添加了缺失的闭合大括号
