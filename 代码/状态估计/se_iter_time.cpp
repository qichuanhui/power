#include <fstream>
#include <iostream>
#include <sys/time.h>
#include "se.h"
#include <math.h>
#include <string.h>

using namespace std;

#define MAX_ITER 50

int se_iter_time()
{
    struct timeval start, finish;
    double H_time = 0, fac_time = 0, rhs_time = 0, back_for_time = 0;

    struct Graph_SE_IN graph_in = se_load("sc_20171208_11100");
    
    gettimeofday(&start, 0);
    struct Ret_H ret_H = se_H_seq(graph_in);
    gettimeofday(&finish, 0);
    H_time = 1000000 * (finish.tv_sec - start.tv_sec) + finish.tv_usec - start.tv_usec;
    H_time /= 1000000;

    gettimeofday(&start, 0);
    ICktSo nicslu_p = nullptr, nicslu_pp = nullptr;
    se_fac_seq(ret_H.Bp_x, ret_H.Bp_i, ret_H.Bp_p,
               ret_H.Bpp_x, ret_H.Bpp_i, ret_H.Bpp_p,
               ret_H.graph.node_num, &nicslu_p, &nicslu_pp);
    gettimeofday(&finish, 0);
    fac_time = 1000000 * (finish.tv_sec - start.tv_sec) + finish.tv_usec - start.tv_usec;
    fac_time /= 1000000;

    struct Node_SE *node = ret_H.graph.node;
    struct Edge_SE *edge = ret_H.graph.edge;
    unsigned *off = ret_H.graph.off;
    int node_num = ret_H.graph.node_num;
    int edge_num = ret_H.graph.edge_num;
    int slackbus = ret_H.slackbus;

    double *H_r_P, *H_r_Q;
    H_r_P = (double *)malloc(node_num * sizeof(double));
    H_r_Q = (double *)malloc(node_num * sizeof(double));

    double max_dVm = 12345678, max_dVa = 12345678;
    int iter;
    for (iter = 0; iter < MAX_ITER && max_dVm > 0.001 && max_dVa > 0.001; ++iter)
    {
        gettimeofday(&start, 0);
        for (int i = 0; i < node_num; ++i)
        {
            double deltaP = 0, deltaQ = 0;
            double tmp_H_r_P = 0, tmp_H_r_Q = 0;
            for (int j = off[i]; j < off[i + 1]; ++j)
            {
                struct Edge_SE e = edge[j];
                int p = e.src;
                double tap_ratio = fabs(e.K / e.Kcount);
                double tap_ratio_square = (e.K / e.Kcount) * (e.K / e.Kcount);

                if (e.K == 0)
                {
                    tmp_H_r_P += e.BIJ * (e.M_P_TLPF - (node[i].Vm * node[i].Vm * e.G - node[i].Vm * node[p].Vm * (e.G * cos(node[i].Va - node[p].Va) + (-e.B) * sin(node[i].Va - node[p].Va)))) * e.Ri_eP;
                    tmp_H_r_Q += (e.B - e.hB) * (e.M_Q_TLPF - (-node[i].Vm * node[i].Vm * (-e.B + 0.5 * e.hB) - node[i].Vm * node[p].Vm * (e.G * sin(node[i].Va - node[p].Va) - (-e.B) * cos(node[i].Va - node[p].Va)))) * e.Ri_eQ;
                }
                else if (e.K > 0)
                {
                    tmp_H_r_P += e.BIJ * (e.M_P_TLPF - (node[i].Vm * node[i].Vm * (e.G / tap_ratio_square) - node[i].Vm * node[p].Vm * ((e.G / tap_ratio) * cos(node[i].Va - node[p].Va) + (-e.B / tap_ratio) * sin(node[i].Va - node[p].Va)))) * e.Ri_eP;
                    tmp_H_r_Q += (e.B / tap_ratio - e.hB) * (e.M_Q_TLPF - (-node[i].Vm * node[i].Vm * (-e.B + 0.5 * e.hB) / tap_ratio_square - node[i].Vm * node[p].Vm * ((e.G / tap_ratio) * sin(node[i].Va - node[p].Va) - (-e.B / tap_ratio) * cos(node[i].Va - node[p].Va)))) * e.Ri_eQ;
                }
                else
                {
                    tmp_H_r_P += e.BIJ * (e.M_P_TLPF - (node[i].Vm * node[i].Vm * e.G - node[i].Vm * node[p].Vm * ((e.G / tap_ratio) * cos(node[i].Va - node[p].Va) + (-e.B / tap_ratio) * sin(node[i].Va - node[p].Va)))) * e.Ri_eP;
                    tmp_H_r_Q += (e.B / tap_ratio - e.hB) * (e.M_Q_TLPF - (-node[i].Vm * node[i].Vm * (-e.B + 0.5 * e.hB) - node[i].Vm * node[p].Vm * ((e.G / tap_ratio) * sin(node[i].Va - node[p].Va) - (-e.B / tap_ratio) * cos(node[i].Va - node[p].Va)))) * e.Ri_eQ;
                }

                if ((-e.K) == 0)
                {
                    tmp_H_r_P += (-1) * e.BIJ * (e.M_P_TLPF_reverse - (node[p].Vm * node[p].Vm * e.G - node[p].Vm * node[i].Vm * (e.G * cos(node[p].Va - node[i].Va) + (-e.B) * sin(node[p].Va - node[i].Va)))) * e.Ri_eP_reverse;
                    tmp_H_r_Q += (-1) * e.B * (e.M_Q_TLPF_reverse - (-node[p].Vm * node[p].Vm * (-e.B + 0.5 * e.hB) - node[p].Vm * node[i].Vm * (e.G * sin(node[p].Va - node[i].Va) - (-e.B) * cos(node[p].Va - node[i].Va)))) * e.Ri_eQ_reverse;
                }
                else if ((-e.K) > 0)
                {
                    tmp_H_r_P += (-1) * e.BIJ * (e.M_P_TLPF_reverse - (node[p].Vm * node[p].Vm * (e.G / tap_ratio_square) - node[p].Vm * node[i].Vm * ((e.G / tap_ratio) * cos(node[p].Va - node[i].Va) + (-e.B / tap_ratio) * sin(node[p].Va - node[i].Va)))) * e.Ri_eP_reverse;
                    tmp_H_r_Q += (-1) * (e.B / tap_ratio) * (e.M_Q_TLPF_reverse - (-node[p].Vm * node[p].Vm * (-e.B + 0.5 * e.hB) / tap_ratio_square - node[p].Vm * node[i].Vm * ((e.G / tap_ratio) * sin(node[p].Va - node[i].Va) - (-e.B / tap_ratio) * cos(node[p].Va - node[i].Va)))) * e.Ri_eQ_reverse;
                }
                else
                {
                    tmp_H_r_P += (-1) * e.BIJ * (e.M_P_TLPF_reverse - (node[p].Vm * node[p].Vm * e.G - node[p].Vm * node[i].Vm * ((e.G / tap_ratio) * cos(node[p].Va - node[i].Va) + (-e.B / tap_ratio) * sin(node[p].Va - node[i].Va)))) * e.Ri_eP_reverse;
                    tmp_H_r_Q += (-1) * (e.B / tap_ratio) * (e.M_Q_TLPF_reverse - (-node[p].Vm * node[p].Vm * (-e.B + 0.5 * e.hB) - node[p].Vm * node[i].Vm * ((e.G / tap_ratio) * sin(node[p].Va - node[i].Va) - (-e.B / tap_ratio) * cos(node[p].Va - node[i].Va)))) * e.Ri_eQ_reverse;
                }
            }
            deltaP = node[i].P - (deltaP + node[i].Vm * node[i].Vm * node[i].sumG);
            tmp_H_r_P += (-1) * node[i].sumBi * deltaP * node[i].Ri_vP;
            if (node[i].type != 3)
            {
                H_r_P[i] = tmp_H_r_P;
            }
            else
            {
                H_r_P[i] = 0;
            }

            deltaQ = node[i].Q - (deltaQ - node[i].Vm * node[i].Vm * node[i].sumB);
            double deltaVm = node[i].M_Vm - node[i].Vm;
            tmp_H_r_Q += (-1) * (2 * node[i].sumB + node[i].sumBedge) * deltaQ * node[i].Ri_vQ + deltaVm * node[i].Ri_V;
            H_r_Q[i] = tmp_H_r_Q;
        }
        gettimeofday(&finish, 0);
        rhs_time += 1000000 * (finish.tv_sec - start.tv_sec) + finish.tv_usec - start.tv_usec;

        if (iter % 2 == 0)
        {
            gettimeofday(&start, 0);
            CKTSO_Solve(nicslu_p, H_r_P, H_r_P, false, false);
            gettimeofday(&finish, 0);
            back_for_time += 1000000 * (finish.tv_sec - start.tv_sec) + finish.tv_usec - start.tv_usec;

            max_dVa = 0;
            for (int i = 0; i < node_num; i++)
            {
                if (fabs(H_r_P[i]) > max_dVa)
                {
                    max_dVa = fabs(H_r_P[i]);
                }
                node[i].Va += H_r_P[i];
            }
        }
        else
        {
            gettimeofday(&start, 0);
            CKTSO_Solve(nicslu_pp, H_r_Q, H_r_Q, false, false);
            gettimeofday(&finish, 0);
            back_for_time += 1000000 * (finish.tv_sec - start.tv_sec) + finish.tv_usec - start.tv_usec;

            max_dVm = 0;
            for (int i = 0; i < node_num; i++)
            {
                if (fabs(H_r_Q[i]) > max_dVm)
                {
                    max_dVm = fabs(H_r_Q[i]);
                }
                node[i].Vm += H_r_Q[i];
            }
            printf("krnl:max_dVa=%lf,max_dVm=%lf\n", max_dVa, max_dVm);
        }
    }

    for (int i = 0; i < 10; i++)
        printf("(%d,%lf,%lf,%lf)\n", i, node[i].Va, node[i].Vm, node[i].M_Vm);

    rhs_time /= 1000000;
    back_for_time /= 1000000;

    printf("Matrix Generation Time: %lf s\n", H_time);
    printf("Matrix Factorization Time: %lf s\n", fac_time);
    printf("Right-Hand Side Calculation Time: %lf s\n", rhs_time);
    printf("Backward-Forward Substitution Time: %lf s\n", back_for_time);

    return 0;
}