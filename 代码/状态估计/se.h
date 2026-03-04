//#ifndef SE_H_INCLUDED
//#define SE_H_INCLUDED

#include "cktso.h"
#include <string>
using namespace std;

#define PI  3.1415926535898
struct Node_SE_In{
    int type;
    double Vm, Va, P, Q, G, B, Ri_V, Ri_vP, Ri_vQ;
    double M_Vm;
};

struct Edge_SE_In{
    int from, to;
    double G, B, hB, K;
    int Kcount;
    double BIJ, M_P_TLPF, M_Q_TLPF, Ri_eP, Ri_eQ;
    double M_P_TLPF_reverse,M_Q_TLPF_reverse,Ri_eP_reverse,Ri_eQ_reverse;
};
/*
struct Edge_PF_In{
    int from,to;
    double G,B,hB,K;
    int Kcount;
    double BIJ;
};
*/

struct Node_SE{
	int type;
    double Vm, Va, P, Q, Ri_V, Ri_vP, Ri_vQ;
    double sumG,sumB,sumBedge,sumBi;
    double M_Vm;
};
struct Edge_SE{
    int src;
    double G, B, hB, K;
    int Kcount;
    double BIJ, M_P_TLPF, M_Q_TLPF, Ri_eP, Ri_eQ;
    double M_P_TLPF_reverse,M_Q_TLPF_reverse,Ri_eP_reverse,Ri_eQ_reverse;
};

struct Graph_SE_IN{
	struct Node_SE_In* node;
	struct Edge_SE_In* edge;
	unsigned * off;
	int node_num;
	int edge_num;
};

struct Graph_SE{
	struct Node_SE* node;
	struct Edge_SE* edge;
	unsigned* off;
	int node_num;
	int edge_num;
};

struct Ret_H{
	double *Bp_x,*Bpp_x;
	unsigned *Bp_i,*Bpp_i;
	unsigned *Bp_p,*Bpp_p;

	struct Graph_SE graph;
	int slackbus;
};

struct Node_SE_V{
	double Vm,Va;
};
struct SE_H_r_PQ{
	double H_r_P;
	double H_r_Q;
};

struct Graph_SE_IN se_load(string casename);
struct Ret_H se_H_seq(struct Graph_SE_IN graph_in);

void *se_fac_seq(double *Bp_x,unsigned *Bp_i,unsigned *Bp_p,
		double *Bpp_x,unsigned *Bpp_i,unsigned *Bpp_p,int node_num,
		ICktSo* nicslu_p,ICktSo* nicslu_pp);
void *se_fac_par(double *Bp_x,unsigned *Bp_i,unsigned *Bp_p,
		double *Bpp_x,unsigned *Bpp_i,unsigned *Bpp_p,int node_num,
		ICktSo* nicslu_p,ICktSo* nicslu_pp);

int se_iter();
int se_iter_time();
void power_show(struct Graph_SE graph,struct Graph_SE_IN graph_in);
void power_show2(struct Graph_SE graph,double precision,int iter,double excusion_time,double *total_time,double total_time_count);


//#endif // SE_H_INCLUDED
