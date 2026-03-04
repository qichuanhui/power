#include <string>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <vector>
#define PI  3.1415926535898
using namespace std;
/*struct Node_PF_In{
    int bus_id;
    int bus_name;
    int area;

    int loss_zone;
    int type;
    double voltage,angle,load_P,load_Q,generation_p,generation_Q;
    double base_kV,desired_volts,MAX_Q,Min_Q,G,B;
    int control_bus_number;
};*/

struct Node_PF_In{
    int type;
    double Vm,Va,P,Q,G,B;
};

/*struct Edge_PF_In{
    int ap_bus,z_bus,area,zone,circuit,type;
    double R,X,B,Line_Q1,Line_Q2,Line_Q3;
    int control_bus,side;
    double transformer_final_turns_ratio,transformer_final_angle,Min_tap,Max_tap;
    int step_size;
    double Min_volt,Max_volt;
};*/
struct Edge_PF_In{
    int from,to;
    double G,B,hB,K;
    int Kcount;
    double BIJ;
};

struct Node_PF{
	double Vm,Va,P,Q;
    int type;
};
struct Edge_PF{
	double G,B;
    int src;
};
struct PF_PQ{
	double deltaP;
	double deltaQ;
};

struct Node_PF_PQ{
	double P,Q;
};

struct Node_PF_V{
	double Vm,Va;
};

struct Graph_PF_IN{
	struct Node_PF_In* node;
	struct Edge_PF_In* edge;
	unsigned * off;
	int node_num;
	int edge_num;
};

struct Graph_PF{
	struct Node_PF* node;
	struct Edge_PF* edge;
	unsigned* off;
	int node_num;
	int edge_num;
};

/*struct Sort_Id_Vertex{
	int Bp_p, Bpp_p, ep, Vm, Va, Pn, Qn, bustype;
};*/
struct Sort_Id_Edge_All{
	double eG, eB;
	int ei;
	double Bp_x, Bpp_x;
};
struct Ret_Ybus{
	double *Bp_x,*Bpp_x;
	unsigned *Bp_i,*Bpp_i;
	unsigned *Bp_p,*Bpp_p;

	struct Graph_PF graph;
};

#define EDGE_BLOCK_SIZE 4
struct Edge_PF_Block{
	int src[EDGE_BLOCK_SIZE];
	double G[EDGE_BLOCK_SIZE];
	double B[EDGE_BLOCK_SIZE];
};

typedef Eigen::SparseLU<Eigen::SparseMatrix<double>> EigenSolver;

void cordic_code_generate();
void cordic(double angle_para,double *cosine, double *sine);

struct Graph_PF_IN pf_load(string casename);
struct Ret_Ybus pf_Ybus_seq(struct Graph_PF_IN graph_in);
struct Ret_Ybus pf_Ybus_par(struct Graph_PF_IN graph_in);
void pf_fac_seq(double *Bp_x,unsigned *Bp_i,unsigned *Bp_p,
		double *Bpp_x,unsigned *Bpp_i,unsigned *Bpp_p,int node_num,
		EigenSolver** nicslu_p,EigenSolver** nicslu_pp);
void pf_fac_par(double *Bp_x,unsigned *Bp_i,unsigned *Bp_p,
		double *Bpp_x,unsigned *Bpp_i,unsigned *Bpp_p,int node_num,
		EigenSolver** nicslu_p,EigenSolver** nicslu_pp);
void rhs_deltaP_par(struct Node_PF *node,struct Edge_PF *edge,unsigned  *off,int node_num,
		double *deltaP,double *deltaQ,double *cosine_array,double *sine_array,
		bool tri_flag,double *maxDeltaP,double *maxDeltaQ);
void rhs_deltaQ_par(struct Node_PF *node,struct Edge_PF *edge,unsigned  *off,int node_num,
		double *deltaP,double *deltaQ,double *cosine_array,double *sine_array,
		bool tri_flag,double *maxDeltaP,double *maxDeltaQ);
void pf_iter();
void pf_iter_MT();

void power_show_csv(struct Graph_PF graph,struct Graph_PF_IN graph_in);
void power_show2(struct Graph_PF graph,double precision,int iter,double excusion_time,double *total_time,double total_time_count);
void program_exit(string casename,double precision,int node_num,int edge_num,int iter,
		double Ybus_time,double fac_time,double rhs_time,double back_for_time,double convert_time,
		struct Ret_Ybus ret_Ybus,struct Graph_PF_IN graph_in,double *total_time,int total_time_count);


//--sp krnl_pf_1.m_axi_gmem0:DDR[0] --sp krnl_pf_1.m_axi_gmem1:DDR[1] --sp krnl_pf_1.m_axi_gmem2:DDR[2] --sp krnl_pf_1.m_axi_gmem3:DDR[3] --profile_kernel data:all:all:all --kernel_frequency 150
