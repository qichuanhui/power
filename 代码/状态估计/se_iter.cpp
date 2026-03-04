#include <fstream>
#include <iostream>
#include <sys/time.h>
#include "se.h"
#include <math.h>

#include <string.h>
//#include "xcl2.hpp"

using namespace std;

#define MAX_ITER 50

int se_iter()
{
    struct timeval start,finish;
    double cost_time;
    double H_time=0,fac_time=0,rhs_time=0,back_for_time=0;
    double convert_time=0;

    struct Graph_SE_IN graph_in=se_load("sc_20171208_11100");
    gettimeofday(&start, 0);//计时开始
    struct Ret_H ret_H=se_H_seq(graph_in);
    gettimeofday(&finish, 0);//计时结束
    cost_time = 1000000 * (finish.tv_sec - start.tv_sec) + finish.tv_usec - start.tv_usec;
    printf("se_H:cost the time is : %lf s.\n",cost_time/(double)1000000);
    H_time=cost_time/(double)1000000;

    gettimeofday(&start, 0);//计时开始
    ICktSo nicslu_p = nullptr ,nicslu_pp = nullptr;
    se_fac_seq(ret_H.Bp_x,ret_H.Bp_i,ret_H.Bp_p,
    		ret_H.Bpp_x,ret_H.Bpp_i,ret_H.Bpp_p,
			ret_H.graph.node_num,&nicslu_p,&nicslu_pp);
    gettimeofday(&finish, 0);//计时结束
    cost_time = 1000000 * (finish.tv_sec - start.tv_sec) + finish.tv_usec - start.tv_usec;
    printf("se_fac:cost the time is : %lf s.\n",cost_time/(double)1000000);
    fac_time=cost_time/(double)1000000;

    struct Node_SE *node=ret_H.graph.node;
    struct Edge_SE *edge=ret_H.graph.edge;
    unsigned  *off=ret_H.graph.off;
    int node_num=ret_H.graph.node_num;
    int edge_num=ret_H.graph.edge_num;

    int slackbus=ret_H.slackbus;

    double *H_r_P,*H_r_Q;
    H_r_P=(double *)malloc(node_num*sizeof(double));
    H_r_Q=(double *)malloc(node_num*sizeof(double));

    double *debug_H_r_P,*debug_H_r_Q;
    debug_H_r_P=(double *)malloc(node_num*sizeof(double));
    debug_H_r_Q=(double *)malloc(node_num*sizeof(double));

    double max_dVm=12345678,max_dVa=12345678;
    int iter;
    for (iter=0; iter<MAX_ITER&&max_dVm>0.001&&max_dVa>0.001; ++iter){
    	gettimeofday(&start, 0);//计时开始
    	/*for(int i=0;i<node_num;i++){
    		debug_H_r_P[i]=0;
    		debug_H_r_Q[i]=0;
    	}*/

    	for (int i=0; i<node_num; ++i){
    		double deltaP=0,deltaQ=0;
    		double tmp_H_r_P=0,tmp_H_r_Q=0;
    		for (int j=off[i]; j<off[i+1]; ++j){
    			struct Edge_SE e=edge[j];
    			int p=e.src;

    			double tap_ratio = fabs(e.K/e.Kcount);
    			double tap_ratio_square = (e.K/e.Kcount)*(e.K/e.Kcount);

    			/*if(e.K == 0 || fabs(e.K) == 1){
    				deltaP += node[i].Vm*node[p].Vm * (-1*e.G*cos(node[i].Va-node[p].Va) + (e.B * sin(node[i].Va - node[p].Va)));
    				deltaQ += node[i].Vm*node[p].Vm * (-1*e.G*sin(node[i].Va-node[p].Va) - (e.B * cos(node[i].Va-node[p].Va)));
    			}
    			else{
    				double newG = e.G/tap_ratio;
    				double newB = e.B/tap_ratio;
    				deltaP += node[i].Vm*node[p].Vm * (-1*newG*cos(node[i].Va-node[p].Va) + (newB * sin(node[i].Va-node[p].Va)));
    				deltaQ += node[i].Vm*node[p].Vm * (-1*newG*sin(node[i].Va-node[p].Va) - (newB * cos(node[i].Va-node[p].Va)));
    			}*/

    			if(e.K == 0){
    				tmp_H_r_P += e.BIJ * (e.M_P_TLPF - (node[i].Vm * node[i].Vm * e.G - node[i].Vm * node[p].Vm * (e.G*cos(node[i].Va-node[p].Va) + (-e.B)*sin(node[i].Va-node[p].Va)))) * e.Ri_eP;
    				tmp_H_r_Q += (e.B - e.hB) * (e.M_Q_TLPF - (- node[i].Vm * node[i].Vm * (-e.B + 0.5*e.hB) - node[i].Vm * node[p].Vm * (e.G*sin(node[i].Va-node[p].Va) - (-e.B)*cos(node[i].Va-node[p].Va)))) * e.Ri_eQ;
    			//if(i==0&&iter==0) printf("(%d<-(%d,%lf))\n",i,p,e.BIJ * (e.M_P_TLPF - (node[i].Vm * node[i].Vm * e.G - node[i].Vm * node[p].Vm * (e.G*cos(node[i].Va-node[p].Va) + (-e.B)*sin(node[i].Va-node[p].Va)))) * e.Ri_eP);
    			}
    			else if(e.K > 0){
    				tmp_H_r_P += e.BIJ * (e.M_P_TLPF - (node[i].Vm * node[i].Vm * (e.G/tap_ratio_square) - node[i].Vm * node[p].Vm * ((e.G/tap_ratio)*cos(node[i].Va-node[p].Va) + (-e.B/tap_ratio)*sin(node[i].Va-node[p].Va)))) * e.Ri_eP;
    				tmp_H_r_Q += (e.B/tap_ratio - e.hB) * (e.M_Q_TLPF - (- node[i].Vm * node[i].Vm * (-e.B + 0.5*e.hB) / tap_ratio_square - node[i].Vm * node[p].Vm * ((e.G/tap_ratio)*sin(node[i].Va-node[p].Va) - (-e.B/tap_ratio)*cos(node[i].Va-node[p].Va)))) * e.Ri_eQ;
    			//if(i==0&&iter==0) printf("(%d<-(%d,%lf))\n",i,p,e.BIJ * (e.M_P_TLPF - (node[i].Vm * node[i].Vm * (e.G/tap_ratio_square) - node[i].Vm * node[p].Vm * ((e.G/tap_ratio)*cos(node[i].Va-node[p].Va) + (-e.B/tap_ratio)*sin(node[i].Va-node[p].Va)))) * e.Ri_eP);
    			}
    			else{
    				tmp_H_r_P += e.BIJ * (e.M_P_TLPF - (node[i].Vm * node[i].Vm * e.G - node[i].Vm * node[p].Vm * ((e.G/tap_ratio)*cos(node[i].Va-node[p].Va) + (-e.B/tap_ratio)*sin(node[i].Va-node[p].Va)))) * e.Ri_eP;
    				tmp_H_r_Q += (e.B/tap_ratio - e.hB) * (e.M_Q_TLPF - (- node[i].Vm * node[i].Vm * (-e.B + 0.5*e.hB) - node[i].Vm * node[p].Vm * ((e.G/tap_ratio)*sin(node[i].Va-node[p].Va) - (-e.B/tap_ratio)*cos(node[i].Va-node[p].Va)))) * e.Ri_eQ;
    			//if(i==0&&iter==0) printf("(%d<-(%d,%lf))\n",i,p,e.BIJ * (e.M_P_TLPF - (node[i].Vm * node[i].Vm * e.G - node[i].Vm * node[p].Vm * ((e.G/tap_ratio)*cos(node[i].Va-node[p].Va) + (-e.B/tap_ratio)*sin(node[i].Va-node[p].Va)))) * e.Ri_eP);
    			}

    			if((-e.K) == 0){
    			    tmp_H_r_P += (-1) * e.BIJ * (e.M_P_TLPF_reverse - (node[p].Vm * node[p].Vm * e.G - node[p].Vm * node[i].Vm * (e.G*cos(node[p].Va-node[i].Va) + (-e.B)*sin(node[p].Va-node[i].Va)))) * e.Ri_eP_reverse;
    			    tmp_H_r_Q += (-1) * e.B * (e.M_Q_TLPF_reverse - (- node[p].Vm * node[p].Vm * (-e.B + 0.5*e.hB) - node[p].Vm * node[i].Vm * (e.G*sin(node[p].Va-node[i].Va) - (-e.B)*cos(node[p].Va-node[i].Va)))) * e.Ri_eQ_reverse;
    			//if(i==0&&iter==0) printf("(%d<-(%d,%lf))\n",i,p,(-1) * e.BIJ * (e.M_P_TLPF_reverse - (node[p].Vm * node[p].Vm * e.G - node[p].Vm * node[i].Vm * (e.G*cos(node[p].Va-node[i].Va) + (-e.B)*sin(node[p].Va-node[i].Va)))) * e.Ri_eP_reverse);
    			}
    			else if((-e.K) > 0){
    			    tmp_H_r_P += (-1) * e.BIJ * (e.M_P_TLPF_reverse - (node[p].Vm * node[p].Vm * (e.G/tap_ratio_square) - node[p].Vm * node[i].Vm * ((e.G/tap_ratio)*cos(node[p].Va-node[i].Va) + (-e.B/tap_ratio)*sin(node[p].Va-node[i].Va)))) * e.Ri_eP_reverse;
    			    tmp_H_r_Q += (-1) * (e.B/tap_ratio) * (e.M_Q_TLPF_reverse - (- node[p].Vm * node[p].Vm * (-e.B + 0.5*e.hB) / tap_ratio_square - node[p].Vm * node[i].Vm * ((e.G/tap_ratio)*sin(node[p].Va-node[i].Va) - (-e.B/tap_ratio)*cos(node[p].Va-node[i].Va)))) * e.Ri_eQ_reverse;
    			//if(i==0&&iter==0) printf("(%d<-(%d,%lf))\n",i,p,(-1) * e.BIJ * (e.M_P_TLPF_reverse - (node[p].Vm * node[p].Vm * (e.G/tap_ratio_square) - node[p].Vm * node[i].Vm * ((e.G/tap_ratio)*cos(node[p].Va-node[i].Va) + (-e.B/tap_ratio)*sin(node[p].Va-node[i].Va)))) * e.Ri_eP_reverse);
    			}
    			else{
    			    tmp_H_r_P += (-1) * e.BIJ * (e.M_P_TLPF_reverse - (node[p].Vm * node[p].Vm * e.G - node[p].Vm * node[i].Vm * ((e.G/tap_ratio)*cos(node[p].Va-node[i].Va) + (-e.B/tap_ratio)*sin(node[p].Va-node[i].Va)))) * e.Ri_eP_reverse;
    			    tmp_H_r_Q += (-1) * (e.B/tap_ratio) * (e.M_Q_TLPF_reverse - (- node[p].Vm * node[p].Vm * (-e.B + 0.5*e.hB) - node[p].Vm * node[i].Vm * ((e.G/tap_ratio)*sin(node[p].Va-node[i].Va) - (-e.B/tap_ratio)*cos(node[p].Va-node[i].Va)))) * e.Ri_eQ_reverse;
    			//if(i==0&&iter==0) printf("(%d<-(%d,%lf))\n",i,p,(-1) * e.BIJ * (e.M_P_TLPF_reverse - (node[p].Vm * node[p].Vm * e.G - node[p].Vm * node[i].Vm * ((e.G/tap_ratio)*cos(node[p].Va-node[i].Va) + (-e.B/tap_ratio)*sin(node[p].Va-node[i].Va)))) * e.Ri_eP_reverse);
    			}

    			//debug
    			/*if(e.K == 0){
    			    debug_H_r_P[i] += e.BIJ * (e.M_P_TLPF - (node[i].Vm * node[i].Vm * e.G - node[i].Vm * node[p].Vm * (e.G*cos(node[i].Va-node[p].Va) + (-e.B)*sin(node[i].Va-node[p].Va)))) * e.Ri_eP;
    			    debug_H_r_Q[i] += (e.B - e.hB) * (e.M_Q_TLPF - (- node[i].Vm * node[i].Vm * (-e.B + 0.5*e.hB) - node[i].Vm * node[p].Vm * (e.G*sin(node[i].Va-node[p].Va) - (-e.B)*cos(node[i].Va-node[p].Va)))) * e.Ri_eQ;
    			    debug_H_r_P[p] += (-1) * e.BIJ * (e.M_P_TLPF - (node[i].Vm * node[i].Vm * e.G - node[i].Vm * node[p].Vm * (e.G*cos(node[i].Va-node[p].Va) + (-e.B)*sin(node[i].Va-node[p].Va)))) * e.Ri_eP;
    			    debug_H_r_Q[p] += (-1) * e.B * (e.M_Q_TLPF - (- node[i].Vm * node[i].Vm * (-e.B + 0.5*e.hB) - node[i].Vm * node[p].Vm * (e.G*sin(node[i].Va-node[p].Va) - (-e.B)*cos(node[i].Va-node[p].Va)))) * e.Ri_eQ;
    			//if(i==0&&iter==0) printf("debug:(%d<-(%d,%lf))\n",i,p,e.BIJ * (e.M_P_TLPF - (node[i].Vm * node[i].Vm * e.G - node[i].Vm * node[p].Vm * (e.G*cos(node[i].Va-node[p].Va) + (-e.B)*sin(node[i].Va-node[p].Va)))) * e.Ri_eP);
    			//if(p==0&&iter==0) printf("debug:(%d<-(%d,%lf))\n",p,i,(-1) * e.BIJ * (e.M_P_TLPF - (node[i].Vm * node[i].Vm * e.G - node[i].Vm * node[p].Vm * (e.G*cos(node[i].Va-node[p].Va) + (-e.B)*sin(node[i].Va-node[p].Va)))) * e.Ri_eP);
    			}
    			else if(e.K > 0){
    			    debug_H_r_P[i] += e.BIJ * (e.M_P_TLPF - (node[i].Vm * node[i].Vm * (e.G/tap_ratio_square) - node[i].Vm * node[p].Vm * ((e.G/tap_ratio)*cos(node[i].Va-node[p].Va) + (-e.B/tap_ratio)*sin(node[i].Va-node[p].Va)))) * e.Ri_eP;
    			    debug_H_r_Q[i] += (e.B/tap_ratio - e.hB) * (e.M_Q_TLPF - (- node[i].Vm * node[i].Vm * (-e.B + 0.5*e.hB) / tap_ratio_square - node[i].Vm * node[p].Vm * ((e.G/tap_ratio)*sin(node[i].Va-node[p].Va) - (-e.B/tap_ratio)*cos(node[i].Va-node[p].Va)))) * e.Ri_eQ;
    			    debug_H_r_P[p] += (-1) * e.BIJ * (e.M_P_TLPF - (node[i].Vm * node[i].Vm * (e.G/tap_ratio_square) - node[i].Vm * node[p].Vm * ((e.G/tap_ratio)*cos(node[i].Va-node[p].Va) + (-e.B/tap_ratio)*sin(node[i].Va-node[p].Va)))) * e.Ri_eP;
    			    debug_H_r_Q[p] += (-1) * (e.B/tap_ratio) * (e.M_Q_TLPF - (- node[i].Vm * node[i].Vm * (-e.B + 0.5*e.hB) / tap_ratio_square - node[i].Vm * node[p].Vm * ((e.G/tap_ratio)*sin(node[i].Va-node[p].Va) - (-e.B/tap_ratio)*cos(node[i].Va-node[p].Va)))) * e.Ri_eQ;
    			//if(i==0&&iter==0) printf("debug:(%d<-(%d,%lf))\n",i,p,e.BIJ * (e.M_P_TLPF - (node[i].Vm * node[i].Vm * (e.G/tap_ratio_square) - node[i].Vm * node[p].Vm * ((e.G/tap_ratio)*cos(node[i].Va-node[p].Va) + (-e.B/tap_ratio)*sin(node[i].Va-node[p].Va)))) * e.Ri_eP);
    			//if(p==0&&iter==0) printf("debug:(%d<-(%d,%lf))\n",p,i,(-1) * e.BIJ * (e.M_P_TLPF - (node[i].Vm * node[i].Vm * (e.G/tap_ratio_square) - node[i].Vm * node[p].Vm * ((e.G/tap_ratio)*cos(node[i].Va-node[p].Va) + (-e.B/tap_ratio)*sin(node[i].Va-node[p].Va)))) * e.Ri_eP);
    			}
    			else{
    			    debug_H_r_P[i] += e.BIJ * (e.M_P_TLPF - (node[i].Vm * node[i].Vm * e.G - node[i].Vm * node[p].Vm * ((e.G/tap_ratio)*cos(node[i].Va-node[p].Va) + (-e.B/tap_ratio)*sin(node[i].Va-node[p].Va)))) * e.Ri_eP;
    			    debug_H_r_Q[i] += (e.B/tap_ratio - e.hB) * (e.M_Q_TLPF - (- node[i].Vm * node[i].Vm * (-e.B + 0.5*e.hB) - node[i].Vm * node[p].Vm * ((e.G/tap_ratio)*sin(node[i].Va-node[p].Va) - (-e.B/tap_ratio)*cos(node[i].Va-node[p].Va)))) * e.Ri_eQ;
    			    debug_H_r_P[p] += (-1) * e.BIJ * (e.M_P_TLPF - (node[i].Vm * node[i].Vm * e.G - node[i].Vm * node[p].Vm * ((e.G/tap_ratio)*cos(node[i].Va-node[p].Va) + (-e.B/tap_ratio)*sin(node[i].Va-node[p].Va)))) * e.Ri_eP;
    			    debug_H_r_Q[p] += (-1) * (e.B/tap_ratio) * (e.M_Q_TLPF - (- node[i].Vm * node[i].Vm * (-e.B + 0.5*e.hB) - node[i].Vm * node[p].Vm * ((e.G/tap_ratio)*sin(node[i].Va-node[p].Va) - (-e.B/tap_ratio)*cos(node[i].Va-node[p].Va)))) * e.Ri_eQ;
    			//if(i==0&&iter==0) printf("debug:(%d<-(%d,%lf))\n",i,p,e.BIJ * (e.M_P_TLPF - (node[i].Vm * node[i].Vm * e.G - node[i].Vm * node[p].Vm * ((e.G/tap_ratio)*cos(node[i].Va-node[p].Va) + (-e.B/tap_ratio)*sin(node[i].Va-node[p].Va)))) * e.Ri_eP);
    			//if(p==0&&iter==0) printf("debug:(%d<-(%d,%lf))\n",p,i,(-1) * e.BIJ * (e.M_P_TLPF - (node[i].Vm * node[i].Vm * e.G - node[i].Vm * node[p].Vm * ((e.G/tap_ratio)*cos(node[i].Va-node[p].Va) + (-e.B/tap_ratio)*sin(node[i].Va-node[p].Va)))) * e.Ri_eP);
    			}*/
    		}
    		deltaP = node[i].P - (deltaP + node[i].Vm*node[i].Vm*node[i].sumG);
    		tmp_H_r_P += (-1) * node[i].sumBi * deltaP * node[i].Ri_vP;
    		if(node[i].type!=3){
    			H_r_P[i] = tmp_H_r_P;
    		}
    		else{
    			H_r_P[i]=0;
    		}

    		deltaQ = node[i].Q - (deltaQ - node[i].Vm*node[i].Vm*node[i].sumB);
    		double deltaVm = node[i].M_Vm - node[i].Vm;
    		tmp_H_r_Q += (-1) * (2*node[i].sumB + node[i].sumBedge)* deltaQ * node[i].Ri_vQ + deltaVm * node[i].Ri_V;
    		H_r_Q[i] = tmp_H_r_Q;

    		//debug
    		//debug_H_r_P[i]+=(-1) * node[i].sumBi * deltaP * node[i].Ri_vP;
    		//debug_H_r_Q[i]+=(-1) * (2*node[i].sumB + node[i].sumBedge)* deltaQ * node[i].Ri_vQ + deltaVm * node[i].Ri_V;

    		//if(i<10&&iter==0) printf("deltaP=%lf,deltaQ=%lf\n",deltaP,deltaQ);
    	}
    	gettimeofday(&finish, 0);//计时结束
    	cost_time = 1000000 * (finish.tv_sec - start.tv_sec) + finish.tv_usec - start.tv_usec;
    	printf("rhs:cost the time is : %lf s.\n",cost_time/(double)1000000);

    	/*for(int i=0;i<node_num;i++){
    		if(node[i].type==3) debug_H_r_P[i]=0;
    		if(fabs(debug_H_r_P[i]-H_r_P[i])>0.0001||fabs(debug_H_r_Q[i]-H_r_Q[i])>0.0001){
    			if(i<10) printf("debug_H_r_P[%d]=%lf,H_r_P[%d]=%lf,debug_H_r_Q[%d]=%lf,H_r_Q[%d]=%lf\n",
    					i,debug_H_r_P[i],i,H_r_P[i],i,debug_H_r_Q[i],i,H_r_Q[i]);
    		}
    		//H_r_P[i]=debug_H_r_P[i];H_r_Q[i]=debug_H_r_Q[i];
    	}*/



    if(iter%2==0){
    	//for(int i=0;i<10;i++) printf("(%d,%lf,%lf)\n",i,H_r_P[i],H_r_Q[i]);

    	/*for(int i=slackbus;i<node_num;i++){
    		H_r_P[i]=H_r_P[i+1];
    	}*/

    	CKTSO_Solve(nicslu_p, H_r_P,H_r_P,false,false);

    	/*for(int i=node_num-1;i>slackbus;i--){
    	    H_r_P[i]=H_r_P[i-1];
    	}
    	H_r_P[slackbus]=0;*/

    	//for(int i=0;i<10;i++) printf("(%d,%lf,%lf)\n",i,H_r_P[i],H_r_Q[i]);
    	max_dVa=0;
    	for(int i=0;i<node_num;i++){
    		if(fabs(H_r_P[i])>max_dVa){
    	    	max_dVa=fabs(H_r_P[i]);
    	    }
    		node[i].Va+=H_r_P[i];
    	}
    	//printf("krnl:max_dVa=%lf,max_dVm=%lf\n",max_dVa,max_dVm);
    }
    else{
    	//for(int i=0;i<10;i++) printf("(%d,%lf,%lf)\n",i,H_r_P[i],H_r_Q[i]);
    	CKTSO_Solve(nicslu_pp, H_r_Q,H_r_Q,false,false );
    	//for(int i=0;i<10;i++) printf("(%d,%lf,%lf)\n",i,H_r_P[i],H_r_Q[i]);
    	max_dVm=0;
    	for(int i=0;i<node_num;i++){
    		if(fabs(H_r_Q[i])>max_dVm){
    	        max_dVm=fabs(H_r_Q[i]);
    	    }
    		node[i].Vm+=H_r_Q[i];
    	    //if(node[i].Vm>node[i].M_Vm) node[i].Vm=node[i].M_Vm;
    	    //if(node[i].Vm<-node[i].M_Vm) node[i].Vm=-node[i].M_Vm;
    	}
    	printf("krnl:max_dVa=%lf,max_dVm=%lf\n",max_dVa,max_dVm);
    }
    }

    for(int i=0;i<10;i++) printf("(%d,%lf,%lf,%lf)\n",i,node[i].Va,node[i].Vm,node[i].M_Vm);
    cout << "Hello world!" << endl;
    return 0;
}
