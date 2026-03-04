#include "pf2.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include <algorithm>
#include <iostream>
#include <sys/time.h>
#include <pthread.h>

struct Ret_Ybus pf_Ybus_seq(struct Graph_PF_IN graph_in){
	int node_num=graph_in.node_num;
	int edge_num=graph_in.edge_num;
	//struct Sort_Id_Vertex* vertex;
	struct Sort_Id_Edge_All *edge_all;
	double *sumG,*sumB,*sumBi;

	double *Bp_x,*Bpp_x;
	unsigned *Bp_i,*Bpp_i;
	unsigned *Bp_p,*Bpp_p;
	struct Graph_PF graph;

	//vertex=(struct Sort_Id_Vertex*)malloc(graph_in.node_num*sizeof(struct Sort_Id_Vertex));
	edge_all=(struct Sort_Id_Edge_All*)malloc((graph_in.edge_num+graph_in.node_num)*sizeof(struct Sort_Id_Edge_All));
	sumG=(double*)malloc(graph_in.node_num*sizeof(double));
	sumB=(double*)malloc(graph_in.node_num*sizeof(double));
	sumBi=(double*)malloc(graph_in.node_num*sizeof(double));

	Bp_x=(double*)malloc((graph_in.edge_num+graph_in.node_num)*sizeof(double));
	Bpp_x=(double*)malloc((graph_in.edge_num+graph_in.node_num)*sizeof(double));
	Bp_i=(unsigned*)malloc((graph_in.edge_num+graph_in.node_num)*sizeof(unsigned));
	Bpp_i=(unsigned*)malloc((graph_in.edge_num+graph_in.node_num)*sizeof(unsigned));
	Bp_p=(unsigned*)malloc((graph_in.node_num+1)*sizeof(unsigned));
	Bpp_p=(unsigned*)malloc((graph_in.node_num+1)*sizeof(unsigned));

	graph.node=(struct Node_PF *)malloc(graph_in.node_num*sizeof(struct Node_PF));
	graph.edge=(struct Edge_PF *)malloc((graph_in.edge_num+graph_in.node_num)*sizeof(struct Edge_PF));
	graph.off=(unsigned *)malloc((node_num+1)*sizeof(unsigned));
	graph.node_num=graph_in.node_num;
	graph.edge_num=graph_in.edge_num+graph_in.node_num;

	for(int v=0;v<node_num;v++){
		sumG[v]=0;
		sumB[v]=0;
		sumBi[v]=0;
		Bp_p[v+1]=1;
		Bpp_p[v+1]=1;
		double edge_all_eG,edge_all_eB;
		int edge_all_ei;
		double edge_all_Bp_x,edge_all_Bpp_x;
		int diag_off=graph_in.off[v];
		for(int j=graph_in.off[v];j<graph_in.off[v+1];j++){
			struct Edge_PF_In e= graph_in.edge[j];
			int t=e.to;
			edge_all_ei=t;
			if(graph_in.node[v].type!=3){
				Bp_p[v+1]+=1;
				edge_all_Bp_x=e.BIJ;//if(e.BIJ==0) printf("~~~~~~error1~~~~\n");
			}
			else{
				edge_all_Bp_x=0;
			}
			if(e.K==0){
				sumG[v] += e.G;
				sumB[v] += -1*e.B + 0.5*e.hB;
				sumBi[v] += -1*e.BIJ;
				edge_all_eG=-e.G;
				edge_all_eB=e.B;

				if((graph_in.node[v].type==0||graph_in.node[v].type==1)
						&&(graph_in.node[t].type==0||graph_in.node[t].type==1)){
					Bpp_p[v+1]+=1;
					edge_all_Bpp_x=e.B;
				}
				else{
					edge_all_Bpp_x=0;
				}
			}
			else if(e.K>0){
				double tap_ratio = e.K/e.Kcount;
				double tap_ratio_square = (e.K/e.Kcount)*(e.K/e.Kcount);
				sumG[v] += 1/(tap_ratio_square)*e.G;
				sumB[v] += 1/(tap_ratio_square)*(-1*e.B) + 0.5*e.hB; // sqrt
				sumBi[v] += -1*e.BIJ;
				edge_all_eG=-e.G/tap_ratio;
				edge_all_eB=e.B/tap_ratio;

				if((graph_in.node[v].type==0||graph_in.node[v].type==1)
						&&(graph_in.node[t].type==0||graph_in.node[t].type==1)){
					Bpp_p[v+1]+=1;
					edge_all_Bpp_x=e.B/tap_ratio;
				}
				else{
					edge_all_Bpp_x=0;
				}
			}
			else{
				double tap_ratio = fabs(e.K/e.Kcount);
				sumG[v] += e.G;
				sumB[v] += -1*e.B + 0.5*e.hB;
				sumBi[v] += -1*e.BIJ;
				edge_all_eG=-e.G/tap_ratio;
				edge_all_eB=e.B/tap_ratio;

				if((graph_in.node[v].type==0||graph_in.node[v].type==1)
					&&(graph_in.node[t].type==0||graph_in.node[t].type==1)){
					Bpp_p[v+1]+=1;
					edge_all_Bpp_x=e.B/tap_ratio;
				}
				else{
					edge_all_Bpp_x=0;
				}
			}
			if(t<v){
				edge_all[v+j]=(struct Sort_Id_Edge_All){edge_all_eG,edge_all_eB,edge_all_ei,edge_all_Bp_x,edge_all_Bpp_x};
				diag_off=j+1;
				//if(v+j<20) printf("~~(%d,%d),off=%d\n",v,t,v+j);
			}
			else{
				edge_all[v+j+1]=(struct Sort_Id_Edge_All){edge_all_eG,edge_all_eB,edge_all_ei,edge_all_Bp_x,edge_all_Bpp_x};
				//if(v+j+1<20) printf("~~(%d,%d),off=%d\n",v,t,v+j+1);
			}
		}
		sumG[v] += graph_in.node[v].G;
		sumB[v] += graph_in.node[v].B;
		if(graph_in.node[v].type==0||graph_in.node[v].type==1){
			edge_all[v+diag_off]=(struct Sort_Id_Edge_All){sumG[v],sumB[v],v,sumBi[v],sumB[v]};
			//if(sumBi[v]==0||sumB[v]==0) printf("~~~~~~error2~~~~sumBi[%d]=%lf,sumB[%d]=%lf\n",v,sumBi[v],v,sumB[v]);
		}
		else if(graph_in.node[v].type==2){
			edge_all[v+diag_off]=(struct Sort_Id_Edge_All){sumG[v],sumB[v],v,sumBi[v],1};
			//if(sumBi[v]==0) printf("~~~~~~error3~~~~\n");
		}
		else{
			edge_all[v+diag_off]=(struct Sort_Id_Edge_All){sumG[v],sumB[v],v,1,1};
		}
		//if(v+diag_off<20) printf("~~(%d,%d),off=%d\n",v,v,v+diag_off);

		graph.node[v].type=graph_in.node[v].type;
		graph.node[v].Vm=graph_in.node[v].Vm;
		graph.node[v].Va=graph_in.node[v].Va*PI/180;
		graph.node[v].P=graph_in.node[v].P;
		graph.node[v].Q=graph_in.node[v].Q;

		//graph.node[v].Vm=1;graph.node[v].Va=0;

	}


	for(int i=0,Bp_k=0,Bpp_k=0;i<(graph_in.edge_num+graph_in.node_num);i++){
		if(edge_all[i].Bp_x!=0){
			Bp_x[Bp_k]=edge_all[i].Bp_x;
			Bp_i[Bp_k]=edge_all[i].ei;
			Bp_k++;
		}
		if(edge_all[i].Bpp_x!=0){
			Bpp_x[Bpp_k]=edge_all[i].Bpp_x;
			Bpp_i[Bpp_k]=edge_all[i].ei;
			Bpp_k++;
		}
		//if(i==(graph_in.edge_num+graph_in.node_num)-1) printf("Bp_k:%d\nBpp_k:%d\n",Bp_k,Bpp_k);
		graph.edge[i].src=edge_all[i].ei;
		graph.edge[i].G=edge_all[i].eG;
		graph.edge[i].B=edge_all[i].eB;
	}
	Bp_p[0]=0;
	Bpp_p[0]=0;
	graph.off[0]=0;
	for(int i=1;i<graph_in.node_num+1;i++){
		Bp_p[i]+=Bp_p[i-1];
		Bpp_p[i]+=Bpp_p[i-1];
		graph.off[i]=graph_in.off[i]+i;
	}

	free(edge_all);
	free(sumG);
	free(sumB);
	free(sumBi);
	//free(graph_in.node);
	//free(graph_in.edge);
	//free(graph_in.off);

//debug
/*	printf("graph_in:node_num=%d\tedge_num=%d\n",graph.node_num,graph.edge_num);
	printf("node:\n");
	for(int i=0;i<graph.node_num;i++){
		if(i>=5&&graph.node_num-i>5) continue;
		printf("%d\t%d\t",i,graph.node[i].type);
		printf("%lf\t%lf\t%lf\t%lf\n",graph.node[i].Vm,graph.node[i].Va,graph.node[i].P,graph.node[i].Q);
	}
	printf("edge:\n");
	for(int i=0;i<graph.edge_num;i++){
		if(i>=5&&graph.edge_num-i>5) continue;
		printf("%d\t%d\t",i,graph.edge[i].src);
		printf("%lf\t%lf\n",graph.edge[i].G,graph.edge[i].B);
	}
	printf("off:\n");
	for(int i=0;i<graph.node_num+1;i++){
		if(i>=5&&graph.node_num+1-i>5) continue;
		printf("%d\t%d\n",i,graph.off[i]);
	}
*/

	struct Ret_Ybus ret_Ybus;
	ret_Ybus.graph=graph;
	ret_Ybus.Bp_x=Bp_x;
	ret_Ybus.Bp_i=Bp_i;
	ret_Ybus.Bp_p=Bp_p;
	ret_Ybus.Bpp_x=Bpp_x;
	ret_Ybus.Bpp_i=Bpp_i;
	ret_Ybus.Bpp_p=Bpp_p;
	return ret_Ybus;
}


#define YBUS_THREAD_NUM 8
struct Ybus_Para_Mt{
	struct Graph_PF_IN graph_in;
	struct Graph_PF graph;
	struct Sort_Id_Edge_All *edge_all;
	double *sumG,*sumB,*sumBi;
	unsigned *Bp_p,*Bpp_p;
	int local_id;
	int local_size;
};

void *pf_Ybus_mt(void *arg){
	struct Ybus_Para_Mt *para=(struct Ybus_Para_Mt*)arg;
	struct Graph_PF_IN graph_in=para->graph_in;
	struct Graph_PF graph=para->graph;
	struct Sort_Id_Edge_All *edge_all=para->edge_all;
	double *sumG=para->sumG;
	double *sumB=para->sumB;
	double *sumBi=para->sumBi;
	unsigned *Bp_p=para->Bp_p;
	unsigned *Bpp_p=para->Bpp_p;
	int local_id=para->local_id;
	int local_size=para->local_size;
	int node_num=para->graph_in.node_num;
	int begin_node=node_num*local_id/local_size;
	int end_node=node_num*(local_id+1)/local_size;

	for(int v=begin_node;v<end_node;v++){
		sumG[v]=0;
		sumB[v]=0;
		sumBi[v]=0;
		Bp_p[v+1]=1;
		Bpp_p[v+1]=1;
		double edge_all_eG,edge_all_eB;
		int edge_all_ei;
		double edge_all_Bp_x,edge_all_Bpp_x;
		int diag_off=graph_in.off[v];
		for(int j=graph_in.off[v];j<graph_in.off[v+1];j++){
			struct Edge_PF_In e= graph_in.edge[j];
			int t=e.to;
			edge_all_ei=t;
			if(graph_in.node[v].type!=3){
				Bp_p[v+1]+=1;
				edge_all_Bp_x=e.BIJ;//if(e.BIJ==0) printf("~~~~~~error1~~~~\n");
			}
			else{
				edge_all_Bp_x=0;
			}
			if(e.K==0){
				sumG[v] += e.G;
				sumB[v] += -1*e.B + 0.5*e.hB;
				sumBi[v] += -1*e.BIJ;
				edge_all_eG=-e.G;
				edge_all_eB=e.B;

				if((graph_in.node[v].type==0||graph_in.node[v].type==1)
						&&(graph_in.node[t].type==0||graph_in.node[t].type==1)){
					Bpp_p[v+1]+=1;
					edge_all_Bpp_x=e.B;
				}
				else{
					edge_all_Bpp_x=0;
				}
			}
			else if(e.K>0){
				double tap_ratio = e.K/e.Kcount;
				double tap_ratio_square = (e.K/e.Kcount)*(e.K/e.Kcount);
				sumG[v] += 1/(tap_ratio_square)*e.G;
				sumB[v] += 1/(tap_ratio_square)*(-1*e.B) + 0.5*e.hB; // sqrt
				sumBi[v] += -1*e.BIJ;
				edge_all_eG=-e.G/tap_ratio;
				edge_all_eB=e.B/tap_ratio;

				if((graph_in.node[v].type==0||graph_in.node[v].type==1)
						&&(graph_in.node[t].type==0||graph_in.node[t].type==1)){
					Bpp_p[v+1]+=1;
					edge_all_Bpp_x=e.B/tap_ratio;
				}
				else{
					edge_all_Bpp_x=0;
				}
			}
			else{
				double tap_ratio = fabs(e.K/e.Kcount);
				sumG[v] += e.G;
				sumB[v] += -1*e.B + 0.5*e.hB;
				sumBi[v] += -1*e.BIJ;
				edge_all_eG=-e.G/tap_ratio;
				edge_all_eB=e.B/tap_ratio;

				if((graph_in.node[v].type==0||graph_in.node[v].type==1)
					&&(graph_in.node[t].type==0||graph_in.node[t].type==1)){
					Bpp_p[v+1]+=1;
					edge_all_Bpp_x=e.B/tap_ratio;
				}
				else{
					edge_all_Bpp_x=0;
				}
			}
			if(t<v){
				edge_all[v+j]=(struct Sort_Id_Edge_All){edge_all_eG,edge_all_eB,edge_all_ei,edge_all_Bp_x,edge_all_Bpp_x};
				diag_off=j+1;
				//if(v+j<20) printf("~~(%d,%d),off=%d\n",v,t,v+j);
			}
			else{
				edge_all[v+j+1]=(struct Sort_Id_Edge_All){edge_all_eG,edge_all_eB,edge_all_ei,edge_all_Bp_x,edge_all_Bpp_x};
				//if(v+j+1<20) printf("~~(%d,%d),off=%d\n",v,t,v+j+1);
			}
		}
		sumG[v] += graph_in.node[v].G;
		sumB[v] += graph_in.node[v].B;
		if(graph_in.node[v].type==0||graph_in.node[v].type==1){
			edge_all[v+diag_off]=(struct Sort_Id_Edge_All){sumG[v],sumB[v],v,sumBi[v],sumB[v]};
			//if(sumBi[v]==0||sumB[v]==0) printf("~~~~~~error2~~~~sumBi[%d]=%lf,sumB[%d]=%lf\n",v,sumBi[v],v,sumB[v]);
		}
		else if(graph_in.node[v].type==2){
			edge_all[v+diag_off]=(struct Sort_Id_Edge_All){sumG[v],sumB[v],v,sumBi[v],1};
			//if(sumBi[v]==0) printf("~~~~~~error3~~~~\n");
		}
		else{
			edge_all[v+diag_off]=(struct Sort_Id_Edge_All){sumG[v],sumB[v],v,1,1};
		}
		//if(v+diag_off<20) printf("~~(%d,%d),off=%d\n",v,v,v+diag_off);

		graph.node[v].type=graph_in.node[v].type;
		graph.node[v].Vm=graph_in.node[v].Vm;
		graph.node[v].Va=graph_in.node[v].Va*PI/180;
		graph.node[v].P=graph_in.node[v].P;
		graph.node[v].Q=graph_in.node[v].Q;


		//graph.node[v].Vm=1;graph.node[v].Va=0;

		}
        return nullptr;
}



struct Ret_Ybus pf_Ybus_par(struct Graph_PF_IN graph_in){
	int node_num=graph_in.node_num;
	int edge_num=graph_in.edge_num;
	//struct Sort_Id_Vertex* vertex;
	struct Sort_Id_Edge_All *edge_all;
	double *sumG,*sumB,*sumBi;

	double *Bp_x,*Bpp_x;
	unsigned *Bp_i,*Bpp_i;
	unsigned *Bp_p,*Bpp_p;
	struct Graph_PF graph;

	//vertex=(struct Sort_Id_Vertex*)malloc(graph_in.node_num*sizeof(struct Sort_Id_Vertex));
	edge_all=(struct Sort_Id_Edge_All*)malloc((graph_in.edge_num+graph_in.node_num)*sizeof(struct Sort_Id_Edge_All));
	sumG=(double*)malloc(graph_in.node_num*sizeof(double));
	sumB=(double*)malloc(graph_in.node_num*sizeof(double));
	sumBi=(double*)malloc(graph_in.node_num*sizeof(double));

	Bp_x=(double*)malloc((graph_in.edge_num+graph_in.node_num)*sizeof(double));
	Bpp_x=(double*)malloc((graph_in.edge_num+graph_in.node_num)*sizeof(double));
	Bp_i=(unsigned*)malloc((graph_in.edge_num+graph_in.node_num)*sizeof(unsigned));
	Bpp_i=(unsigned*)malloc((graph_in.edge_num+graph_in.node_num)*sizeof(unsigned));
	Bp_p=(unsigned*)malloc((graph_in.node_num+1)*sizeof(unsigned));
	Bpp_p=(unsigned*)malloc((graph_in.node_num+1)*sizeof(unsigned));

	graph.node=(struct Node_PF *)malloc(graph_in.node_num*sizeof(struct Node_PF));
	graph.edge=(struct Edge_PF *)malloc((graph_in.edge_num+graph_in.node_num)*sizeof(struct Edge_PF));
	graph.off=(unsigned *)malloc((node_num+1)*sizeof(unsigned));
	graph.node_num=graph_in.node_num;
	graph.edge_num=graph_in.edge_num+graph_in.node_num;


	struct Ybus_Para_Mt Ybus_para_mt[16];
	pthread_t Ybus_p_arr[16];
	for(int tid=0;tid<YBUS_THREAD_NUM;tid++){
		Ybus_para_mt[tid].graph_in=graph_in;
		Ybus_para_mt[tid].graph=graph;
		Ybus_para_mt[tid].edge_all=edge_all;
		Ybus_para_mt[tid].sumG=sumG;
		Ybus_para_mt[tid].sumB=sumB;
		Ybus_para_mt[tid].sumBi=sumBi;
		Ybus_para_mt[tid].Bp_p=Bp_p;
		Ybus_para_mt[tid].Bpp_p=Bpp_p;
		Ybus_para_mt[tid].local_id=tid;
		Ybus_para_mt[tid].local_size=YBUS_THREAD_NUM;
	}
	for(int tid=0;tid<YBUS_THREAD_NUM;tid++){
		int err=pthread_create(&Ybus_p_arr[tid],NULL,pf_Ybus_mt,&Ybus_para_mt[tid]);
	}
	for(int tid=0;tid<YBUS_THREAD_NUM;tid++){
		pthread_join(Ybus_p_arr[tid],NULL);
	}



	/*for(int v=0;v<node_num;v++){
		sumG[v]=0;
		sumB[v]=0;
		sumBi[v]=0;
		Bp_p[v+1]=1;
		Bpp_p[v+1]=1;
		double edge_all_eG,edge_all_eB;
		int edge_all_ei;
		double edge_all_Bp_x,edge_all_Bpp_x;
		int diag_off=graph_in.off[v];
		for(int j=graph_in.off[v];j<graph_in.off[v+1];j++){
			struct Edge_PF_In e= graph_in.edge[j];
			int t=e.to;
			edge_all_ei=t;
			if(graph_in.node[v].type!=3){
				Bp_p[v+1]+=1;
				edge_all_Bp_x=e.BIJ;//if(e.BIJ==0) printf("~~~~~~error1~~~~\n");
			}
			else{
				edge_all_Bp_x=0;
			}
			if(e.K==0){
				sumG[v] += e.G;
				sumB[v] += -1*e.B + 0.5*e.hB;
				sumBi[v] += -1*e.BIJ;
				edge_all_eG=-e.G;
				edge_all_eB=e.B;

				if((graph_in.node[v].type==0||graph_in.node[v].type==1)
						&&(graph_in.node[t].type==0||graph_in.node[t].type==1)){
					Bpp_p[v+1]+=1;
					edge_all_Bpp_x=e.B;
				}
				else{
					edge_all_Bpp_x=0;
				}
			}
			else if(e.K>0){
				double tap_ratio = e.K/e.Kcount;
				double tap_ratio_square = (e.K/e.Kcount)*(e.K/e.Kcount);
				sumG[v] += 1/(tap_ratio_square)*e.G;
				sumB[v] += 1/(tap_ratio_square)*(-1*e.B) + 0.5*e.hB; // sqrt
				sumBi[v] += -1*e.BIJ;
				edge_all_eG=-e.G/tap_ratio;
				edge_all_eB=e.B/tap_ratio;

				if((graph_in.node[v].type==0||graph_in.node[v].type==1)
						&&(graph_in.node[t].type==0||graph_in.node[t].type==1)){
					Bpp_p[v+1]+=1;
					edge_all_Bpp_x=e.B/tap_ratio;
				}
				else{
					edge_all_Bpp_x=0;
				}
			}
			else{
				double tap_ratio = fabs(e.K/e.Kcount);
				sumG[v] += e.G;
				sumB[v] += -1*e.B + 0.5*e.hB;
				sumBi[v] += -1*e.BIJ;
				edge_all_eG=-e.G/tap_ratio;
				edge_all_eB=e.B/tap_ratio;

				if((graph_in.node[v].type==0||graph_in.node[v].type==1)
					&&(graph_in.node[t].type==0||graph_in.node[t].type==1)){
					Bpp_p[v+1]+=1;
					edge_all_Bpp_x=e.B/tap_ratio;
				}
				else{
					edge_all_Bpp_x=0;
				}
			}
			if(t<v){
				edge_all[v+j]=(struct Sort_Id_Edge_All){edge_all_eG,edge_all_eB,edge_all_ei,edge_all_Bp_x,edge_all_Bpp_x};
				diag_off=j+1;
				//if(v+j<20) printf("~~(%d,%d),off=%d\n",v,t,v+j);
			}
			else{
				edge_all[v+j+1]=(struct Sort_Id_Edge_All){edge_all_eG,edge_all_eB,edge_all_ei,edge_all_Bp_x,edge_all_Bpp_x};
				//if(v+j+1<20) printf("~~(%d,%d),off=%d\n",v,t,v+j+1);
			}
		}
		sumG[v] += graph_in.node[v].G;
		sumB[v] += graph_in.node[v].B;
		if(graph_in.node[v].type==0||graph_in.node[v].type==1){
			edge_all[v+diag_off]=(struct Sort_Id_Edge_All){sumG[v],sumB[v],v,sumBi[v],sumB[v]};
			//if(sumBi[v]==0||sumB[v]==0) printf("~~~~~~error2~~~~sumBi[%d]=%lf,sumB[%d]=%lf\n",v,sumBi[v],v,sumB[v]);
		}
		else if(graph_in.node[v].type==2){
			edge_all[v+diag_off]=(struct Sort_Id_Edge_All){sumG[v],sumB[v],v,sumBi[v],1};
			//if(sumBi[v]==0) printf("~~~~~~error3~~~~\n");
		}
		else{
			edge_all[v+diag_off]=(struct Sort_Id_Edge_All){sumG[v],sumB[v],v,1,1};
		}
		//if(v+diag_off<20) printf("~~(%d,%d),off=%d\n",v,v,v+diag_off);

		graph.node[v].type=graph_in.node[v].type;
		graph.node[v].Vm=graph_in.node[v].Vm;
		graph.node[v].Va=graph_in.node[v].Va*PI/180;
		graph.node[v].P=graph_in.node[v].P;
		graph.node[v].Q=graph_in.node[v].Q;


		//graph.node[v].Vm=1;graph.node[v].Va=0;

	}*/


	for(int i=0,Bp_k=0,Bpp_k=0;i<(graph_in.edge_num+graph_in.node_num);i++){
		if(edge_all[i].Bp_x!=0){
			Bp_x[Bp_k]=edge_all[i].Bp_x;
			Bp_i[Bp_k]=edge_all[i].ei;
			Bp_k++;
		}
		if(edge_all[i].Bpp_x!=0){
			Bpp_x[Bpp_k]=edge_all[i].Bpp_x;
			Bpp_i[Bpp_k]=edge_all[i].ei;
			Bpp_k++;
		}
		//if(i==(graph_in.edge_num+graph_in.node_num)-1) printf("Bp_k:%d\nBpp_k:%d\n",Bp_k,Bpp_k);
		graph.edge[i].src=edge_all[i].ei;
		graph.edge[i].G=edge_all[i].eG;
		graph.edge[i].B=edge_all[i].eB;
	}
	Bp_p[0]=0;
	Bpp_p[0]=0;
	graph.off[0]=0;
	for(int i=1;i<graph_in.node_num+1;i++){
		Bp_p[i]+=Bp_p[i-1];
		Bpp_p[i]+=Bpp_p[i-1];
		graph.off[i]=graph_in.off[i]+i;
	}

	free(edge_all);
	free(sumG);
	free(sumB);
	free(sumBi);
	//free(graph_in.node);
	//free(graph_in.edge);
	//free(graph_in.off);

	struct Ret_Ybus ret_Ybus;
	ret_Ybus.graph=graph;
	ret_Ybus.Bp_x=Bp_x;
	ret_Ybus.Bp_i=Bp_i;
	ret_Ybus.Bp_p=Bp_p;
	ret_Ybus.Bpp_x=Bpp_x;
	ret_Ybus.Bpp_i=Bpp_i;
	ret_Ybus.Bpp_p=Bpp_p;
	return ret_Ybus;
}
