#include "se.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include <algorithm>
#include <iostream>
#include <sys/time.h>
#include <pthread.h>

#include <vector>
#include <queue>
#include <fstream>
#include <sstream>
#include <string>

using namespace std;
struct Color{
	int r,g,b;
};
struct Color color[3]={{236,81,72},{236,81,72},{236,81,72}};
void power_show(struct Graph_SE graph,struct Graph_SE_IN graph_in){
	char outfile_node[128]="/home/qch/se_src/src/se_out/se_node.csv";
	FILE *fp_node=fopen(outfile_node,"w");
	char outfile_edge[128]="/home/qch/se_src/src/se_out/se_edge.csv";
	FILE *fp_edge=fopen(outfile_edge,"w");

	int node_num=graph.node_num;
	int edge_num=graph.edge_num;
	struct Node_SE* node=graph.node;
	struct Edge_SE* edge=graph.edge;
	unsigned* off=graph.off;

	bool *node_is_show=(bool *)malloc(node_num*sizeof(bool));
	int *node_new_id=(int *)malloc(node_num*sizeof(int));
	double *x=(double *)malloc(node_num*sizeof(double));
	double *y=(double *)malloc(node_num*sizeof(double));
	queue<int> vec_node;
	int show_num=300;
	int slackbus=0;
	for(int i=0;i<node_num;i++){
		if(node[i].type==3){
			node_is_show[i]=true;
			vec_node.push(i);
			x[i]=0;
			y[i]=0;
			slackbus=i;
		}
		else{
			node_is_show[i]=false;
		}
	}

	int show_count=1;
	while(show_count<show_num&&!vec_node.empty()){
		int i=vec_node.front();
		vec_node.pop();
		for(int j=off[i];j<off[i+1];j++){
			if(node_is_show[edge[j].src]==false){
				node_is_show[edge[j].src]=true;

				int locate=rand()%4;
				if(locate==0){
					x[edge[j].src]=15+(rand()%40000)/1000.0+x[i];
					y[edge[j].src]=10+(rand()%30000)/1000.0+y[i];
				}
				else if(locate==1){
					x[edge[j].src]=-(15+(rand()%40000)/1000.0)+x[i];
					y[edge[j].src]=10+(rand()%30000)/1000.0+y[i];
				}
				else if(locate==2){
					x[edge[j].src]=-(15+(rand()%40000)/1000.0)+x[i];
					y[edge[j].src]=-(10+(rand()%30000)/1000.0)+y[i];
				}
				else{
					x[edge[j].src]=15+(rand()%40000)/1000.0+x[i];
					y[edge[j].src]=-(10+(rand()%30000)/1000.0)+y[i];
				}

				if(x[edge[j].src]>390) x[edge[j].src]=390;
				if(x[edge[j].src]<-390) x[edge[j].src]=-390;
				if(y[edge[j].src]>290) y[edge[j].src]=290;
				if(y[edge[j].src]<-290) y[edge[j].src]=-290;
				vec_node.push(edge[j].src);
				show_count++;
			}
		}
	}
	//printf("slackbus=%d,show_count=%d\n",slackbus,show_count);

	for(int i=0,k=0;i<node_num;i++){
		if(node_is_show[i]==true){
			node_new_id[i]=k;
			k++;
		}
		else{
			node_new_id[i]=2333333;
		}
	}

	for(int i=0,node_count=0;i<node_num;i++){
		if(node_is_show[i]==false) continue;

		int type_value;
		if(node[i].type==0||node[i].type==1) type_value=0;//PQ
		else if(node[i].type==2) type_value=1;//PV
		else type_value=2;//slave

		//fprintf(fp,"          <attvalue for=\"node_info_class\" value=\"\\n量测电压幅值:%f\\n量测电压相位:%f\\n量测有功功率:%f\\n量测无功功率:%f\\n电压幅值权重:%f\\n有功功率权重:%f\\n无功功率权重:%f\\n收敛电压幅值:%f\\n收敛电压相位:%f\"></attvalue>\n",
						//graph_in.node[i].Vm,graph_in.node[i].Va*PI/180,graph_in.node[i].P,graph_in.node[i].Q,graph_in.node[i].Ri_V,graph_in.node[i].Ri_vP,graph_in.node[i].Ri_vQ,node[i].Vm,node[i].Va);
		if(node_count) fprintf(fp_node,"\n");
		fprintf(fp_node,"%d,%d,%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f",node_new_id[i],i,type_value,x[i],y[i],graph_in.node[i].Vm,graph_in.node[i].Va*PI/180,graph_in.node[i].P,graph_in.node[i].Q,graph_in.node[i].Ri_V,graph_in.node[i].Ri_vP,graph_in.node[i].Ri_vQ,node[i].Vm,node[i].Va);
		node_count++;
	}

	for(int i=0,edge_count=0;i<node_num;i++){
		for(int j=off[i];j<off[i+1];j++){
			if(node_is_show[i]==false||node_is_show[edge[j].src]==false) continue;
			if(edge_count) fprintf(fp_edge,"\n");
			fprintf(fp_edge,"%d,%d,%d,%d,%d,%f,%f,%f,%f,%d,%f,%f,%f,%f,%f,%f,%f,%f",
								edge_count,node_new_id[edge[j].src],node_new_id[i],edge[j].src,i,
								edge[j].G,edge[j].B,edge[j].hB,edge[j].K,edge[j].Kcount,edge[j].M_P_TLPF,edge[j].M_Q_TLPF,edge[j].Ri_eP,edge[j].Ri_eQ,edge[j].M_P_TLPF_reverse,edge[j].M_Q_TLPF_reverse,edge[j].Ri_eP_reverse,edge[j].Ri_eQ_reverse);
			//fprintf(fp,"      <edge id=\"%d\" source=\"%d\" target=\"%d\" weight=\"\\n电导:%f\\n电纳:%f\\n补偿电纳:%f\\n变压器变比:%f\\n支路数:%d\\n量测有功功率:%f\\n量测无功功率:%f\\n线路有功权重:%f\\n线路无功权重:%f\">\n",
					//edge_count,node_new_id[edge[j].src],node_new_id[i],edge[j].G,edge[j].B,edge[j].hB,edge[j].K,edge[j].Kcount,edge[j].M_P_TLPF,edge[j].M_Q_TLPF,edge[j].Ri_eP,edge[j].Ri_eQ);
			edge_count++;
		}
	}

	fflush(fp_node);
	fflush(fp_edge);
	fclose(fp_node);
	fclose(fp_edge);
}

void power_show2(struct Graph_SE graph,double precision,int iter,double excusion_time,double *total_time,double total_time_count){
	/*char outfile[128]="/home/zyh/power_data/se_out/se.json";
	FILE *fp=fopen(outfile,"w");

	fprintf(fp,"\\n节点数:%d\\n有向边数:%d\\n精度:%f\\n迭代次数:%d\\n总执行时间:%fms",
			graph.node_num,graph.edge_num,precision,iter,excusion_time);

	fclose(fp);*/

	char outfile2[128]="/home/qch/se_src/src/se_out/se_time.csv";
	FILE *fp2=fopen(outfile2,"w");

	double excusion_time2=0;
	for(int i=0;i<total_time_count;i++){
		fprintf(fp2,"%f",total_time[i]);
		if(i!=total_time_count-1){
			fprintf(fp2,",");
		}
		excusion_time2+=total_time[i];
	}
	fflush(fp2);
	fclose(fp2);

	//printf("excusion_time=%f,excusion_time2=%f\n",excusion_time,excusion_time2);
}
