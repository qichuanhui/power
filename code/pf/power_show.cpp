#include "pf.h"

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
void power_show(struct Graph_PF graph,struct Graph_PF_IN graph_in){
	char outfile[128]="/home/qch/src/pf_out/Case10790.gexf";
	FILE *fp=fopen(outfile,"w");
	fprintf(fp,"<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
	fprintf(fp,"<gexf xmlns=\"http://www.gexf.net/1.2draft\" version=\"1.2\" xmlns:viz=\"http://www.gexf.net/1.2draft/viz\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:schemaLocation=\"http://www.gexf.net/1.2draft http://www.gexf.net/1.2draft/gexf.xsd\">\n");
	fprintf(fp,"  <meta lastmodifieddate=\"2014-01-30\">\n");
	fprintf(fp,"    <creator>Gephi 0.8.1</creator>\n");
	fprintf(fp,"    <description></description>\n");
	fprintf(fp,"  </meta>\n");
	fprintf(fp,"  <graph defaultedgetype=\"directed\" idtype=\"static\" mode=\"static\">\n");
	fprintf(fp,"    <attributes class=\"node\" mode=\"static\">\n");
	fprintf(fp,"      <attribute id=\"modularity_class\" title=\"Modularity Class\" type=\"integer\"></attribute>\n");
	//fprintf(fp,"      <attribute id=\"node_Vm_class\" title=\"Node_Vm Class\" type=\"float\"></attribute>\n");
	//fprintf(fp,"      <attribute id=\"node_Va_class\" title=\"Node_Va Class\" type=\"float\"></attribute>\n");
	fprintf(fp,"      <attribute id=\"node_info_class\" title=\"Node_info Class\" type=\"string\"></attribute>\n");
	fprintf(fp,"    </attributes>\n");

	//fprintf(fp,"    <attributes class=\"edge\" mode=\"static\">\n");
	//fprintf(fp,"      <attribute id=\"edge_info_class\" title=\"Edge_info Class\" type=\"string\"></attribute>\n");
	//fprintf(fp,"    </attributes>\n");

	fprintf(fp,"    <nodes>\n");

	int node_num=graph.node_num;
	int edge_num=graph.edge_num;
	struct Node_PF* node=graph.node;
	struct Edge_PF* edge=graph.edge;
	unsigned* off=graph.off;

	bool *node_is_show=(bool *)malloc(node_num*sizeof(bool));
	int *node_new_id=(int *)malloc(node_num*sizeof(int));
	int show_num=1000;

	for(int i=0;i<node_num;i++){
		if(i<1000/*||node_num-i<100||node[i].type==3||(rand()%node_num)<show_num*/){
			node_is_show[i]=true;
		}
		else{
			node_is_show[i]=false;
		}
	}
	/*for(int i=0,edge_count=0;i<node_num;i++){
		if(node_is_show[i]==false) continue;
		int count=0;
		for(int j=off[i];j<off[i+1];j++){
			if(edge[j].src!=i&&node_is_show[edge[j].src]==true) count++;
		}
		if(count==0){
			node_is_show[i]=false;
		}
	}*/

	for(int i=0,k=0;i<node_num;i++){
		if(node_is_show[i]==true){
			node_new_id[i]=k;
			k++;
		}
		else{
			node_new_id[i]=2333333;
		}
	}

	for(int i=0;i<node_num;i++){
		if(node_is_show[i]==false) continue;

		double x=-400.0+(rand()%800000)/1000.0;
		double y=-300.0+(rand()%600000)/1000.0;
		int type_value;
		if(node[i].type==0||node[i].type==1) type_value=0;//PQ
		else if(node[i].type==2) type_value=1;//PV
		else type_value=2;//slave
		fprintf(fp,"      <node id=\"%d\" label=\"node%d\">\n",node_new_id[i],i);
		fprintf(fp,"        <attvalues>\n");
		fprintf(fp,"          <attvalue for=\"modularity_class\" value=\"%d\"></attvalue>\n",type_value);
		//fprintf(fp,"          <attvalue for=\"node_Vm_class\" value=\"%f\"></attvalue>\n",node[i].Vm);
		//fprintf(fp,"          <attvalue for=\"node_Va_class\" value=\"%f\"></attvalue>\n",node[i].Va);

		fprintf(fp,"          <attvalue for=\"node_info_class\" value=\"\\ninput voltage:%f\\ninput angle:%f\\ninput active power:%f\\ninput reactive power:%f\\noutput voltage:%f\\noutput angle:%f\"></attvalue>\n",
						graph_in.node[i].Vm,graph_in.node[i].Va*PI/180,graph_in.node[i].P,graph_in.node[i].Q,node[i].Vm,node[i].Va);

		fprintf(fp,"        </attvalues>\n");
		fprintf(fp,"        <viz:size value=\"%d\"></viz:size>\n",10);
		fprintf(fp,"        <viz:position x=\"%f\" y=\"%f\" z=\"0.0\"></viz:position>\n",x,y);
		fprintf(fp,"        <viz:color r=\"%d\" g=\"%d\" b=\"%d\"></viz:color>\n",
				color[type_value].r,color[type_value].g,color[type_value].b);
		fprintf(fp,"      </node>\n");
	}
	fprintf(fp,"    </nodes>\n");
	fprintf(fp,"    <edges>\n");
	for(int i=0,edge_count=0;i<node_num;i++){
		for(int j=off[i];j<off[i+1];j++){
			if(node_is_show[i]==false||node_is_show[edge[j].src]==false) continue;

			fprintf(fp,"      <edge id=\"%d\" source=\"%d\" target=\"%d\" weight=\"\\nG:%f\\nB:%f\">\n",
					edge_count,node_new_id[edge[j].src],node_new_id[i],edge[j].G,edge[j].B);
			edge_count++;
			fprintf(fp,"        <attvalues></attvalues>\n");
			//fprintf(fp,"        <attvalues>\n");
			//fprintf(fp,"          <attvalue for=\"edge_info_class\" value=\"\\nG:%f\\nB:%f\"></attvalue>\n",
			//		edge[j].G,edge[j].B);
			//fprintf(fp,"        </attvalues>\n");

			fprintf(fp,"      </edge>\n");
		}
	}
	fprintf(fp,"    </edges>\n");
	fprintf(fp,"  </graph>\n");
	fprintf(fp,"</gexf>\n");

	fclose(fp);
}

queue<int> vec_node;
void power_show_csv(struct Graph_PF graph,struct Graph_PF_IN graph_in){
	char outfile_node[128]="/home/qch/src/pf_out/pf_node.csv";
	FILE *fp_node=fopen(outfile_node,"w");
	char outfile_edge[128]="/home/qch/src/pf_out/pf_edge.csv";
	FILE *fp_edge=fopen(outfile_edge,"w");

	int node_num=graph.node_num;
	int edge_num=graph.edge_num;
	struct Node_PF* node=graph.node;
	struct Edge_PF* edge=graph.edge;
	unsigned* off=graph.off;

	bool *node_is_show=(bool *)malloc(node_num*sizeof(bool));
	int *node_new_id=(int *)malloc(node_num*sizeof(int));
	double *x=(double *)malloc(node_num*sizeof(double));
	double *y=(double *)malloc(node_num*sizeof(double));

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

	/*for(int iter=0;iter<10;iter++){
		for(int i=0;i<node_num;i++){
			if(node_is_show[i]==false) continue;
			double sum_x=0,sum_y=0;
			int count=0;
			for(int j=off[i];j<off[i+1];j++){
				if(node_is_show[edge[j].src]==false) continue;
				sum_x+=x[edge[j].src];
				sum_y+=y[edge[j].src];
				count++;
			}
			if((x[i]-sum_x/count)*(x[i]-sum_x/count)+(y[i]-sum_y/count)*(y[i]-sum_y/count)<10*10){
				continue;
			}
			x[i]=0.2*x[i]+0.8*sum_x/count;
			y[i]=0.15*y[i]+0.85*sum_y/count;
		}
	}*/

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

		//double x=-400.0+(rand()%800000)/1000.0;
		//double y=-300.0+(rand()%600000)/1000.0;
		int type_value;
		if(node[i].type==0||node[i].type==1) type_value=0;//PQ
		else if(node[i].type==2) type_value=1;//PV
		else type_value=2;//slave

		if(node_count) fprintf(fp_node,"\n");
		fprintf(fp_node,"%d,%d,%d,%f,%f,%f,%f,%f,%f,%f,%f",node_new_id[i],i,type_value,x[i],y[i],graph_in.node[i].Vm,graph_in.node[i].Va*PI/180,graph_in.node[i].P,graph_in.node[i].Q,node[i].Vm,node[i].Va);
		node_count++;
	}
	for(int i=0,edge_count=0;i<node_num;i++){
		for(int j=off[i];j<off[i+1];j++){
			if(node_is_show[i]==false||node_is_show[edge[j].src]==false) continue;
			if(edge_count) fprintf(fp_edge,"\n");
			fprintf(fp_edge,"%d,%d,%d,%d,%d,%f,%f",
					edge_count,node_new_id[edge[j].src],node_new_id[i],edge[j].src,i,
					edge[j].G,edge[j].B);
			edge_count++;
		}
	}
	fflush(fp_node);
	fflush(fp_edge);
	fclose(fp_node);
	fclose(fp_edge);
}

void power_show2(struct Graph_PF graph,double precision,int iter,double excusion_time,double *total_time,double total_time_count){
	/*char outfile[128]="/home/zyh/power_data/pf_out/Case10790.json";
	FILE *fp=fopen(outfile,"w");

	fprintf(fp,"\\nnode_num:%d\\nedge_num:%d\\nprecision:%f\\niter:%d\\nexcusion_time:%fms",
			graph.node_num,graph.edge_num,precision,iter,excusion_time);

	fclose(fp);*/

	char outfile2[128]="/home/qch/src/pf_out/pf_time.csv";
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
