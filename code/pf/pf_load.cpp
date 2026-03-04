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
#include <fstream>
#include <sstream>
#include <string>

using namespace std;

void read_file(const char *filename,vector< vector<string> > &strArray){
    ifstream inFile(filename, ios::in);
    string lineStr;
    //vector< vector<string> > strArray;
    while (getline(inFile, lineStr))
    {
        //cout << lineStr << endl;
        stringstream ss(lineStr);
        string str;
        vector<string> lineArray;
        while (getline(ss, str, ','))
            lineArray.push_back(str);
        if(lineArray.at(0).empty()==false)
        	strArray.push_back(lineArray);
    }
}

vector< vector<string> > strNode;
vector< vector<string> > strEdge;
struct Graph_PF_IN pf_load(string casename){
	int node_num;
	int edge_lines;
	int edge_num=0;

	string nodefile="/home/qch/powers/first_power/power/power_data/pf_data/"+casename+"_nodeinfo.csv";
	string edgefile="/home/qch/powers/first_power/power/power_data/pf_data/"+casename+"_edgeinfo.csv";

	//char nodefile[128]="/home/zyh/power_data/pf_data/Case10790_Nodeinfo_YT.csv";
	//char edgefile[128]="/home/zyh/power_data/pf_data/Case10790_Edgeinfo_YT.csv";

	//char nodefile[128]="/home/zyh/power_data/pf_data/IEEE118_nodeinfo.csv";
	//char edgefile[128]="/home/zyh/power_data/pf_data/IEEE118_edgeinfo.csv";

	//char nodefile[128]="/home/zyh/power_data/pf_data/IEEE30_nodeinfo.csv";
	//char edgefile[128]="/home/zyh/power_data/pf_data/IEEE30_edgeinfo.csv";

	//char nodefile[128]="/home/zyh/power_data/pf_data/IEEE14_node.csv";
	//char edgefile[128]="/home/zyh/power_data/pf_data/IEEE14_branchs.csv";



	//string nodefile="/home/zyh/power_data/pf_data/sc_20160818_1120.QS_nodeinfo_mod.csv";
	//string edgefile="/home/zyh/power_data/pf_data/sc_20160818_1120.QS_edgeinfo_mod.csv";

	//string nodefile="/home/zyh/power_data/se_data/pf_sc_20171207_11100nodeinfo.csv";
	//string edgefile="/home/zyh/power_data/se_data/pf_sc_20171207_11100edgeinfo.csv";


	read_file(nodefile.c_str(),strNode);
	read_file(edgefile.c_str(),strEdge);

	node_num=strNode.size()-1;
	edge_lines=strEdge.size()-1;

	struct Node_PF_In *node=(struct Node_PF_In *)malloc(node_num*sizeof(struct Node_PF_In));
	struct Edge_PF_In *edge=(struct Edge_PF_In *)malloc((edge_lines*2)*sizeof(struct Edge_PF_In));
	struct Edge_PF_In *edge_pre=(struct Edge_PF_In *)malloc((edge_lines*2)*sizeof(struct Edge_PF_In));
	unsigned *off=(unsigned *)malloc((node_num+1)*sizeof(unsigned));

	//printf("node_num=%d\tedge_lines=%d\n",node_num,edge_lines);

	//FILE *fp_node=fopen(nodefile,"r");
	//char temp[512];
	//fgets(temp,512,fp_node);
	for(int i=0;i<node_num;i++){
		int bus_id,bus_name,area,loss_zone,type;
		double voltage,angle,load_P,load_Q,generation_p,generation_Q;
		double base_kV,desired_volts,MAX_Q,Min_Q,G,B;
		int control_bus_number;

		string str;
		str=strNode.at(i+1).at(4);
		type=atoi(str.c_str());
		str=strNode.at(i+1).at(5);
		voltage=atof(str.c_str());
		str=strNode.at(i+1).at(6);
		angle=atof(str.c_str());
		str=strNode.at(i+1).at(7);
		load_P=atof(str.c_str());
		str=strNode.at(i+1).at(8);
		load_Q=atof(str.c_str());
		str=strNode.at(i+1).at(9);
		generation_p=atof(str.c_str());
		str=strNode.at(i+1).at(10);
		generation_Q=atof(str.c_str());
		str=strNode.at(i+1).at(15);
		G=atof(str.c_str());
		str=strNode.at(i+1).at(16);
		B=atof(str.c_str());

		/*fscanf(fp_node,"%d",&bus_id);fseek(fp_node, 1L, SEEK_CUR);
		fscanf(fp_node,"%d",&bus_name);fseek(fp_node, 1L, SEEK_CUR);
		fscanf(fp_node,"%d",&area);fseek(fp_node, 1L, SEEK_CUR);
		fscanf(fp_node,"%d",&loss_zone);fseek(fp_node, 1L, SEEK_CUR);
		fscanf(fp_node,"%d",&type);fseek(fp_node, 1L, SEEK_CUR);

		fscanf(fp_node,"%lf",&voltage);fseek(fp_node, 1L, SEEK_CUR);
		fscanf(fp_node,"%lf",&angle);fseek(fp_node, 1L, SEEK_CUR);
		fscanf(fp_node,"%lf",&load_P);fseek(fp_node, 1L, SEEK_CUR);
		fscanf(fp_node,"%lf",&load_Q);fseek(fp_node, 1L, SEEK_CUR);
		fscanf(fp_node,"%lf",&generation_p);fseek(fp_node, 1L, SEEK_CUR);
		fscanf(fp_node,"%lf",&generation_Q);fseek(fp_node, 1L, SEEK_CUR);

		fscanf(fp_node,"%lf",&base_kV);fseek(fp_node, 1L, SEEK_CUR);
		fscanf(fp_node,"%lf",&desired_volts);fseek(fp_node, 1L, SEEK_CUR);
		fscanf(fp_node,"%lf",&MAX_Q);fseek(fp_node, 1L, SEEK_CUR);
		fscanf(fp_node,"%lf",&Min_Q);fseek(fp_node, 1L, SEEK_CUR);
		fscanf(fp_node,"%lf",&G);fseek(fp_node, 1L, SEEK_CUR);
		fscanf(fp_node,"%lf",&B);fseek(fp_node, 1L, SEEK_CUR);

		fscanf(fp_node,"%d",&control_bus_number);fseek(fp_node, 1L, SEEK_CUR);*/

		node[i].type=type;
		node[i].Vm=voltage;
		node[i].Va=angle;
		node[i].P=generation_p-load_P;//P-Pl=$3-$10=$9-$7
		node[i].Q=generation_Q-load_Q;//Q-Ql=$4-$11=$10-$8
		node[i].G=G;
		node[i].B=B;

		//if(angle>=180||angle<-180) {printf("node %d:angle>=180||angle<-180",i);exit(-1);}
	}
	//fclose(fp_node);

	//FILE *fp_edge=fopen(edgefile,"r");
	//fgets(temp,512,fp_edge);
	for(int i = 0 ;i < edge_lines; i++){
		 int ap_bus,z_bus,area,zone,circuit,type;
		 double R,X,B,Line_Q1,Line_Q2,Line_Q3;
		 int control_bus,side;
		 double transformer_final_turns_ratio,transformer_final_angle,Min_tap,Max_tap;
		 int step_size;
		 double Min_volt,Max_volt;

		 string str;
		 str=strEdge.at(i+1).at(0);
		 ap_bus=atoi(str.c_str());
		 str=strEdge.at(i+1).at(1);
		 z_bus=atoi(str.c_str());

		 str=strEdge.at(i+1).at(6);
		 R=atof(str.c_str());
		 str=strEdge.at(i+1).at(7);
		 X=atof(str.c_str());
		 str=strEdge.at(i+1).at(8);
		 B=atof(str.c_str());

		 str=strEdge.at(i+1).at(14);
		 transformer_final_turns_ratio=atof(str.c_str());

		 /*fscanf(fp_edge,"%d",&ap_bus);fseek(fp_edge, 1L, SEEK_CUR);
	    fscanf(fp_edge,"%d",&z_bus);fseek(fp_edge, 1L, SEEK_CUR);
	    fscanf(fp_edge,"%d",&area);fseek(fp_edge, 1L, SEEK_CUR);
	    fscanf(fp_edge,"%d",&zone);fseek(fp_edge, 1L, SEEK_CUR);
	    fscanf(fp_edge,"%d",&circuit);fseek(fp_edge, 1L, SEEK_CUR);
	    fscanf(fp_edge,"%d",&type);fseek(fp_edge, 1L, SEEK_CUR);

	    fscanf(fp_edge,"%lf",&R);fseek(fp_edge, 1L, SEEK_CUR);
	    fscanf(fp_edge,"%lf",&X);fseek(fp_edge, 1L, SEEK_CUR);
	    fscanf(fp_edge,"%lf",&B);fseek(fp_edge, 1L, SEEK_CUR);
	    fscanf(fp_edge,"%lf",&Line_Q1);fseek(fp_edge, 1L, SEEK_CUR);
	    fscanf(fp_edge,"%lf",&Line_Q2);fseek(fp_edge, 1L, SEEK_CUR);
	    fscanf(fp_edge,"%lf",&Line_Q3);fseek(fp_edge, 1L, SEEK_CUR);

	    fscanf(fp_edge,"%d",&control_bus);fseek(fp_edge, 1L, SEEK_CUR);
	    fscanf(fp_edge,"%d",&side);fseek(fp_edge, 1L, SEEK_CUR);
	    fscanf(fp_edge,"%lf",&transformer_final_turns_ratio);fseek(fp_edge, 1L, SEEK_CUR);
	    fscanf(fp_edge,"%lf",&transformer_final_angle);fseek(fp_edge, 1L, SEEK_CUR);
	    fscanf(fp_edge,"%lf",&Min_tap);fseek(fp_edge, 1L, SEEK_CUR);
	    fscanf(fp_edge,"%lf",&Max_tap);fseek(fp_edge, 1L, SEEK_CUR);

	    fscanf(fp_edge,"%d",&step_size);fseek(fp_edge, 1L, SEEK_CUR);

	    fscanf(fp_edge,"%lf",&Min_volt);fseek(fp_edge, 1L, SEEK_CUR);
	    fscanf(fp_edge,"%lf",&Max_volt);fseek(fp_edge, 1L, SEEK_CUR);*/

	    edge[i].from=ap_bus-1;
	    edge[i].to=z_bus-1;
	    edge[i].G=R/(R*R+X*X);
	    edge[i].B=X/(R*R+X*X);
	    edge[i].hB=B;
	    edge[i].K=transformer_final_turns_ratio;
	    edge[i].Kcount=1;
	    edge[i].BIJ=1/X;

	    edge[i+edge_lines].from=z_bus-1;
	    edge[i+edge_lines].to=ap_bus-1;
	    edge[i+edge_lines].G=R/(R*R+X*X);
	    edge[i+edge_lines].B=X/(R*R+X*X);
	    edge[i+edge_lines].hB=B;
	    edge[i+edge_lines].K=-transformer_final_turns_ratio;
	    edge[i+edge_lines].Kcount=1;
	    edge[i+edge_lines].BIJ=1/X;
	}
	//fclose(fp_edge);

/*
	printf("graph_in:node_num=%d\tedge_num=%d\n",node_num,edge_num);
	printf("node:\n");
	for(int i=0;i<node_num;i++){
		if(i>=5&&node_num-i>5) continue;
		printf("%d\t%d\t",i,node[i].type);
		printf("%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",node[i].Vm,node[i].Va,node[i].P,node[i].Q,node[i].G,node[i].B);
	}
	printf("edge:\n");
	for(int i=0;i<2*edge_lines;i++){
	    if(i>=5&&2*edge_lines-i>5) continue;
	    printf("%d\t%d\t%d\t",i,edge[i].from,edge[i].to);
	    printf("%lf\t%lf\t%lf\t%lf\t",edge[i].G,edge[i].B,edge[i].hB,edge[i].K);
	    printf("%d\t",edge[i].Kcount);
	    printf("%lf\n",edge[i].BIJ);
	}
*/

	struct Edge_sort
	{
		inline bool operator() (const Edge_PF_In& tuple1, const Edge_PF_In& tuple2)
	    {
			return (tuple1.from < tuple2.from
					||(tuple1.from == tuple2.from&&tuple1.to < tuple2.to));
	    }
	};
	std::sort(edge, edge+2*edge_lines, Edge_sort());

	edge_num=0;
	for(int i=0,cur_from=-1,cur_to=-1;i<2*edge_lines;i++){
		if(edge[i].from!=cur_from||edge[i].to!=cur_to){
			edge_num++;
			edge_pre[edge_num-1]=edge[i];
			cur_from=edge[i].from;
			cur_to=edge[i].to;
		}
		else{//printf("~~~~~%d\n",i-edge_num);
			edge_pre[edge_num-1].G+=edge[i].G;
			edge_pre[edge_num-1].B+=edge[i].B;
			edge_pre[edge_num-1].hB+=edge[i].hB;
			edge_pre[edge_num-1].K+=edge[i].K;
			edge_pre[edge_num-1].Kcount+=edge[i].Kcount;
			edge_pre[edge_num-1].BIJ+=edge[i].BIJ;
		}
	}
	free(edge);

	for(int i=0,cur=-1;i<edge_num;i++){
	    while(edge_pre[i].from!=cur){
	        assert(edge_pre[i].from-cur==1);
	        cur++;
	        off[cur]=i;
	    }
	}
	off[node_num]=edge_num;
/*
//debug
	printf("graph_in:node_num=%d\tedge_num=%d\n",node_num,edge_num);
	printf("node:\n");
	for(int i=0;i<node_num;i++){
		if(i>=5&&node_num-i>5) continue;
		printf("%d\t%d\t",i,node[i].type);
		printf("%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",node[i].Vm,node[i].Va,node[i].P,node[i].Q,node[i].G,node[i].B);
	}
	printf("edge:\n");
	for(int i=0;i<edge_num;i++){
	    if(i>=5&&edge_num-i>5) continue;
	    printf("%d\t%d\t%d\t",i,edge_pre[i].from,edge_pre[i].to);
	    printf("%lf\t%lf\t%lf\t%lf\t",edge_pre[i].G,edge_pre[i].B,edge_pre[i].hB,edge_pre[i].K);
	    printf("%d\t",edge_pre[i].Kcount);
	    printf("%lf\n",edge_pre[i].BIJ);
	}
	printf("off:\n");
	for(int i=0;i<node_num+1;i++){
		if(i>=5&&node_num+1-i>5) continue;
		printf("%d\t%d\n",i,off[i]);
	}*/

	struct Graph_PF_IN graph_in;
	graph_in.node=node;
	graph_in.edge=edge_pre;
	graph_in.off=off;
	graph_in.node_num=node_num;
	graph_in.edge_num=edge_num;
	return graph_in;
}
