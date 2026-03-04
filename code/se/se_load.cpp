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

struct Graph_SE_IN se_load(string casename){
	int node_num;
	int edge_lines;
	int edge_num=0;

	string nodefile="/home/qch/powers/first_power/power/power_show/se_data/"+casename+"nodeinfo.csv";
	string edgefile="/home/qch/powers/first_power/power/power_show/se_data/"+casename+"edgeinfo.csv";

	//char nodefile[64]="/home/zyh/power_data/se_data/sc_20171208_11100nodeinfo.csv";
	//char edgefile[64]="/home/zyh/power_data/se_data/sc_20171208_11100edgeinfo.csv";
	//char nodefile[64]="/home/zyh/power_data/se_data/sc_20171207_06300nodeinfo.csv";
	//char edgefile[64]="/home/zyh/power_data/se_data/sc_20171207_06300edgeinfo.csv";
	//char nodefile[64]="/home/zyh/power_data/se_data/sc_20171207_11100nodeinfo.csv";
	//char edgefile[64]="/home/zyh/power_data/se_data/sc_20171207_11100edgeinfo.csv";
	//char nodefile[64]="/home/zyh/power_data/se_data/sc_20171207_10000nodeinfo.csv";
	//char edgefile[64]="/home/zyh/power_data/se_data/sc_20171207_10000edgeinfo.csv";
	//char nodefile[64]="/home/zyh/power_data/se_data/sc_20171207_09150nodeinfo.csv";
	//char edgefile[64]="/home/zyh/power_data/se_data/sc_20171207_09150edgeinfo.csv";
	//char nodefile[64]="/home/zyh/power_data/se_data/sc_20171207_08000nodeinfo.csv";
	//char edgefile[64]="/home/zyh/power_data/se_data/sc_20171207_08000edgeinfo.csv";

	//char nodefile[64]="/home/zyh/power_data/se_data/pf_sc_20171207_11100nodeinfo.csv";
	//char edgefile[64]="/home/zyh/power_data/se_data/pf_sc_20171207_11100edgeinfo.csv";

	//char nodefile[64]="/home/zyh/power_data/se_data/sc_20171128_174550nodeinfo.csv";
	//char edgefile[64]="/home/zyh/power_data/se_data/sc_20171128_174550edgeinfo.csv";

	//char nodefile[128]="/home/zyh/power_data/se_data/sc_20160818_1120.QS_nodeinfo_mod_se.csv";
	//char edgefile[128]="/home/zyh/power_data/se_data/sc_20160818_1120.QS_edgeinfo_mod_se.csv";
	//char nodefile[128]="/home/zyh/power_data/se_data/fj_20160805_1500_QS_nodeinfo_se.csv";
	//char edgefile[128]="/home/zyh/power_data/se_data/fj_20160805_1500_QS_edgeinfo_se.csv";

    vector< vector<string> > strNode;
    read_file(nodefile.c_str(),strNode);
    vector< vector<string> > strEdge;
    read_file(edgefile.c_str(),strEdge);

    node_num=strNode.size()-1;
    edge_lines=strEdge.size()-1;
    struct Node_SE_In *node=(struct Node_SE_In *)malloc(node_num*sizeof(struct Node_SE_In));
	struct Edge_SE_In *edge=(struct Edge_SE_In *)malloc((edge_lines*2)*sizeof(struct Edge_SE_In));
	struct Edge_SE_In *edge_pre=(struct Edge_SE_In *)malloc((edge_lines*2)*sizeof(struct Edge_SE_In));
	unsigned *off=(unsigned *)malloc((node_num+1)*sizeof(unsigned));

    //printf("node_num=%d\tedge_lines=%d\n",node_num,edge_lines);
	for(int i=0;i<node_num;i++){
        /*int bus_id;//$0
		char bus_name[64];//$1
		int area,loss_zone,type;//$2~$4
		double voltage,angle,load_P,load_Q,generation_p,generation_Q;//$5~$10
		double base_kV,desired_volts,MAX_Q,Min_Q,G,B;//$11~$16
		int control_bus_number;//$17
		double M_Vm,M_Va,M_P,M_Q,Ri_V,Ri_vP,Ri_vQ;//$18~$24  */
		int type;//$4
		double M_Vm;//$18
		double M_Va;//$19
        double M_P;//$20
        double M_Q;//$21
        double G;//$15
        double B;//$16
        double Ri_V;//$22
		double Ri_vP;//$23
		double Ri_vQ;//$24

		string str;
		str=strNode.at(i+1).at(4);
		type=atoi(str.c_str());
		str=strNode.at(i+1).at(18);
		M_Vm=atof(str.c_str());
		str=strNode.at(i+1).at(19);
		M_Va=atof(str.c_str());
		str=strNode.at(i+1).at(20);
		M_P=atof(str.c_str());
		str=strNode.at(i+1).at(21);
		M_Q=atof(str.c_str());
		str=strNode.at(i+1).at(15);
		G=atof(str.c_str());
		str=strNode.at(i+1).at(16);
		B=atof(str.c_str());
		str=strNode.at(i+1).at(22);
		Ri_V=atof(str.c_str());
		str=strNode.at(i+1).at(23);
		Ri_vP=atof(str.c_str());
		str=strNode.at(i+1).at(24);//if(i==2) cout << str << endl;
		Ri_vQ=atof(str.c_str());

		node[i].type=type;
		node[i].Vm=M_Vm;
		node[i].Va=M_Va;
		node[i].P=M_P;
		node[i].Q=M_Q;
		node[i].G=G;
		node[i].B=B;
		node[i].Ri_V=Ri_V;
		node[i].Ri_vP=Ri_vP;
		node[i].Ri_vQ=Ri_vQ;

		node[i].M_Vm=M_Vm;
	}

    for(int i=0;i<edge_lines;i++){
        /*int ap_bus,z_bus,area,zone,circuit,type;//$0~$5
        double R,X,B,Line_Q1,Line_Q2,Line_Q3;//$6~$11
        int control_bus,side;//$12~$13
        double transformer_final_turns_ratio,transformer_final_angle,Min_tap,Max_tap;//$14~$17
        int step_size;//$18
        double Min_volt,Max_volt;//$19~$20
        double M_P_TLPF,M_Q_TLPF,M_P_TLPF_reverse,M_Q_TLPF_reverse,Ri_eP,Ri_eQ,Ri_eP_reverse,Ri_eQ_reverse;//$21~$28  */

        int ap_bus,z_bus;
        double R,X,B;
        double transformer_final_turns_ratio;
        double M_P_TLPF,M_Q_TLPF,M_P_TLPF_reverse,M_Q_TLPF_reverse,Ri_eP,Ri_eQ,Ri_eP_reverse,Ri_eQ_reverse;

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

        str=strEdge.at(i+1).at(21);//if(i==0) cout << str << atof(str.c_str()) << endl;
		M_P_TLPF=atof(str.c_str());
		str=strEdge.at(i+1).at(22);//if(i==0) cout << str << atof(str.c_str()) << endl;
		M_Q_TLPF=atof(str.c_str());
		str=strEdge.at(i+1).at(23);//if(i==0) cout << str << atof(str.c_str()) << endl;
		M_P_TLPF_reverse=atof(str.c_str());
		str=strEdge.at(i+1).at(24);//if(i==0) cout << str << atof(str.c_str()) << endl;
		M_Q_TLPF_reverse=atof(str.c_str());
		str=strEdge.at(i+1).at(25);
		Ri_eP=atof(str.c_str());
		str=strEdge.at(i+1).at(26);
		Ri_eQ=atof(str.c_str());
		str=strEdge.at(i+1).at(27);
		Ri_eP_reverse=atof(str.c_str());
		str=strEdge.at(i+1).at(28);
		Ri_eQ_reverse=atof(str.c_str());

        edge[i].from=ap_bus-1;
	    edge[i].to=z_bus-1;
	    edge[i].G=R/(R*R+X*X);
	    edge[i].B=X/(R*R+X*X);
	    edge[i].hB=B;
	    edge[i].K=transformer_final_turns_ratio;
	    edge[i].Kcount=1;
	    edge[i].BIJ=1/X;
	    edge[i].M_P_TLPF=M_P_TLPF;
	    edge[i].M_Q_TLPF=M_Q_TLPF;
	    edge[i].Ri_eP=Ri_eP;
	    edge[i].Ri_eQ=Ri_eQ;
	    edge[i].M_P_TLPF_reverse=M_P_TLPF_reverse;
	    edge[i].M_Q_TLPF_reverse=M_Q_TLPF_reverse;
	    edge[i].Ri_eP_reverse=Ri_eP_reverse;
	    edge[i].Ri_eQ_reverse=Ri_eQ_reverse;


	    edge[i+edge_lines].from=z_bus-1;
	    edge[i+edge_lines].to=ap_bus-1;
	    edge[i+edge_lines].G=R/(R*R+X*X);
	    edge[i+edge_lines].B=X/(R*R+X*X);
	    edge[i+edge_lines].hB=B;
	    edge[i+edge_lines].K=-transformer_final_turns_ratio;
	    edge[i+edge_lines].Kcount=1;
	    edge[i+edge_lines].BIJ=1/X;
	    edge[i+edge_lines].M_P_TLPF=M_P_TLPF_reverse;
	    edge[i+edge_lines].M_Q_TLPF=M_Q_TLPF_reverse;
	    edge[i+edge_lines].Ri_eP=Ri_eP_reverse;
	    edge[i+edge_lines].Ri_eQ=Ri_eQ_reverse;
	    edge[i+edge_lines].M_P_TLPF_reverse=M_P_TLPF;
	    edge[i+edge_lines].M_Q_TLPF_reverse=M_Q_TLPF;
	    edge[i+edge_lines].Ri_eP_reverse=Ri_eP;
	    edge[i+edge_lines].Ri_eQ_reverse=Ri_eQ;
    }

#define print_node_info(node,node_num) \
	printf("graph_in:node_num=%d\n",node_num);\
	printf("node:\n");\
	for(int i=0;i<node_num;i++){\
		if(i>=5&&node_num-i>5) continue;\
		printf("%d\t%d\t",i,node[i].type);\
		printf("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",\
         node[i].Vm,node[i].Va,node[i].P,node[i].Q,node[i].G,node[i].B,node[i].Ri_V,node[i].Ri_vP,node[i].Ri_vQ);\
	}

#define print_edge_info(edge_pre,edge_num) \
	printf("graph_in:edge_num:%d\n",edge_num);\
	for(int i=0;i<edge_num;i++){\
	    if(i>=5&&edge_num-i>5) continue;\
	    printf("%d\t%d\t%d\t",i,edge_pre[i].from,edge_pre[i].to);\
	    printf("%f\t%f\t%f\t%f\t",edge_pre[i].G,edge_pre[i].B,edge_pre[i].hB,edge_pre[i].K);\
	    printf("%d\t",edge_pre[i].Kcount);\
	    printf("%f\t",edge_pre[i].BIJ);\
	    printf("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",\
            edge_pre[i].M_P_TLPF,edge_pre[i].M_Q_TLPF,edge_pre[i].Ri_eP,edge_pre[i].Ri_eQ,\
            edge_pre[i].M_P_TLPF_reverse,edge_pre[i].M_Q_TLPF_reverse,edge_pre[i].Ri_eP_reverse,edge_pre[i].Ri_eQ_reverse);\
	}

	//print_node_info(node,node_num)
	//print_edge_info(edge,2*edge_lines)

	struct Edge_sort
	{
		inline bool operator() (const Edge_SE_In& tuple1, const Edge_SE_In& tuple2)
	    {
			return (tuple1.from < tuple2.from
					||(tuple1.from == tuple2.from&&tuple1.to < tuple2.to));
	    }
	};
	std::sort(edge, edge+2*edge_lines, Edge_sort());

	//print_edge_info(edge,2*edge_lines)

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

			edge_pre[edge_num-1].M_P_TLPF+=edge[i].M_P_TLPF;
			edge_pre[edge_num-1].M_Q_TLPF+=edge[i].M_Q_TLPF;
			edge_pre[edge_num-1].M_P_TLPF_reverse+=edge[i].M_P_TLPF_reverse;
			edge_pre[edge_num-1].M_Q_TLPF_reverse+=edge[i].M_Q_TLPF_reverse;
		}
	}

	for(int i=0,cur=-1;i<edge_num;i++){
	    while(edge_pre[i].from!=cur){
	        assert(edge_pre[i].from-cur==1);
	        cur++;
	        off[cur]=i;
	    }
	}
	off[node_num]=edge_num;

//debug
#define print_off(off,node_num) \
	printf("off:\n");\
	for(int i=0;i<node_num+1;i++){\
		if(i>=5&&node_num+1-i>5) continue;\
		printf("%d\t%d\n",i,off[i]);\
	}

	//print_edge_info(edge_pre,edge_num)
	//print_off(off,node_num)


	free(edge);
	struct Graph_SE_IN graph_in;
	graph_in.node=node;
	graph_in.edge=edge_pre;
	graph_in.off=off;
	graph_in.node_num=node_num;
	graph_in.edge_num=edge_num;
	return graph_in;

	/*FILE *fp_node=fopen(nodefile,"r");
	char temp[512];
	fgets(temp,512,fp_node);
	for(int i=0;i<node_num;i++){
		int bus_id;
		char bus_name[64];
		int area,loss_zone,type;
		double voltage,angle,load_P,load_Q,generation_p,generation_Q;
		double base_kV,desired_volts,MAX_Q,Min_Q,G,B;
		int control_bus_number;
		double M_Vm,M_Va,M_P,M_Q,Ri_V,Ri_vP,Ri_vQ;

		fscanf(fp_node,"%d",&bus_id);fseek(fp_node, 1L, SEEK_CUR);
		fscanf(fp_node,"%s",bus_name);fseek(fp_node, 1L, SEEK_CUR);
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

		fscanf(fp_node,"%d",&control_bus_number);fseek(fp_node, 1L, SEEK_CUR);

		fscanf(fp_node,"%lf",&M_Vm);fseek(fp_node, 1L, SEEK_CUR);
		fscanf(fp_node,"%lf",&M_Va);fseek(fp_node, 1L, SEEK_CUR);
		fscanf(fp_node,"%lf",&M_P);fseek(fp_node, 1L, SEEK_CUR);
		fscanf(fp_node,"%lf",&M_Q);fseek(fp_node, 1L, SEEK_CUR);
		fscanf(fp_node,"%lf",&Ri_V);fseek(fp_node, 1L, SEEK_CUR);
		fscanf(fp_node,"%lf",&Ri_vP);fseek(fp_node, 1L, SEEK_CUR);
		fscanf(fp_node,"%lf",&Ri_vQ);fseek(fp_node, 1L, SEEK_CUR);

		node[i].type=type;
		node[i].Vm=M_Vm;
		node[i].Va=M_Va;
		node[i].P=M_P;//P-Pl=$3-$10=$9-$7
		node[i].Q=M_Q;//Q-Ql=$4-$11=$10-$8
		node[i].G=G;
		node[i].B=B;
		node[i].Ri_vP=Ri_vP;
		node[i].Ri_vQ=Ri_vQ;
	}
	fclose(fp_node);

	FILE *fp_edge=fopen(edgefile,"r");
	fgets(temp,512,fp_edge);
	for(int i = 0 ;i < edge_lines; i++){
		 int ap_bus,z_bus,area,zone,circuit,type;
		 double R,X,B,Line_Q1,Line_Q2,Line_Q3;
		 int control_bus,side;
		 double transformer_final_turns_ratio,transformer_final_angle,Min_tap,Max_tap;
		 int step_size;
		 double Min_volt,Max_volt;
		 double M_P_TLPF,M_Q_TLPF,M_P_TLPF_reverse,M_Q_TLPF_reverse,Ri_eP,Ri_eQ,Ri_eP_reverse,Ri_eQ_reverse;



		fscanf(fp_edge,"%d",&ap_bus);fseek(fp_edge, 1L, SEEK_CUR);
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
	    fscanf(fp_edge,"%lf",&Max_volt);fseek(fp_edge, 1L, SEEK_CUR);

	    fscanf(fp_edge,"%lf",&M_P_TLPF);fseek(fp_edge, 1L, SEEK_CUR);
	    fscanf(fp_edge,"%lf",&M_Q_TLPF);fseek(fp_edge, 1L, SEEK_CUR);
	    fscanf(fp_edge,"%lf",&M_P_TLPF_reverse);fseek(fp_edge, 1L, SEEK_CUR);
	    fscanf(fp_edge,"%lf",&M_Q_TLPF_reverse);fseek(fp_edge, 1L, SEEK_CUR);
	    fscanf(fp_edge,"%lf",&Ri_eP);fseek(fp_edge, 1L, SEEK_CUR);
	    fscanf(fp_edge,"%lf",&Ri_eQ);fseek(fp_edge, 1L, SEEK_CUR);
	    fscanf(fp_edge,"%lf",&Ri_eP_reverse);fseek(fp_edge, 1L, SEEK_CUR);
	    fscanf(fp_edge,"%lf",&Ri_eQ_reverse);fseek(fp_edge, 1L, SEEK_CUR);

	    edge[i].from=ap_bus-1;
	    edge[i].to=z_bus-1;
	    edge[i].G=R/(R*R+X*X);
	    edge[i].B=X/(R*R+X*X);
	    edge[i].hB=B;
	    edge[i].K=transformer_final_turns_ratio;
	    edge[i].Kcount=1;
	    edge[i].BIJ=1/X;
	    edge[i].M_P_TLPF=M_P_TLPF;
	    edge[i].M_Q_TLPF=M_Q_TLPF;
	    edge[i].Ri_eP=Ri_eP;
	    edge[i].Ri_eQ=Ri_eQ;
	    edge[i].M_P_TLPF_reverse=M_P_TLPF_reverse;
	    edge[i].M_Q_TLPF_reverse=M_Q_TLPF_reverse;
	    edge[i].Ri_eP_reverse=Ri_eP_reverse;
	    edge[i].Ri_eQ_reverse=Ri_eQ_reverse;


	    edge[i+edge_lines].from=z_bus-1;
	    edge[i+edge_lines].to=ap_bus-1;se_load
	    edge[i+edge_lines].G=R/(R*R+X*X);
	    edge[i+edge_lines].B=X/(R*R+X*X);
	    edge[i+edge_lines].hB=B;
	    edge[i+edge_lines].K=-transformer_final_turns_ratio;
	    edge[i+edge_lines].Kcount=1;
	    edge[i+edge_lines].BIJ=1/X;
	    edge[i].M_P_TLPF=M_P_TLPF_reverse;
	    edge[i].M_Q_TLPF=M_Q_TLPF_reverse;
	    edge[i].Ri_eP=Ri_eP_reverse;
	    edge[i].Ri_eQ=Ri_eQ_reverse;
	    edge[i].M_P_TLPF_reverse=M_P_TLPF;
	    edge[i].M_Q_TLPF_reverse=M_Q_TLPFe;
	    edge[i].Ri_eP_reverse=Ri_eP;
	    edge[i].Ri_eQ_reverse=Ri_eQ;
	}
	fclose(fp_edge);*/
}

