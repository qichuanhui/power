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
using namespace std;

struct Neighbor{
    int src;
    double B,BIJ;
    double Pweight,Qweight;
};
struct Matrix_H{
    int src;
    double value;
};

struct Ret_H se_H_seq(struct Graph_SE_IN graph_in){
	int node_num=graph_in.node_num;
	int edge_num=graph_in.edge_num;
	//struct Sort_Id_Vertex* vertex;
	struct Neighbor *neighbor;
	double *sumG,*sumB,*sumBedge,*sumBi;
	double *sumGii_P,*sumGii_Q;
	int slackbus;

	neighbor=(struct Neighbor*)malloc((graph_in.edge_num)*sizeof(struct Neighbor));
	sumG=(double*)malloc(graph_in.node_num*sizeof(double));
	sumB=(double*)malloc(graph_in.node_num*sizeof(double));
	sumBedge=(double*)malloc(graph_in.node_num*sizeof(double));
	sumBi=(double*)malloc(graph_in.node_num*sizeof(double));
	sumGii_P=(double*)malloc(graph_in.node_num*sizeof(double));
	sumGii_Q=(double*)malloc(graph_in.node_num*sizeof(double));

	for(int i=0;i<node_num;i++){
        sumGii_P[i]=0;sumGii_Q[i]=0;
        sumG[i]=0;
        sumB[i]=0;
        sumBedge[i]=0;
        sumBi[i]=0;
	}

	for(int v=0;v<node_num;v++){
		double neighbor_B,neighbor_BIJ;
		for(int j=graph_in.off[v];j<graph_in.off[v+1];j++){
			struct Edge_SE_In e= graph_in.edge[j];
			int t=e.to;

			if(e.K==0){
				sumG[v] += e.G;
				sumB[v] += -1*e.B + 0.5*e.hB;
				sumBedge[v] += e.B;
				sumBi[v] += -1*e.BIJ;

				sumGii_Q[v] += (e.B - e.hB)*(e.B - e.hB)*e.Ri_eQ;
				sumGii_Q[t] += (e.B * e.B)*graph_in.node[v].Ri_vQ + (e.B*e.B)*e.Ri_eQ;

				neighbor_B = e.B;
			}
			else if(e.K>0){
				double tap_ratio = e.K/e.Kcount;
				double tap_ratio_square = (e.K/e.Kcount)*(e.K/e.Kcount);
				sumG[v] += 1/(tap_ratio_square)*e.G;
				sumB[v] += 1/(tap_ratio_square)*(-1*e.B + 0.5*e.hB); // sqrt
				sumBedge[v] += e.B/tap_ratio;
				sumBi[v] += -1*e.BIJ;

				sumGii_Q[v] += (e.B/tap_ratio - e.hB)*(e.B/tap_ratio - e.hB)*e.Ri_eQ;
                sumGii_Q[t] += ( e.B * e.B / tap_ratio_square)*graph_in.node[v].Ri_vQ + ( e.B * e.B / tap_ratio_square)*e.Ri_eQ;

                neighbor_B = e.B/tap_ratio;
			}
			else{
				double tap_ratio = fabs(e.K/e.Kcount);
				double tap_ratio_square = (e.K/e.Kcount)*(e.K/e.Kcount);
				sumG[v] += e.G;
				sumB[v] += -1*e.B + 0.5*e.hB;
				sumBedge[v] += e.B/tap_ratio;
				sumBi[v] += -1*e.BIJ;

				sumGii_Q[v] += (e.B/tap_ratio - e.hB)*(e.B/tap_ratio - e.hB)*e.Ri_eQ;
                sumGii_Q[t] += (e.B * e.B / tap_ratio_square)*graph_in.node[v].Ri_vQ + ( e.B * e.B / tap_ratio_square)*e.Ri_eQ;

                neighbor_B = e.B / tap_ratio;
			}
            sumGii_P[v] += (e.BIJ * e.BIJ)*e.Ri_eP; // diagonal element in row i of system Gain matrix_P: Hii'*Hii + sum(Hji'*Hji)
            sumGii_P[t] += (e.BIJ * e.BIJ)*e.Ri_eP + (e.BIJ * e.BIJ)*graph_in.node[v].Ri_vP;

            neighbor_BIJ = e.BIJ;
            neighbor[j].src=t;
            neighbor[j].B=neighbor_B;neighbor[j].BIJ=neighbor_BIJ;
            neighbor[j].Pweight=graph_in.node[t].Ri_vP;
            neighbor[j].Qweight=graph_in.node[t].Ri_vQ;

		}
		sumGii_P[v] += (sumBi[v] * sumBi[v])*graph_in.node[v].Ri_vP;
		sumG[v] += graph_in.node[v].G;
		sumB[v] += graph_in.node[v].B;
		sumGii_Q[v] += (- 2 * sumB[v] - sumBedge[v]) * (- 2 * sumB[v] - sumBedge[v]) * graph_in.node[v].Ri_vQ + 1 * graph_in.node[v].Ri_V;

		if(graph_in.node[v].type==3) slackbus=v;
		/*
		if s.flag != 3 then
			s.@Gip += sort_id (s.exId, s.@sumGii_P)
		end,
		s.@Giq += sort_id (s.exId, s.@sumGii_Q),
		//@@B[s.exId-1][s.exId-1] += s.@sumB,  // FOR TESTING...
		//@@G[s.exId-1][s.exId-1] += s.@sumG,  // FOR TESTING...
		CASE WHEN (s.flag == 3) THEN
			@@slackbus += s.exId
		END;
		*/
	}


	vector< vector<struct Matrix_H> > G_offdiag_B,G_offdiag_BIJ;
	for(int v=0;v<node_num;v++){
        vector<struct Matrix_H> tmp_B,tmp_BIJ;
        if(graph_in.node[v].type!=3){
            tmp_BIJ.push_back({v,sumGii_P[v]});
		}
		else{
            tmp_BIJ.push_back({v,1.0});//black line
		}
		tmp_B.push_back({v,sumGii_Q[v]});
        G_offdiag_B.push_back(tmp_B);
        G_offdiag_BIJ.push_back(tmp_BIJ);
	}
    for(int v=0;v<node_num;v++){
        for(int j=graph_in.off[v];j<graph_in.off[v+1];j++){
            Edge_SE_In e=graph_in.edge[j];
            int t=e.to;
            if(v!=slackbus&&t!=slackbus){
                /*for(int jj=graph_in.off[v];jj<graph_in.off[v+1];jj++){
                    for(int jjj=graph_in.off[t];jjj<graph_in.off[t+1];jjj++){
                        if(neighbor[jj].src==neighbor[jjj].src){
                            G_offdiag_BIJ.at(v).push_back({t,(neighbor[jj].BIJ*neighbor[jjj].BIJ*neighbor[jj].Pweight)});
                        }
                    }
                }
                G_offdiag_BIJ.at(v).push_back({t,(sumBi[v] * e.BIJ * graph_in.node[v].Ri_vP
                                                  + sumBi[t] * e.BIJ * graph_in.node[t].Ri_vP - e.BIJ * e.BIJ * e.Ri_eP)});
                G_offdiag_BIJ.at(t).push_back({v,(- e.BIJ * e.BIJ * e.Ri_eP)});*/

            	/*G_offdiag_BIJ.at(v).push_back({t,(- e.BIJ * e.BIJ)});
            	G_offdiag_BIJ.at(t).push_back({v,(- e.BIJ * e.BIJ)});*/
            	G_offdiag_BIJ.at(v).push_back({t,2*(- e.BIJ * e.BIJ)});
            	//if(v<10) printf("G_offdiag_BIJ:(%d,%d,%f,%d,%f)\n",v,t,2*(- e.BIJ * e.BIJ),\
            			G_offdiag_BIJ.at(v).at(G_offdiag_BIJ.at(v).size()-1).src,\
						G_offdiag_BIJ.at(v).at(G_offdiag_BIJ.at(v).size()-1).value);
            }

            /*for(int jj=graph_in.off[v];jj<graph_in.off[v+1];jj++){
                for(int jjj=graph_in.off[t];jjj<graph_in.off[t+1];jjj++){
                    if(neighbor[jj].src==neighbor[jjj].src){
                        G_offdiag_B.at(v).push_back({t,(neighbor[jj].B*neighbor[jjj].B*neighbor[jj].Qweight)});
                    }
                }
            }*/
            if(e.K==0){
                /*G_offdiag_B.at(v).push_back({t,((2 * sumB[v] + sumBedge[v]) * e.B * graph_in.node[v].Ri_vQ
                                            + (2 * sumB[t] + sumBedge[t]) * e.B * graph_in.node[t].Ri_vQ - (e.B - e.hB) * e.B * e.Ri_eQ)});
                G_offdiag_B.at(t).push_back({v,(- (e.B - e.hB) * e.B * e.Ri_eQ)});*/

            	/*G_offdiag_B.at(v).push_back({t,(- (e.B - e.hB) * e.B)});
            	G_offdiag_B.at(t).push_back({v,(- (e.B - e.hB) * e.B)});*/
            	G_offdiag_B.at(v).push_back({t,2*(- (e.B - e.hB) * e.B)});
            	//if(v<10) printf("G_offdiag_B:(%d,%u,%f,%d,%f)\n",v,t,2*(- (e.B - e.hB) * e.B),\
            			G_offdiag_B.at(v).at(G_offdiag_B.at(v).size()-1).src,\
						G_offdiag_B.at(v).at(G_offdiag_B.at(v).size()-1).value);
            }
            else{
                double tap_ratio = fabs(e.K/e.Kcount);
                /*G_offdiag_B.at(v).push_back({t,((2 * sumB[v] + sumBedge[v]) * e.B / tap_ratio * graph_in.node[v].Ri_vQ
                                            + (2 * sumB[t] + sumBedge[t]) * e.B / tap_ratio * graph_in.node[t].Ri_vQ - (e.B / tap_ratio - e.hB) * e.B / tap_ratio * e.Ri_eQ)});
                G_offdiag_B.at(t).push_back({v,(- (e.B / tap_ratio - e.hB) * e.B / tap_ratio * e.Ri_eQ)});*/

                /*G_offdiag_B.at(v).push_back({t,(- (e.B / tap_ratio - e.hB) * e.B / tap_ratio)});
                G_offdiag_B.at(t).push_back({v,(- (e.B / tap_ratio - e.hB) * e.B / tap_ratio)});*/
                G_offdiag_B.at(v).push_back({t,2*(- (e.B / tap_ratio - e.hB) * e.B / tap_ratio)});
                //if(v<10) printf("G_offdiag_B:(%d,%u,%f,%d,%f)\n",v,t,2*(- (e.B / tap_ratio - e.hB) * e.B / tap_ratio),\
                		G_offdiag_B.at(v).at(G_offdiag_B.at(v).size()-1).src,\
                		G_offdiag_B.at(v).at(G_offdiag_B.at(v).size()-1).value);
            }

            // Find 2-step-neighbor info and save
            /*if(v!=slackbus){
                for(int jjj=graph_in.off[t];jjj<graph_in.off[t+1];jjj++){
                    bool flag=false;
                    if(neighbor[jjj].src!=v||neighbor[jjj].src!=slackbus){
                        flag=true;
                        for(int jj=graph_in.off[v];jj<graph_in.off[v+1];jj++){
                            if(neighbor[jjj].src==neighbor[jj].src){
                                flag=false;
                            }
                        }
                    }
                    if(flag){
                        G_offdiag_BIJ.at(v).push_back({neighbor[jjj].src,(e.BIJ * neighbor[jjj].BIJ * graph_in.node[t].Ri_vP)});
                    }
                }
            }*/

            /*for(int jjj=graph_in.off[t];jjj<graph_in.off[t+1];jjj++){
                bool flag=false;
                if(neighbor[jjj].src!=v){
                    flag=true;
                    for(int jj=graph_in.off[v];jj<graph_in.off[v+1];jj++){
                        if(neighbor[jjj].src==neighbor[jj].src){
                            flag=false;
                        }
                    }
                }
                if(flag){
                    if(e.K==0){
                        G_offdiag_B.at(v).push_back({neighbor[jjj].src,(e.B * neighbor[jjj].B * graph_in.node[t].Ri_vQ)});
                    }
                    else{
                        double tap_ratio = abs(e.K/e.Kcount);
                        G_offdiag_B.at(v).push_back({neighbor[jjj].src,(e.B / tap_ratio * neighbor[jjj].B * graph_in.node[t].Ri_vQ)});
                    }
                }
            }*/
        }
    }

    struct Gain_sort
	{
            inline bool operator() (const struct Matrix_H& tuple1, const struct Matrix_H& tuple2)
            {
                return (tuple1.src < tuple2.src);
            }
	};

    //vector< vector<struct Matrix_H> > G_offdiag_B2,G_offdiag_BIJ2;
    //for(int v=0;v<node_num;v++){
        //sort( G_offdiag_BIJ.at(v).begin(), G_offdiag_BIJ.at(v).end(),Gain_sort());
        //sort( G_offdiag_B.at(v).begin(), G_offdiag_B.at(v).end(),Gain_sort());

        /*vector<struct Matrix_H> tmp_BIJ,tmp_B;
        for(int i=0;i<G_offdiag_BIJ.at(v).size();i++){
        	if(i==0){
        		tmp_BIJ.push_back(G_offdiag_BIJ.at(v).at(i));
        	}
        	else if(tmp_BIJ.at(tmp_BIJ.size()-1).src!=G_offdiag_BIJ.at(v).at(i).src){
        		tmp_BIJ.push_back(G_offdiag_BIJ.at(v).at(i));
        	}
        	else{
        		tmp_BIJ.at(tmp_BIJ.size()-1).value+=G_offdiag_BIJ.at(v).at(i).value;
        	}
        }
        for(int i=0;i<G_offdiag_B.at(v).size();i++){
        	if(i==0){
                tmp_B.push_back(G_offdiag_B.at(v).at(i));
            }
            else if(tmp_B.at(tmp_B.size()-1).src!=G_offdiag_B.at(v).at(i).src){
                tmp_B.push_back(G_offdiag_B.at(v).at(i));
            }
            else{
            	tmp_B.at(tmp_B.size()-1).value+=G_offdiag_B.at(v).at(i).value;
            }
        }
        G_offdiag_BIJ2.push_back(tmp_BIJ);
        G_offdiag_B2.push_back(tmp_B);*/
        //G_offdiag_BIJ2.push_back(G_offdiag_BIJ.at(v));
        //G_offdiag_B2.push_back(G_offdiag_B.at(v));

        /*if(v<10){
        	for(int j=0;j<G_offdiag_BIJ.at(v).size();j++)
        		printf("(%d,%d,%f) ",v,G_offdiag_BIJ.at(v).at(j).src,G_offdiag_BIJ.at(v).at(j).value);
        	printf("\n~~~~~~~\n");

            for(int j=0;j<G_offdiag_BIJ2.at(v).size();j++)
                printf("(%d,%d,%f) ",v,G_offdiag_BIJ2.at(v).at(j).src,G_offdiag_BIJ2.at(v).at(j).value);
            printf("\n~~~~~~~\n");
        }*/
    //}

	double *Bp_x,*Bpp_x;
	unsigned *Bp_i,*Bpp_i;
	unsigned *Bp_p,*Bpp_p;
	struct Graph_SE graph;

	Bp_p=(unsigned*)malloc((graph_in.node_num+1)*sizeof(unsigned));
	Bpp_p=(unsigned*)malloc((graph_in.node_num+1)*sizeof(unsigned));
	Bp_p[0]=0;
	Bpp_p[0]=0;
	for(int v=0;v<node_num;v++){
        Bp_p[v+1]=Bp_p[v]+G_offdiag_BIJ.at(v).size();
        Bpp_p[v+1]=Bpp_p[v]+G_offdiag_B.at(v).size();
	}

	Bp_x=(double*)malloc((Bp_p[graph_in.node_num])*sizeof(double));
	Bpp_x=(double*)malloc((Bpp_p[graph_in.node_num])*sizeof(double));
	Bp_i=(unsigned*)malloc((Bp_p[graph_in.node_num])*sizeof(unsigned));
	Bpp_i=(unsigned*)malloc((Bpp_p[graph_in.node_num])*sizeof(unsigned));

	int p_k=0;
	int pp_k=0;
    for(int v=0;v<node_num;v++){
        for(int j=0;j<G_offdiag_BIJ.at(v).size();j++){
            struct Matrix_H tmp=G_offdiag_BIJ.at(v).at(j);
            Bp_i[p_k]=tmp.src;
            Bp_x[p_k]=tmp.value;//if(v<10) printf("(%d,%d,%f,%d,%f)\n",v,tmp.src,tmp.value,Bp_i[p_k],Bp_x[p_k]);
            p_k++;
        }
        for(int j=0;j<G_offdiag_B.at(v).size();j++){
            struct Matrix_H tmp=G_offdiag_B.at(v).at(j);
            Bpp_i[pp_k]=tmp.src;
            Bpp_x[pp_k]=tmp.value;//if(v<10) printf("(%d,%d,%f,%d,%f)\n",v,tmp.src,tmp.value,Bpp_i[pp_k],Bpp_x[pp_k]);
            pp_k++;
        }
    }

    /*Bp_p[0]=0;
    for(int v=0;v<node_num-1;v++){
    	if(v<slackbus)
    		Bp_p[v+1]=Bp_p[v]+G_offdiag_BIJ2.at(v).size();
    	else
    		Bp_p[v+1]=Bp_p[v]+G_offdiag_BIJ2.at(v+1).size();
    }
    p_k=0;
    for(int v=0;v<node_num;v++){
    	if(v==slackbus) continue;
        for(int j=0;j<G_offdiag_BIJ2.at(v).size();j++){
            struct Matrix_H tmp=G_offdiag_BIJ2.at(v).at(j);
            Bp_i[p_k]=tmp.src;if(tmp.src>slackbus) Bp_i[p_k]-=1;
            Bp_x[p_k]=tmp.value;
            p_k++;//if(v<10) printf("(%d,%u,%f)\t",v,Bp_i[p_k],Bp_x[p_k]);if(v==10&&j==0) printf("\n");
        }
    }*/

    //printf("p_k=%d\tpp_k=%d\t%d\t%d\n",p_k,pp_k,Bp_p[node_num],Bpp_p[node_num]);

	graph.node=(struct Node_SE *)malloc(graph_in.node_num*sizeof(struct Node_SE));
	graph.edge=(struct Edge_SE *)malloc((graph_in.edge_num)*sizeof(struct Edge_SE));
	graph.off=graph_in.off;
	graph.node_num=graph_in.node_num;
	graph.edge_num=graph_in.edge_num;

	for(int i=0;i<graph_in.node_num;i++){
        graph.node[i].type=graph_in.node[i].type;
        graph.node[i].Vm=graph_in.node[i].Vm;
        graph.node[i].Va=graph_in.node[i].Va*PI/180;
        graph.node[i].P=graph_in.node[i].P;
        graph.node[i].Q=graph_in.node[i].Q;
        graph.node[i].Ri_V=graph_in.node[i].Ri_V;
        graph.node[i].Ri_vP=graph_in.node[i].Ri_vP;
        graph.node[i].Ri_vQ=graph_in.node[i].Ri_vQ;
        graph.node[i].sumG=sumG[i];
        graph.node[i].sumB=sumB[i];
        graph.node[i].sumBedge=sumBedge[i];
        graph.node[i].sumBi=sumBi[i];

        graph.node[i].M_Vm=graph_in.node[i].M_Vm;

        graph.node[i].Va=0;//must initialize to 0


        //flatstart
        if(graph.node[i].type==0||graph.node[i].type==1||graph.node[i].type==2){
        	graph.node[i].P=graph_in.node[i].P;
            graph.node[i].Q=graph_in.node[i].Q;
            graph.node[i].Vm=1;
            graph.node[i].Va=0;
        }
        else if(graph.node[i].type==3){
            graph.node[i].P=graph_in.node[i].P;
            graph.node[i].Q=graph_in.node[i].Q;
            graph.node[i].Vm=graph_in.node[i].M_Vm;
            graph.node[i].Va=0;
        }
        else{
            graph.node[i].P=0;
            graph.node[i].Q=0;
            graph.node[i].Vm=0;
            graph.node[i].Va=0;
        }

	}
    for(int i=0;i<graph_in.edge_num;i++){
        graph.edge[i].src=graph_in.edge[i].to;
        graph.edge[i].G=graph_in.edge[i].G;
        graph.edge[i].B=graph_in.edge[i].B;
        graph.edge[i].hB=graph_in.edge[i].hB;
        graph.edge[i].K=graph_in.edge[i].K;
        graph.edge[i].Kcount=graph_in.edge[i].Kcount;
        graph.edge[i].BIJ=graph_in.edge[i].BIJ;
        graph.edge[i].M_P_TLPF=graph_in.edge[i].M_P_TLPF;
        graph.edge[i].M_Q_TLPF=graph_in.edge[i].M_Q_TLPF;
        graph.edge[i].Ri_eP=graph_in.edge[i].Ri_eP;
        graph.edge[i].Ri_eQ=graph_in.edge[i].Ri_eQ;
        graph.edge[i].M_P_TLPF_reverse=graph_in.edge[i].M_P_TLPF_reverse;
        graph.edge[i].M_Q_TLPF_reverse=graph_in.edge[i].M_Q_TLPF_reverse;
        graph.edge[i].Ri_eP_reverse=graph_in.edge[i].Ri_eP_reverse;
        graph.edge[i].Ri_eQ_reverse=graph_in.edge[i].Ri_eQ_reverse;
    }


	/*printf("graph:node_num=%d\n",graph.node_num);\
	printf("node:\n");\
	for(int i=0;i<graph.node_num;i++){\
		if(i>=5&&graph.node_num-i>5) continue;\
		printf("%d\t%d\t",i,graph.node[i].type);\
		printf("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",\
				graph.node[i].Vm,graph.node[i].Va,graph.node[i].P,graph.node[i].Q,\
				graph.node[i].Ri_V,graph.node[i].Ri_vP,graph.node[i].Ri_vQ,\
				sumG[i],sumB[i],sumBedge[i],sumBi[i],graph.node[i].M_Vm);\
	}


	("graph_in:edge_num:%d\n",edge_num);\
	for(int i=0;i<edge_num;i++){\
	    if(i>=10&&edge_num-i>10) continue;\
	    printf("%d\t%d\t%d\t",i,graph_in.edge[i].from,graph.edge[i].src);\
	    printf("%f\t%f\t%f\t%f\t",graph.edge[i].G,graph.edge[i].B,graph.edge[i].hB,graph.edge[i].K);\
	    printf("%d\t",graph.edge[i].Kcount);\
	    printf("%f\t",graph.edge[i].BIJ);\
	    printf("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",\
	    		graph.edge[i].M_P_TLPF,graph.edge[i].M_Q_TLPF,graph.edge[i].Ri_eP,graph.edge[i].Ri_eQ,\
				graph.edge[i].M_P_TLPF_reverse,graph.edge[i].M_Q_TLPF_reverse,graph.edge[i].Ri_eP_reverse,graph.edge[i].Ri_eQ_reverse);\
	}*/


	free(neighbor);
	free(sumG);
	free(sumB);
	free(sumBedge);
	free(sumBi);
	free(sumGii_P);
	free(sumGii_Q);
	//free(graph_in.node);
	//free(graph_in.edge);

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

	struct Ret_H ret_H;
	ret_H.graph=graph;
	ret_H.Bp_x=Bp_x;
	ret_H.Bp_i=Bp_i;
	ret_H.Bp_p=Bp_p;
	ret_H.Bpp_x=Bpp_x;
	ret_H.Bpp_i=Bpp_i;
	ret_H.Bpp_p=Bpp_p;
	ret_H.slackbus=slackbus;
	//printf("slackbus=%d\n",slackbus);
	//printf("make_H finished\n");
	return ret_H;
}
