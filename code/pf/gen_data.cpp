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

void pf_iter(){
	struct timeval start,finish;
	double cost_time;
	double Ybus_time=0,fac_time=0,rhs_time=0,back_for_time=0;

	struct Graph_PF_IN graph_in=pf_load("Case10790");
	gettimeofday(&start, 0);//计时开始
	struct Ret_Ybus ret_Ybus=pf_Ybus_seq(graph_in);
	gettimeofday(&finish, 0);//计时结束
	cost_time = 1000000 * (finish.tv_sec - start.tv_sec) + finish.tv_usec - start.tv_usec;
	printf("pf_Ybus:cost the time is : %lf s.\n",cost_time/(double)1000000);
	Ybus_time=cost_time/(double)1000000;

	gettimeofday(&start, 0);//计时开始
	ICktSo nicslu_p= nullptr,nicslu_pp =nullptr;
	pf_fac_seq(ret_Ybus.Bp_x,ret_Ybus.Bp_i,ret_Ybus.Bp_p,
		ret_Ybus.Bpp_x,ret_Ybus.Bpp_i,ret_Ybus.Bpp_p,
		ret_Ybus.graph.node_num,&nicslu_p,&nicslu_pp);
	gettimeofday(&finish, 0);//计时结束
	cost_time = 1000000 * (finish.tv_sec - start.tv_sec) + finish.tv_usec - start.tv_usec;
	printf("pf_fac:cost the time is : %lf s.\n",cost_time/(double)1000000);
	fac_time=cost_time/(double)1000000;

	struct Node_PF *node=ret_Ybus.graph.node;
	struct Edge_PF *edge=ret_Ybus.graph.edge;
	unsigned  *off=ret_Ybus.graph.off;
	int node_num=ret_Ybus.graph.node_num;
	int edge_num=ret_Ybus.graph.edge_num;

	double *deltaP=(double *)malloc(node_num*sizeof(double));
	double *deltaQ=(double *)malloc(node_num*sizeof(double));
	double *cosine_array=(double *)malloc(edge_num*sizeof(double));
	double *sine_array=(double *)malloc(edge_num*sizeof(double));
	int iter;
	for (iter=0; iter<10; ++iter){
	  	double maxDeltaP=0;
	  	double maxDeltaQ=0;

	  	gettimeofday(&start, 0);//计时开始
	  	// Calculate deltaP and deltaQ to update Va
	  	for (int i=0; i<node_num; ++i){
	  		deltaP[i] = 0;
	  		deltaQ[i] = 0;
	  		if (node[i].type!=3){
	  			deltaP[i] = node[i].P;
	  			if (node[i].type<2)
	  				deltaQ[i] = node[i].Q;
	  		}

	  		// Calculate network injections
	  		if (node[i].type!=3){
	  			for (int j=off[i]; j<off[i+1]; ++j){
	  				int p=edge[j].src;
	  				if (iter == 0){
	  					//if this is the first iteraton, calculate the cosine and sine terms
	  					cosine_array[j] = cos(node[i].Va-node[p].Va);
	  					sine_array[j] = sin(node[i].Va-node[p].Va);
	  					//cordic(node[i].Va-node[p].Va,&cosine_array[j],&sine_array[j]);

	  				}
	  				deltaP[i] += - node[i].Vm*node[p].Vm*(edge[j].G*cosine_array[j] + edge[j].B*sine_array[j]);
	  				if (node[i].type<2) // calculate Q for PQ buses
	  					deltaQ[i] += - node[i].Vm*node[p].Vm*(edge[j].G*sine_array[j] - edge[j].B*cosine_array[j]);
	  			}
	  		}

	  		// Get max P and max Q
	  		if (fabs(deltaP[i]) > maxDeltaP)
	  			maxDeltaP = fabs(deltaP[i]);
	  		if (fabs(deltaQ[i]) > maxDeltaQ)
	  			maxDeltaQ = fabs(deltaQ[i]);

	  		deltaP[i] /= node[i].Vm;
	  		deltaQ[i] /= node[i].Vm;
	  	}
	  	gettimeofday(&finish, 0);//计时结束
	  	cost_time = 1000000 * (finish.tv_sec - start.tv_sec) + finish.tv_usec - start.tv_usec;
	  	printf("deltaP:cost the time is : %lf s.\n",cost_time/(double)1000000);
	  	rhs_time+=cost_time/(double)1000000;

	  	std::cout << "Iter" << iter << " Updated for Va, maxDeltaP, " << maxDeltaP << ", maxDeltaQ, " << maxDeltaQ << std::endl;

	      // Decide if the program converges
	  		if ( maxDeltaP < 0.001 && maxDeltaQ < 0.001){
	  			break;
	  		}

	      //reset the max mismatche
	      maxDeltaP = 0;
	      maxDeltaQ = 0;

	      gettimeofday(&start, 0);//计时开始
	      // solve for deltaVa
	      // std::cout << " ======================== Solve for V angle ========================"<<std::endl;
	      CKTSO_Solve(nicslu_p, deltaP,deltaP, false, false);
	      // Update V angle (Va)
	      for (int i=0; i<node_num; ++i){
	  			node[i].Va -= deltaP[i];
	        //std::cout << "Iter" << iter << ", Va" << i << ", " << Va[i]<< ", deltaVa" << i << ", " << deltaP[i]<< std::endl;
	      }
	      gettimeofday(&finish, 0);//计时结束
	      cost_time = 1000000 * (finish.tv_sec - start.tv_sec) + finish.tv_usec - start.tv_usec;
	      printf("back_forward and update Va:cost the time is : %lf s.\n",cost_time/(double)1000000);
	      back_for_time+=cost_time/(double)1000000;
	      // std::cout << " ======================== Update deltaP and deltaQ ========================"<<std::endl;

	      gettimeofday(&start, 0);//计时开始
	      // Calculate deltaP and deltaQ to update Vm
	      for (int i=0; i<node_num; ++i){
	      	  deltaP[i] = 0;
	      	  deltaQ[i] = 0;

	      	  if (node[i].type!=3){
	      		  deltaP[i] = node[i].P;
	      		  if (node[i].type<2)
	      			  deltaQ[i] = node[i].Q;
	      	  }

	      	  // Calculate network injections
	      	  if (node[i].type!=3){
	      		  for (int j=off[i]; j<off[i+1]; ++j){
	      			  int p=edge[j].src;
	      			  cosine_array[j] = cos(node[i].Va-node[p].Va);
	      			  sine_array[j] = sin(node[i].Va-node[p].Va);
	      			  //cordic(node[i].Va-node[p].Va,&cosine_array[j],&sine_array[j]);

	      			  deltaP[i] += - node[i].Vm*node[p].Vm*(edge[j].G*cosine_array[j] + edge[j].B*sine_array[j]);
	      			  if (node[i].type<2) // calculate Q for PQ buses
	      				  deltaQ[i] += - node[i].Vm*node[p].Vm*(edge[j].G*sine_array[j] - edge[j].B*cosine_array[j]);
	      		  }
	      	  }

	      	  // Get max P and max Q
	      	  if (fabs(deltaP[i]) > maxDeltaP)
	      		  maxDeltaP = fabs(deltaP[i]);
	      	  	  if (fabs(deltaQ[i]) > maxDeltaQ)
	      	  		  maxDeltaQ = fabs(deltaQ[i]);

	      	  deltaP[i] /= node[i].Vm;
	      	  deltaQ[i] /= node[i].Vm;
	      }
	      gettimeofday(&finish, 0);//计时结束
	      cost_time = 1000000 * (finish.tv_sec - start.tv_sec) + finish.tv_usec - start.tv_usec;
	      printf("deltaQ:cost the time is : %lf s.\n",cost_time/(double)1000000);
	      rhs_time+=cost_time/(double)1000000;
	      std::cout << "Iter" << iter << " Updated for Vm, maxDeltaP, " << maxDeltaP << ", maxDeltaQ, " << maxDeltaQ << std::endl;

	      // Decide if the program converges
	      if ( maxDeltaP < 0.001 && maxDeltaQ < 0.001){
	  			break;
	      }

	      gettimeofday(&start, 0);//计时开始
	      // solve for deltaVm
	      // std::cout << " ======================== Solve for V magnitude ========================"<<std::endl;
	      CKTSO_Solve(nicslu_pp, deltaQ,deltaQ, false, false);

	      // Update V magnitude (Vm)
	      for (int i=0; i<node_num; ++i)
	      {
	    	  node[i].Vm -= deltaQ[i];
	    	  //std::cout << "Iter" << iter << ", Vm" << i << ", " << Vm[i]<< ", deltaVm" << i << ", " << deltaQ[i]<< std::endl;
	      }
	      gettimeofday(&finish, 0);//计时结束
	      cost_time = 1000000 * (finish.tv_sec - start.tv_sec) + finish.tv_usec - start.tv_usec;
	      printf("back_forward and update Vm:cost the time is : %lf s.\n",cost_time/(double)1000000);
	      back_for_time+=cost_time/(double)1000000;
	}
	printf("Ybus_time:%lf\tfac_time:%lf\trhs_time:%lf\tback_for_time:%lf\n",Ybus_time,fac_time,rhs_time,back_for_time);
	printf("rhs_avg:%lf\tback_for_avg:%lf\n",rhs_time/(iter+1),back_for_time/(iter+1));
}

#define RHS_THREAD_NUM 8
struct Rhs_Para_Mt{
	struct Node_PF *node;
	struct Edge_PF *edge;
	unsigned  *off;
	int node_num;
	double *deltaP;
	double *deltaQ;
	double *maxDeltaP;
	double *maxDeltaQ;
	double *cosine_array;
	double *sine_array;
	bool tri_flag;
	int local_id;
	int local_size;
};
void *thr_rhs_deltaP(void *arg){
	struct Rhs_Para_Mt *para=(struct Rhs_Para_Mt *)arg;
	struct Node_PF *node=para->node;
	struct Edge_PF *edge=para->edge;
	unsigned  *off=para->off;
	double *deltaP=para->deltaP;
	double *deltaQ=para->deltaQ;
	double *cosine_array=para->cosine_array;
	double *sine_array=para->sine_array;
	int local_id=para->local_id;
	int local_size=para->local_size;
	int node_num=para->node_num;
	int begin_node=node_num*local_id/local_size;
	int end_node=node_num*(local_id+1)/local_size;
	bool tri_flag=para->tri_flag;
	double maxDeltaP=0;
	for (int i=begin_node; i<end_node; ++i){
		deltaP[i] = 0;
		 if (node[i].type!=3){
			 deltaP[i] = node[i].P/node[i].Vm;
		 }

		 // Calculate network injections
		 if (node[i].type!=3){
		  	for (int j=off[i]; j<off[i+1]; ++j){
		  		int p=edge[j].src;
		  		if (tri_flag){
		  			//if this is the first iteraton, calculate the cosine and sine terms
		  			//cosine_array[j] = cos(node[i].Va-node[p].Va);
		  			//sine_array[j] = sin(node[i].Va-node[p].Va);

		  			cordic(node[i].Va-node[p].Va,&cosine_array[j],&sine_array[j]);
		  		}
		  		deltaP[i] += - node[p].Vm*(edge[j].G*cosine_array[j] + edge[j].B*sine_array[j]);
		  	}
		 }
	}
    return nullptr;
}

void *thr_rhs_deltaQ(void *arg){
	struct Rhs_Para_Mt *para=(struct Rhs_Para_Mt *)arg;
	struct Node_PF *node=para->node;
	struct Edge_PF *edge=para->edge;
	unsigned  *off=para->off;
	double *deltaP=para->deltaP;
	double *deltaQ=para->deltaQ;
	double *cosine_array=para->cosine_array;
	double *sine_array=para->sine_array;
	int local_id=para->local_id;
	int local_size=para->local_size;
	int node_num=para->node_num;
	int begin_node=node_num*local_id/local_size;
	int end_node=node_num*(local_id+1)/local_size;
	bool tri_flag=para->tri_flag;

	double maxDeltaP=0;
	double maxDeltaQ=0;
	for (int i=begin_node; i<end_node; ++i){
		deltaP[i] = 0;
		deltaQ[i] = 0;
		 if (node[i].type!=3){
			 deltaP[i] = node[i].P/node[i].Vm;
			 if (node[i].type<2)
				 deltaQ[i] = node[i].Q/node[i].Vm;
		 }

		 // Calculate network injections
		 if (node[i].type!=3){
		  	for (int j=off[i]; j<off[i+1]; ++j){
		  		int p=edge[j].src;
		  		if (tri_flag){
		  			//if this is the first iteraton, calculate the cosine and sine terms
		  			//cosine_array[j] = cos(node[i].Va-node[p].Va);
		  			//sine_array[j] = sin(node[i].Va-node[p].Va);
		  			cordic(node[i].Va-node[p].Va,&cosine_array[j],&sine_array[j]);
		  		}
		  		deltaP[i] += - node[p].Vm*(edge[j].G*cosine_array[j] + edge[j].B*sine_array[j]);
		  		if (node[i].type<2) // calculate Q for PQ buses
		  			deltaQ[i] += - node[p].Vm*(edge[j].G*sine_array[j] - edge[j].B*cosine_array[j]);
		  	}
		 }

		 // Get max P and max Q
		 if (fabs(deltaP[i]*node[i].Vm) > maxDeltaP)
			 maxDeltaP = fabs(deltaP[i]*node[i].Vm);
		 if (fabs(deltaQ[i]*node[i].Vm) > maxDeltaQ)
			 maxDeltaQ = fabs(deltaQ[i]*node[i].Vm);
	}
	para->maxDeltaP[local_id]=maxDeltaP;
	para->maxDeltaQ[local_id]=maxDeltaQ;
    return nullptr;
}

void rhs_deltaP_par(struct Node_PF *node,struct Edge_PF *edge,unsigned  *off,int node_num,
		double *deltaP,double *deltaQ,double *cosine_array,double *sine_array,
		bool tri_flag,double *maxDeltaP,double *maxDeltaQ){
	struct Rhs_Para_Mt rhs_para_mt[16];
	double l_maxDeltaP[16];
	double l_maxDeltaQ[16];
	pthread_t rhs_p_arr[16];
	for(int tid=0;tid<RHS_THREAD_NUM;tid++){
		rhs_para_mt[tid].node=node;
		rhs_para_mt[tid].edge=edge;
		rhs_para_mt[tid].off=off;
		rhs_para_mt[tid].node_num=node_num;
		rhs_para_mt[tid].deltaP=deltaP;
		rhs_para_mt[tid].deltaQ=deltaQ;
		rhs_para_mt[tid].maxDeltaP=l_maxDeltaP;
		rhs_para_mt[tid].maxDeltaQ=l_maxDeltaQ;
		rhs_para_mt[tid].cosine_array=cosine_array;
		rhs_para_mt[tid].sine_array=sine_array;
		rhs_para_mt[tid].tri_flag=tri_flag;
		rhs_para_mt[tid].local_id=tid;
		rhs_para_mt[tid].local_size=RHS_THREAD_NUM;
	}
	for(int tid=0;tid<RHS_THREAD_NUM;tid++){
		int err=pthread_create(&rhs_p_arr[tid],NULL,thr_rhs_deltaP,&rhs_para_mt[tid]);
		  		//if(err!=0) {printf("create thread:%d failed",tid);exit(-1);}
	}
	for(int tid=0;tid<RHS_THREAD_NUM;tid++){
		pthread_join(rhs_p_arr[tid],NULL);
	}
}

void rhs_deltaQ_par(struct Node_PF *node,struct Edge_PF *edge,unsigned  *off,int node_num,
		double *deltaP,double *deltaQ,double *cosine_array,double *sine_array,
		bool tri_flag,double *maxDeltaP,double *maxDeltaQ){
	struct Rhs_Para_Mt rhs_para_mt[16];
	double l_maxDeltaP[16];
	double l_maxDeltaQ[16];
	pthread_t rhs_p_arr[16];
	for(int tid=0;tid<RHS_THREAD_NUM;tid++){
		rhs_para_mt[tid].node=node;
		rhs_para_mt[tid].edge=edge;
		rhs_para_mt[tid].off=off;
		rhs_para_mt[tid].node_num=node_num;
		rhs_para_mt[tid].deltaP=deltaP;
		rhs_para_mt[tid].deltaQ=deltaQ;
		rhs_para_mt[tid].maxDeltaP=l_maxDeltaP;
		rhs_para_mt[tid].maxDeltaQ=l_maxDeltaQ;
		rhs_para_mt[tid].cosine_array=cosine_array;
		rhs_para_mt[tid].sine_array=sine_array;
		rhs_para_mt[tid].tri_flag=tri_flag;
		rhs_para_mt[tid].local_id=tid;
		rhs_para_mt[tid].local_size=RHS_THREAD_NUM;
	}
	for(int tid=0;tid<RHS_THREAD_NUM;tid++){
		int err=pthread_create(&rhs_p_arr[tid],NULL,thr_rhs_deltaQ,&rhs_para_mt[tid]);
		  		//if(err!=0) {printf("create thread:%d failed",tid);exit(-1);}
	}
	for(int tid=0;tid<RHS_THREAD_NUM;tid++){
		pthread_join(rhs_p_arr[tid],NULL);
	}
	for(int tid=0;tid<RHS_THREAD_NUM;tid++){
		if(*maxDeltaP<rhs_para_mt[tid].maxDeltaP[tid])
			*maxDeltaP=rhs_para_mt[tid].maxDeltaP[tid];
		if(*maxDeltaQ<rhs_para_mt[tid].maxDeltaQ[tid])
			*maxDeltaQ=rhs_para_mt[tid].maxDeltaQ[tid];
	}
}

void pf_iter_MT(){
	struct timeval start,finish;
	double cost_time;
	double Ybus_time=0,fac_time=0,rhs_time=0,back_for_time=0;
	gettimeofday(&start, 0);//计时开始
	struct Graph_PF_IN graph_in=pf_load("Case10790");
	struct Ret_Ybus ret_Ybus=pf_Ybus_par(graph_in);
	gettimeofday(&finish, 0);//计时结束
	cost_time = 1000000 * (finish.tv_sec - start.tv_sec) + finish.tv_usec - start.tv_usec;
	printf("pf_Ybus:cost the time is : %lf s.\n",cost_time/(double)1000000);
	Ybus_time=cost_time/(double)1000000;

	gettimeofday(&start, 0);//计时开始
	ICktSo nicslu_p = nullptr,nicslu_pp = nullptr;
	pf_fac_par(ret_Ybus.Bp_x,ret_Ybus.Bp_i,ret_Ybus.Bp_p,
			ret_Ybus.Bpp_x,ret_Ybus.Bpp_i,ret_Ybus.Bpp_p,
			ret_Ybus.graph.node_num,&nicslu_p,&nicslu_pp);
	gettimeofday(&finish, 0);//计时结束
	cost_time = 1000000 * (finish.tv_sec - start.tv_sec) + finish.tv_usec - start.tv_usec;
	printf("pf_fac:cost the time is : %lf s.\n",cost_time/(double)1000000);
	fac_time=cost_time/(double)1000000;

	struct Node_PF *node=ret_Ybus.graph.node;
	struct Edge_PF *edge=ret_Ybus.graph.edge;
	unsigned  *off=ret_Ybus.graph.off;
	int node_num=ret_Ybus.graph.node_num;
	int edge_num=ret_Ybus.graph.edge_num;

	double *deltaP=(double *)malloc(node_num*sizeof(double));
	double *deltaQ=(double *)malloc(node_num*sizeof(double));
	double *cosine_array=(double *)malloc(edge_num*sizeof(double));
	double *sine_array=(double *)malloc(edge_num*sizeof(double));

	int iter;
	for (iter=0; iter<10; ++iter){
	  	double maxDeltaP=0;
	  	double maxDeltaQ=0;

	  	gettimeofday(&start, 0);//计时开始
	  	{
	  		if(iter==0)
	  			rhs_deltaP_par(node,edge,off,node_num,
	  				deltaP,deltaQ,cosine_array,sine_array,
	  				true,&maxDeltaP,&maxDeltaQ);
	  		else
	  			rhs_deltaP_par(node,edge,off,node_num,
	  				deltaP,deltaQ,cosine_array,sine_array,
	  				false,&maxDeltaP,&maxDeltaQ);

	  	}
	  	// Calculate deltaP and deltaQ to update Va
	  	/*for (int i=0; i<node_num; ++i){
	  		deltaP[i] = 0;
	  		deltaQ[i] = 0;
	  		if (node[i].type!=3){
	  			deltaP[i] = node[i].P;
	  			if (node[i].type<2)
	  				deltaQ[i] = node[i].Q;
	  		}

	  		// Calculate network injections
	  		if (node[i].type!=3){
	  			for (int j=off[i]; j<off[i+1]; ++j){
	  				int p=edge[j].src;
	  				if (iter == 0){
	  					//if this is the first iteraton, calculate the cosine and sine terms
	  					cosine_array[j] = cos(node[i].Va-node[p].Va);
	  					sine_array[j] = sin(node[i].Va-node[p].Va);
	  				}
	  				deltaP[i] += - node[i].Vm*node[p].Vm*(edge[j].G*cosine_array[j] + edge[j].B*sine_array[j]);
	  				if (node[i].type<2) // calculate Q for PQ buses
	  					deltaQ[i] += - node[i].Vm*node[p].Vm*(edge[j].G*sine_array[j] - edge[j].B*cosine_array[j]);
	  			}
	  		}

	  		// Get max P and max Q
	  		if (fabs(deltaP[i]) > maxDeltaP)
	  			maxDeltaP = fabs(deltaP[i]);
	  		if (fabs(deltaQ[i]) > maxDeltaQ)
	  			maxDeltaQ = fabs(deltaQ[i]);

	  		deltaP[i] /= node[i].Vm;
	  		deltaQ[i] /= node[i].Vm;
	  	}*/
	  	gettimeofday(&finish, 0);//计时结束
	  	cost_time = 1000000 * (finish.tv_sec - start.tv_sec) + finish.tv_usec - start.tv_usec;
	  	printf("deltaP:cost the time is : %lf s.\n",cost_time/(double)1000000);
	  	rhs_time+=cost_time/(double)1000000;

	  	std::cout << "Iter" << iter << " Updated for Va, maxDeltaP, " << maxDeltaP << ", maxDeltaQ, " << maxDeltaQ << std::endl;

	      // Decide if the program converges
	  		/*if ( maxDeltaP < 0.001 && maxDeltaQ < 0.001){
	  			break;
	  		}*/

	      //reset the max mismatche
	      maxDeltaP = 0;
	      maxDeltaQ = 0;

	      gettimeofday(&start, 0);//计时开始
	      // solve for deltaVa
	      // std::cout << " ======================== Solve for V angle ========================"<<std::endl;
	      CKTSO_Solve(nicslu_p, deltaP,deltaP,false,false);
	      // Update V angle (Va)
	      for (int i=0; i<node_num; ++i){
	  			node[i].Va -= deltaP[i];
	  			if(node[i].Va>=PI){
	  				node[i].Va-=2*PI;
	  			}
	  			if(node[i].Va<-PI){
	  				node[i].Va+=2*PI;
	  			}
	        //std::cout << "Iter" << iter << ", Va" << i << ", " << Va[i]<< ", deltaVa" << i << ", " << deltaP[i]<< std::endl;
	      }
	      gettimeofday(&finish, 0);//计时结束
	      cost_time = 1000000 * (finish.tv_sec - start.tv_sec) + finish.tv_usec - start.tv_usec;
	      printf("back_forward and update Va:cost the time is : %lf s.\n",cost_time/(double)1000000);
	      back_for_time+=cost_time/(double)1000000;
	      // std::cout << " ======================== Update deltaP and deltaQ ========================"<<std::endl;

	      gettimeofday(&start, 0);//计时开始
	      {
	    	  rhs_deltaQ_par(node,edge,off,node_num,
	    			  deltaP,deltaQ,cosine_array,sine_array,
					  true,&maxDeltaP,&maxDeltaQ);
	      }
	      /*// Calculate deltaP and deltaQ to update Vm
	      for (int i=0; i<node_num; ++i){
	      	  deltaP[i] = 0;
	      	  deltaQ[i] = 0;

	      	  if (node[i].type!=3){
	      		  deltaP[i] = node[i].P;
	      		  if (node[i].type<2)
	      			  deltaQ[i] = node[i].Q;
	      	  }

	      	  // Calculate network injections
	      	  if (node[i].type!=3){
	      		  for (int j=off[i]; j<off[i+1]; ++j){
	      			  int p=edge[j].src;
	      			  cosine_array[j] = cos(node[i].Va-node[p].Va);
	      			  sine_array[j] = sin(node[i].Va-node[p].Va);
	      			  deltaP[i] += - node[i].Vm*node[p].Vm*(edge[j].G*cosine_array[j] + edge[j].B*sine_array[j]);
	      			  if (node[i].type<2) // calculate Q for PQ buses
	      				  deltaQ[i] += - node[i].Vm*node[p].Vm*(edge[j].G*sine_array[j] - edge[j].B*cosine_array[j]);
	      		  }
	      	  }

	      	  // Get max P and max Q
	      	  if (fabs(deltaP[i]) > maxDeltaP)
	      		  maxDeltaP = fabs(deltaP[i]);
	      	  	  if (fabs(deltaQ[i]) > maxDeltaQ)
	      	  		  maxDeltaQ = fabs(deltaQ[i]);

	      	  deltaP[i] /= node[i].Vm;
	      	  deltaQ[i] /= node[i].Vm;
	      }*/
	      gettimeofday(&finish, 0);//计时结束
	      cost_time = 1000000 * (finish.tv_sec - start.tv_sec) + finish.tv_usec - start.tv_usec;
	      printf("deltaQ:cost the time is : %lf s.\n",cost_time/(double)1000000);
	      rhs_time+=cost_time/(double)1000000;
	      std::cout << "Iter" << iter << " Updated for Vm, maxDeltaP, " << maxDeltaP << ", maxDeltaQ, " << maxDeltaQ << std::endl;

	      // Decide if the program converges
	      if ( maxDeltaP < 0.001 && maxDeltaQ < 0.001){
	  			break;
	      }

	      gettimeofday(&start, 0);//计时开始
	      // solve for deltaVm
	      // std::cout << " ======================== Solve for V magnitude ========================"<<std::endl;
	      CKTSO_Solve(nicslu_pp,deltaQ ,deltaQ,false,false);

	      // Update V magnitude (Vm)
	      for (int i=0; i<node_num; ++i)
	      {
	    	  node[i].Vm -= deltaQ[i];
	    	  //std::cout << "Iter" << iter << ", Vm" << i << ", " << Vm[i]<< ", deltaVm" << i << ", " << deltaQ[i]<< std::endl;
	      }
	      gettimeofday(&finish, 0);//计时结束
	      cost_time = 1000000 * (finish.tv_sec - start.tv_sec) + finish.tv_usec - start.tv_usec;
	      printf("back_forward and update Vm:cost the time is : %lf s.\n",cost_time/(double)1000000);
	      back_for_time+=cost_time/(double)1000000;
	}
	printf("Ybus_time:%lf\tfac_time:%lf\trhs_time:%lf\tback_for_time:%lf\n",Ybus_time,fac_time,rhs_time,back_for_time);
	printf("rhs_avg:%lf\tback_for_avg:%lf\n",rhs_time/(iter+1),back_for_time/(iter+1));
}


/*struct Node_PF_In* load_node_pf(int lines){
    struct Node_PF_In *node=(struct Node_PF_In *)malloc(lines*sizeof(struct Node_PF_In));
    char nodefile[64]="/home/zyh/power_data/pf_data/Case10790_Nodeinfo_YT.csv";
    int i;
    //char edgefile[64]="Case10790_Edgeinfo_YT.csv";
    FILE *fp_node=fopen(nodefile,"r");

    char temp[512];
    fgets(temp,512,fp_node);

    for(i = 0 ;i < lines; i++){
        fscanf(fp_node,"%d",&node[i].bus_id);fseek(fp_node, 1L, SEEK_CUR);
        fscanf(fp_node,"%d",&node[i].bus_name);fseek(fp_node, 1L, SEEK_CUR);
        fscanf(fp_node,"%d",&node[i].area);fseek(fp_node, 1L, SEEK_CUR);
        fscanf(fp_node,"%d",&node[i].loss_zone);fseek(fp_node, 1L, SEEK_CUR);
        fscanf(fp_node,"%d",&node[i].type);fseek(fp_node, 1L, SEEK_CUR);

        fscanf(fp_node,"%lf",&node[i].voltage);fseek(fp_node, 1L, SEEK_CUR);
        fscanf(fp_node,"%lf",&node[i].angle);fseek(fp_node, 1L, SEEK_CUR);
        fscanf(fp_node,"%lf",&node[i].load_P);fseek(fp_node, 1L, SEEK_CUR);
        fscanf(fp_node,"%lf",&node[i].load_Q);fseek(fp_node, 1L, SEEK_CUR);
        fscanf(fp_node,"%lf",&node[i].generation_p);fseek(fp_node, 1L, SEEK_CUR);
        fscanf(fp_node,"%lf",&node[i].generation_Q);fseek(fp_node, 1L, SEEK_CUR);

        fscanf(fp_node,"%lf",&node[i].base_kV);fseek(fp_node, 1L, SEEK_CUR);
        fscanf(fp_node,"%lf",&node[i].desired_volts);fseek(fp_node, 1L, SEEK_CUR);
        fscanf(fp_node,"%lf",&node[i].MAX_Q);fseek(fp_node, 1L, SEEK_CUR);
        fscanf(fp_node,"%lf",&node[i].Min_Q);fseek(fp_node, 1L, SEEK_CUR);
        fscanf(fp_node,"%lf",&node[i].G);fseek(fp_node, 1L, SEEK_CUR);
        fscanf(fp_node,"%lf",&node[i].B);fseek(fp_node, 1L, SEEK_CUR);

        fscanf(fp_node,"%d",&node[i].control_bus_number);fseek(fp_node, 1L, SEEK_CUR);
    }
    fclose(fp_node);

    printf("\nlines:%d\n",lines);
    for(i=0;i<lines;i++){
        if(i>=5&&linnodefilees-i>5) continue;
        printf("%d\t%d\t%d\t%d\t%d\t",node[i].bus_id,node[i].bus_name,node[i].area,node[i].loss_zone,node[i].type);
        printf("%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t",node[i].voltage,node[i].angle,node[i].load_P,node[i].load_Q,node[i].generation_p,node[i].generation_Q);
        printf("%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t",node[i].base_kV,node[i].desired_volts,node[i].MAX_Q,node[i].Min_Q,node[i].G,node[i].B);
        printf("%d\n",node[i].control_bus_number);
    }
    return node;
}

const double pi_value = 3.141592653589793;
double rand_Vm(){
    double ret;
    if(rand()%2==0){
        ret=1.0+(0.2*rand()/(RAif(j<10) ND_MAX+1.0));
    }
    else{
        ret=1.0-(0.2*rand()/(RAND_MAX+1.0));
    }
    return ret;
}
double rand_Va(){
    double ret;
    if(rand()%2==0){
        ret=0.1000*pi_value*rand()/(RAND_MAX+1.0);
    }
    else{
        ret=-0.1000*pi_value*rand()/(RAND_MAX+1.0);
    }if(j<10)
    return ret;
}
struct Node_PF_Out* get_node_out(struct Node_PF_In *node_in,int lines){
    //srand((int)time(0));
    srand(23333);

    struct Node_PF_Out *node_out=(struct Node_PF_Out *)malloc((lines+32)*sizeof(struct Node_PF_Out));
    int i;
    for(i=0;i<lines;i++){
        node_out[i].type=(unsigned char)node_in[i].type;
        //P:generator_p Q:generator_Q Pl:LoadP Ql:loadQ
        if(node_out[i].type==0||node_out[i].type==1||node_out[i].type==2){
            node_out[i].P=node_in[i].generation_p-node_in[i].load_P;
            node_out[i].Q=node_in[i].generation_Q-node_in[i].load_Q;
        }
        else{G
            node_out[i].P=0.0;//random
            node_out[i].Q=0.0;//random
        }
        if(node_out[i].type==0||node_out[i].type==1||node_out[i].type==2||node_out[i].type==3){
            node_out[i].Vm=rand_Vm();//1.0;
            node_out[i].Va=rand_Va();//0.0;
        }
        else{
            node_out[i].Vm=rand_Vm();//0.0;
            node_out[i].Va=rand_Va();//0.0;
        }
    }
Node_PF_Out
    printf("\nnode_num:%d\n",lines);
    for(i=0;i<lines;i++){
        if(i>=5&&lines-i>5) continue;
        printf("%d\t",node_out[i].type);
        printf("%lf\t%lf\t%lf\t%lf\n",node_out[i].P,node_out[i].Q,node_out[i].Vm,node_out[i].Va);
    }

    return node_out;
}



struct Edge_PF_In* load_edge_pf(int lines){
    struct Edge_PF_In *edge=(struct Edge_PF_In *)malloc(lines*sizeof(struct Edge_PF_In));
    char edgefile[64]="/home/zyh/power_data/pf_data/Case10790_Edgeinfo_YT.csv";
    int i;

    FILE *fp_edge=fopen(edgefile,"r");
    fseek(fp_edge, 18L, SEEK_SET);   // ŽÓÎÄŒþµÚ¶þÐÐ¿ªÊŒ¶ÁÈ¡

    char temp[512];
    fgets(temp,512,fp_edge);

    for(i = 0 ;i < lines; i++){
        fscanf(fp_edge,"%d",&edge[i].ap_bus);fseek(fp_edge, 1L, SEEK_CUR);
        fscanf(fp_edge,"%d",&edge[i].z_bus);fseek(fp_edge, 1L, SEEK_CUR);
        fscanf(fp_edge,"%d",&edge[i].area);fseek(fp_edge, 1L, SEEK_CUR);
        fscanf(fp_edge,"%d",&edge[i].zone);fseek(fp_edge, 1L, SEEK_CUR);
        fscanf(fp_edge,"%d",&edge[i].circuit);fseek(fp_edge, 1L, SEEK_CUR);
        fscanf(fp_edge,"%d",&edge[i].type);fseek(fp_edge, 1L, SEEK_CUR);

        fscanf(fp_edge,"%lf",&edge[i].R);fseek(fp_edge, 1L, SEEK_CUR);
        fscanf(fp_edge,"%lf",&edge[i].X);fseek(fp_edge, 1L, SEEK_CUR);
        fscanf(fp_edge,"%lf",&edge[i].B);fseek(fp_edge, 1L, SEEK_CUR);
        fscanf(fp_edge,"%lf",&edge[i].Line_Q1);fseek(fp_edge, 1L, SEEK_CUR);
        fscanf(fp_edge,"%lf",&edge[i].Line_Q2);fseek(fp_edge, 1L, SEEK_CUR);
        fscanf(fp_edge,"%lf",&edge[i].Line_Q3);fseek(fp_edge, 1L, SEEK_CUR);

        fscanf(fp_edge,"%d",&edge[i].control_bus);fseek(fp_edge, 1L, SEEK_CUR);
        fscanf(fp_edge,"%d",&edge[i].side);fseek(fp_edge, 1L, SEEK_CUR);
        fscanf(fp_edge,"%lf",&edge[i].transformer_final_turns_ratio);fseek(fp_edge, 1L, SEEK_CUR);
        fscanf(fp_edge,"%lf",&edge[i].transformer_final_angle);fseek(fp_edge, 1L, SEEK_CUR);
        fscanf(fp_edge,"%lf",&edge[i].Min_tap);fseek(fp_edge, 1L, SEEK_CUR);
        fscanf(fp_edge,"%lf",&edge[i].Max_tap);fseek(fp_edge, 1L, SEEK_CUR);

        fscanf(fp_edge,"%d",&edge[i].step_size);fseek(fp_edge, 1L, SEEK_CUR);

        fscanf(fp_edge,"%lf",&edge[i].Min_volt);fseek(fp_edge, 1L, SEEK_CUR);
        fscanf(fp_edge,"%lf",&edge[i].Max_volt);fseek(fp_edge, 1L, SEEK_CUR);
    }
    fclose(fp_edge);

    printf("\n");
    for(i=0;i<lines;i++){
        if(i>=5&&lines-i>5) continue;
        printf("%d\t%d\t%d\t%d\t%d\t%d\t",edge[i].ap_bus,edge[i].z_bus,edge[i].area,edge[i].zone,edge[i].circuit,edge[i].type);
        printf("%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t",edge[i].R,edge[i].X,edge[i].B,edge[i].Line_Q1,edge[i].Line_Q2,edge[i].Line_Q3);
        printf("%d\t%d\t",edge[i].control_bus,edge[i].side);
        printf("%lf\t%lf\t%lf\t%lf\t",edge[i].transformer_final_turns_ratio,edge[i].transformer_final_angle,edge[i].Min_tap,edge[i].Max_tap);
        printf("%d\t",edge[i].step_size);
        printf("%lf\t%lf\n",edge[i].Min_volt,edge[i].Max_volt);
    }
    return edge;
}

struct Edge_PF_In_Pre* pre_edge_pf(struct Edge_PF_In *edge_in,int lines,int *new_lines){
    struct Edge_PF_In *edge_sort;
    struct Edge_PF_In_Pre *edge_pre;
    int i,k;
    edge_sort=(struct Edge_PF_In *)malloc(2*lines*sizeof(struct Edge_PF_In));
    for(i=0;i<lines;i++){
        edge_sort[i]=edge_in[i];
        edge_sort[i+lines]=edge_in[i];
        int temp=edge_sort[i+lines].ap_bus;
        edge_sort[i+lines].ap_bus=edge_sort[i+lines].z_bus;
        edge_sort[i+lines].z_bus=temp;
        edge_sort[i+lines].transformer_final_turns_ratio=-edge_sort[i+lines].transformer_final_turns_ratio;
    }
    for(int i = 1;i<2*lines;i++){
        struct Edge_PF_In temp = edge_sort[i];
        int j;
        for(j = i-1;j>=0;j--){
        //œ«ŽóÓÚtempµÄÊýÏòºóÒÆ¶¯Ò»²œ
            if(edge_sort[j].ap_bus>temp.ap_bus||(edge_sort[j].ap_bus==temp.ap_bus&&edge_sort[j].z_bus>temp.z_bus)){
                edge_sort[j+1] = edge_sort[j];//ŒÇÂŒjµÄÖµÒ²ŸÍÊÇtempÒª²åÈëµÄÎ»ÖÃ
            }else{
                break;
            }
        }
        edge_sort[j+1] = temp;
    }

    edge_pre=(struct Edge_PF_In_Pre *)malloc(2*lines*sizeof(struct Edge_PF_In_Pre));
    int cur_ap_bus=0,cur_z_bus=0;
    k=0;
    for(int i = 1;i<2*lines;i++){
        if(cur_ap_bus!=edge_sort[i].ap_bus||cur_z_bus!=edge_sort[i].z_bus){
            ++k;
            edge_pre[k-1].from=edge_sort[i].ap_bus-1;
            edge_pre[k-1].to=edge_sort[i].z_bus-1;
            edge_pre[k-1].G=0.0;
            edge_pre[k-1].B=0.0;
            edge_pre[k-1].hB=0.0;
            edge_pre[k-1].K=0.0;
            edge_pre[k-1].Kcount=0;
            edge_pre[k-1].BIJ=0.0;

            cur_ap_bus=edge_sort[i].ap_bus;
            cur_z_bus=edge_sort[i].z_bus;
        }
        edge_pre[k-1].G+=edge_sort[i].R/(edge_sort[i].R*edge_sort[i].R+edge_sort[i].X*edge_sort[i].X);
        edge_pre[k-1].B+=edge_sort[i].X/(edge_sort[i].R*edge_sort[i].R+edge_sort[i].X*edge_sort[i].X);
        edge_pre[k-1].hB+=edge_sort[i].B;
        edge_pre[k-1].K+=edge_sort[i].transformer_final_turns_ratio;
        edge_pre[k-1].Kcount+=1;
        edge_pre[k-1].BIJ+=1.0/edge_sort[i].X;
    }
    *new_lines=k;

    printf("\nlines:%d\tnew_lines:%d\n",lines,*new_lines);
    for(i=0;i<*new_lines;i++){
        if(i>=5&&*new_lines-i>20) continue;
        printf("%d\t%d\t",edge_pre[i].from,edge_pre[i].to);
        printf("%lf\t%lf\t%lf\t%lf\t",edge_pre[i].G,edge_pre[i].B,edge_pre[i].hB,edge_pre[i].K);
        printf("%d\t",edge_pre[i].Kcount);
        printf("%lf\n",edge_pre[i].BIJ);
    }

    return edge_pre;
}

unsigned short *get_edge_off(struct Edge_PF_In_Pre *edge_pre, int node_num, int edge_num){
    unsigned short *edge_off=(unsigned short *)malloc((node_num+1+32)*sizeof(unsigned short));
    int cur=-1;
    int i;
    for(i=0;i<edge_num;i++){
        while(edge_pre[i].from!=cur){
            assert(edge_pre[i].from-cur==1);
            cur++;
            edge_off[cur]=(unsigned short)i;
        }
    }
    edge_off[node_num]=(unsigned short)edge_num;

    printf("cur:%d\n",cur);
    printf("\nnode_num:%d\tedge_num:%d\n",node_num,edge_num);
    for(i=0;i<node_num+1;i++){
        if(i>=5&&node_num+1-i>10) continue;
        printf("%hu\n",edge_off[i]);
    }

    return edge_off;
}

struct Edge_PF_Out *get_edge_info(struct Edge_PF_In_Pre *edge_pre, int edge_num){
    struct Edge_PF_Out *edge_info=(struct Edge_PF_Out *)malloc((edge_num+32)*sizeof(struct Edge_PF_Out));
    int i;
    for(i=0;i<edge_num;i++){
        edge_info[i].src=(unsigned short)edge_pre[i].to;
        if(edge_pre[i].K==0){
            edge_info[i].G=-edge_pre[i].G;
            edge_info[i].B=edge_pre[i].B;
        }
        else{
            double tap_ratio = edge_pre[i].K/edge_pre[i].Kcount;
            edge_info[i].G=-edge_pre[i].G/tap_ratio;
            edge_info[i].B=edge_pre[i].B/tap_ratio;
        }
    }

    printf("\nedge_num:%d\n",edge_num);
    for(i=0;i<edge_num;i++){
        if(i>=5&&edge_num-i>20) continue;
        printf("%d\t",edge_info[i].src);
        printf("%lf\t%lf\n",edge_info[i].G,edge_info[i].B);
    }

    return edge_info;
}



#include <sys/time.h>
void run_pf(struct Node_PF_Out *node,unsigned short *edge_off,struct Edge_PF_Out *edge_info,int node_num,int edge_num){
    int i,j;
    double *deltaP=(double *)malloc(node_num*sizeof(double));
    double *deltaQ=(double *)malloc(node_num*sizeof(double));
    double maxDeltaP=0.0,maxDeltaQ=0.0;

    struct timeval start,finish;
    double cost_time;
    gettimeofday(&start, 0);//计时开始
    for (i=0; i<node_num; ++i){
  			deltaP[i] = 0;
  			deltaQ[i] = 0;

  			if (node[i].type!=3){
  				deltaP[i] = node[i].P;
  				if (node[i].type<2)
  					deltaQ[i] = node[i].Q;
  			}

  			// Calculate network injections
  			if (node[i].type!=3){
  				for (j=edge_off[i]; j<edge_off[i+1]; ++j){
  					unsigned short p=edge_info[j].src;
                    //if this is the first iteraton, calculate the cosine and sine terms
                    double cosine_array = cos(node[i].Va-node[p].Va);
                    double sine_array = sin(node[i].Va-node[p].Va);
                    deltaP[i] += - node[i].Vm*node[p].Vm*(edge_info[j].G*cosine_array + edge_info[j].B*sine_array);
  					if (node[i].type<2) // calculate Q for PQ buses
                        deltaQ[i] += - node[i].Vm*node[p].Vm*(edge_info[j].G*sine_array - edge_info[j].B*cosine_array);
                    }
                }

                // Get max P and max Q
                if (fabs(deltaP[i]) > maxDeltaP)
                    maxDeltaP = fabs(deltaP[i]);
                if (fabs(deltaQ[i]) > maxDeltaQ)
                    maxDeltaQ = fabs(deltaQ[i]);

                deltaP[i] /= node[i].Vm;
                deltaQ[i] /= node[i].Vm;

  			if(node[i].type==3){
                deltaP[i]=0.0;
            }
            if(node[i].type==2||node[i].type==3){
                deltaQ[i]=0.0;
            }
            if(i<10||node_num-i<=10){
                printf("i:%d\tdeltaP:%lf\tdeltaQ:%lf\n",i,deltaP[i],deltaQ[i]);
            }
  		}
    	gettimeofday(&finish, 0);//计时结束
        cost_time = 1000000 * (finish.tv_sec - start.tv_sec) + finish.tv_usec - start.tv_usec;
        printf("run_pf:cost the time is : %lf s.\n",cost_time/(double)1000000);

  		//printf("maxDeltaP:%lf\tmaxDeltaQ:%lf\n",maxDeltaP,maxDeltaQ);
}

void run_pf2(struct Node_PF_Out *node,unsigned short *edge_off,struct Edge_PF_Out *edge_info,int node_num,int edge_num){
    int i,j;
    double *deltaP=(double *)malloc(node_num*sizeof(double));
    double *deltaQ=(double *)malloc(node_num*sizeof(double));
    double maxDeltaP=0.0,maxDeltaQ=0.0;

    struct timeval start,finish;
    double cost_time;
    gettimeofday(&start, 0);//计时开始
    for (i=0; i<node_num; ++i){
  			deltaP[i] = node[i].P;
  			deltaQ[i] = node[i].Q;

  			// Calculate network injections
  			//if (node[i].type!=3){
  			for (j=edge_off[i]; j<edge_off[i+1]; ++j){
  				unsigned short p=edge_info[j].src;
                //if this is the first iteraton, calculate the cosine and sine terms
                double cosine_array = cos(node[i].Va-node[p].Va);
                double sine_array = sin(node[i].Va-node[p].Va);
                deltaP[i] += - node[i].Vm*node[p].Vm*(edge_info[j].G*cosine_array + edge_info[j].B*sine_array);
                deltaQ[i] += - node[i].Vm*node[p].Vm*(edge_info[j].G*sine_array - edge_info[j].B*cosine_array);
            }

            // Get max P and max Q
            if (fabs(deltaP[i]) > maxDeltaP)
                maxDeltaP = fabs(deltaP[i]);
            if (fabs(deltaQ[i]) > maxDeltaQ)
                maxDeltaQ = fabs(deltaQ[i]);

            deltaP[i] /= node[i].Vm;
            deltaQ[i] /= node[i].Vm;

  			if(node[i].type==3){
                deltaP[i]=0.0;
            }
            if(node[i].type==2||node[i].type==3){
                deltaQ[i]=0.0;
            }
            if(i<10||node_num-i<=10){
                printf("i:%d\tdeltaP:%lf\tdeltaQ:%lf\n",i,deltaP[i],deltaQ[i]);
            }
  		}
    	gettimeofday(&finish, 0);//计时结束
        cost_time = 1000000 * (finish.tv_sec - start.tv_sec) + finish.tv_usec - start.tv_usec;
        printf("run_pf:cost the time is : %lf s.\n",cost_time/(double)1000000);

  		//printf("maxDeltaP:%lf\tmaxDeltaQ:%lf\n",maxDeltaP,maxDeltaQ);
}

void pf_debug(struct Node_PF_Out *node,unsigned short *edge_off,struct Edge_PF_Out *edge_info,int node_num,int edge_num,struct PF_PQ *pf_pq){
    int i,j;
    double *deltaP=(double *)malloc(node_num*sizeof(double));
    double *deltaQ=(double *)malloc(node_num*sizeof(double));
    double maxDeltaP=0.0,maxDeltaQ=0.0;

    for (i=0; i<node_num; ++i){
  			deltaP[i] = 0;
  			deltaQ[i] = 0;

  			if (node[i].type!=3){
  				deltaP[i] = node[i].P;
  				if (node[i].type<2)
  					deltaQ[i] = node[i].Q;
  			}

  			// Calculate network injections
  			if (node[i].type!=3){
  				for (j=edge_off[i]; j<edge_off[i+1]; ++j){
  					unsigned short p=edge_info[j].src;
                    //if this is the first iteraton, calculate the cosine and sine terms
                    double cosine_array = cos(node[i].Va-node[p].Va);
                    double sine_array = sin(node[i].Va-node[p].Va);
                    deltaP[i] += - node[i].Vm*node[p].Vm*(edge_info[j].G*cosine_array + edge_info[j].B*sine_array);
  					if (node[i].type<2) // calculate Q for PQ buses
                        deltaQ[i] += - node[i].Vm*node[p].Vm*(edge_info[j].G*sine_array - edge_info[j].B*cosine_array);
                    }Node_PF_Out
                }

                // Get max P and max Q
                if (fabs(deltaP[i]) > maxDeltaP)
                    maxDeltaP = fabs(deltaP[i]);
                if (fabs(deltaQ[i]) > maxDeltaQ)
                    maxDeltaQ = fabs(deltaQ[i]);

                deltaP[i] /= node[i].Vm;
                deltaQ[i] /= node[i].Vm;

  			if(node[i].type==3){
                deltaP[i]=0.0;
            }
            if(node[i].type==2||node[i].type==3){
                deltaQ[i]=0.0;
            }
            if(i<10||node_numNode_PF_Out-i<=10){
                printf("host:i:%d\tdeltaP:%lf\tdeltaQ:%lf\n",i,deltaP[i],deltaQ[i]);
                printf("device:i:%d\tdeltaP:%lf\tdeltaQ:%lf\n",i,pf_pq[i].deltaP
                		,pf_pq[i].deltaQ);
            }
  		}

  		printf("maxDeltaP:%lf\tmaxDeltaQ:%lf\n",maxDeltaP,maxDeltaQ);
}


void print_double(double in){
    unsigned char *p=(unsigned char *)&in;
    int i;
    for(i=0;i<8;i++){
        printf("%02X",(unsigned char)p[i]);
    }
}Node_PF_Out
void print_int(unsigned int in){
    unsigned char *p=(unsigned char *)&in;
    int i;
    for(i=0;i<4;i++){
        printf("%02X",(unsigned char)p[i]);
    }
}
void print_short(unsigned short in){
    unsigned char *p=(unsigned char *)&in;
    int i;
    for(i=0;i<2;i++){
        printf("%02X",(unsigned char)p[i]);
    }
}
void print_char(unsigned char in){
    unsigned char *p=(unsigned char *)&in;
    int i;
    for(i=0;i<1;i++){
        printf("%02X",(unsigned Node_PF_Outchar)p[i]);
    }
}
void test(){
    union Node{
        long long lv;
        double dv;
    };

    int type=0x02;
    union Node P,Q,Vm,Va;
    P.lv=0x3FEFFBE76C8B4396;
    Q.lv=0x3FE4D42C3C9EECC0;
    Vm.lv=0x3FF0000000000000;
    Va.lv=0x0000000000000000;

    printf("%d %lf %lf %lf %lf\n",type,P.dv,Q.dv,Vm.dv,Va.dv);
    printf("%02X%016I64X%016I64X%016I64X%016I64X\n",type,P.dv,Q.dv,Vm.dv,Va.dv);

    double x1=3.14;
    int x2=6837823;Node_PF_Out
    short x3=23232;Node_PF_Out
    char x4=134;
    print_double(x1);printf("\t%016I64X\n",x1);
    print_int(x2);printf("\t%08X\n",x2);
    print_short(x3);printf("\t%04X\n",(short)x3);
    print_char(x4);printf("\t%02X\n",(unsigned char)x4);
}
*/
