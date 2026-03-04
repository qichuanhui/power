#include "hls_math.h"
#include <ap_int.h>
#include <string.h>
#include <stdio.h>


struct Node_SE_V{
	double Vm,Va;
};
struct SE_H_r_PQ{
	double H_r_P;
	double H_r_Q;
};



struct Node_SE_CONST{//32byte
    double P, Q, Ri_V, M_Vm;
};
struct Edge_SE_CONST{//128byte
	int src,Kcount;
	double K;
    double G, B, hB, BIJ;
    double M_P_TLPF, M_Q_TLPF, Ri_eP, Ri_eQ;
    double M_P_TLPF_reverse,M_Q_TLPF_reverse,Ri_eP_reverse,Ri_eQ_reverse;
    //double filling1,filling2;
};

#define MAX_NODE_NUM (6*1024)
#define MAX_EDGE_NUM (8*1024)


union Int2_or_Double{
	int data_i[2];
	double data_d;
};
union Int64_or_Double{
	long data_i;
	double data_d;
};

#define DATAWIDTH 512
typedef ap_int<DATAWIDTH> TYPE;
typedef struct Block_Type_D{
	double data[8];
}TYPE_D;
typedef struct Block_Type_I{
	int data[16];
}TYPE_I;

void init_bank_memcpy(const TYPE* const  node_const,struct Node_SE_CONST*  l_node_const,
				const TYPE* const  edge_off,int*  l_off_0,int*  l_off_1,
                const TYPE* const  edge_const,struct Edge_SE_CONST*  l_edge_const,
				const int begin,
				const int end
				){
	//#pragma HLS INLINE
	#pragma HLS allocation instances=init_bank_memcpy limit=4 function
	const int node_num=end-begin;

	for(int i=0;i<node_num;i+=8){
		#pragma HLS PIPELINE II=1

		TYPE tmp[4];
		tmp[0]=node_const[i/2];
		tmp[1]=node_const[i/2+1];
		tmp[2]=node_const[i/2+2];
		tmp[3]=node_const[i/2+3];
		for(int j=0;j<8;j++){
			#pragma HLS unroll factor=8

			long data_l_P,data_l_Q,data_l_Ri_V,data_l_M_Vm;
			data_l_P=tmp[j/2].range(256*(j%2)+63,256*(j%2)).to_long();
			data_l_Q=tmp[j/2].range(256*(j%2)+127,256*(j%2)+64).to_long();
			data_l_Ri_V=tmp[j/2].range(256*(j%2)+191,256*(j%2)+128).to_long();
			data_l_M_Vm=tmp[j/2].range(256*(j%2)+255,256*(j%2)+192).to_long();

			union Int64_or_Double data_tmp_P,data_tmp_Q,data_tmp_Ri_V,data_tmp_M_Vm;
			data_tmp_P.data_i=data_l_P;
			data_tmp_Q.data_i=data_l_Q;
			data_tmp_Ri_V.data_i=data_l_Ri_V;
			data_tmp_M_Vm.data_i=data_l_M_Vm;

			double data_P,data_Q,data_Ri_V,data_M_Vm;
			data_P=data_tmp_P.data_d;
			data_Q=data_tmp_Q.data_d;
			data_Ri_V=data_tmp_Ri_V.data_d;
			data_M_Vm=data_tmp_M_Vm.data_d;

			l_node_const[i+j].P=data_P;
			l_node_const[i+j].Q=data_Q;
			l_node_const[i+j].Ri_V=data_Ri_V;
			l_node_const[i+j].M_Vm=data_M_Vm;
		}
	}

	const int off_length=node_num+1;
	for(int i=0;i<off_length;i+=16){
		#pragma HLS PIPELINE II=1

		TYPE tmp=edge_off[i/16];
		for(int j=0;j<8;j++){
			#pragma HLS unroll factor=8

			int data;
			data=tmp.range(32*j+31,32*j).to_int();

			l_off_0[i+j]=data;
			l_off_1[i+j]=data;
		}
		for(int j=8;j<16;j++){
			#pragma HLS unroll factor=8

			int data;
			data=tmp.range(32*j+31,32*j).to_int();

			l_off_0[i+j]=data;
			l_off_1[i+j]=data;
		}
	}

	const int edge_num=l_off_0[node_num]-l_off_0[0];

	for(int i=0;i<edge_num;i+=16){
		#pragma HLS PIPELINE II=1

		TYPE tmp[32];
		for(int j=0;j<32;j++){
			#pragma HLS unroll factor=32
			tmp[j]=edge_const[i*2+j];
		}

		int src,Kcount;
		double K;
	    double G, B, hB, BIJ;
	    double M_P_TLPF, M_Q_TLPF, Ri_eP, Ri_eQ;
	    double M_P_TLPF_reverse,M_Q_TLPF_reverse,Ri_eP_reverse,Ri_eQ_reverse;
	    double filling1,filling2;

		double data[16];
		for(int j=0;j<16;j++){
			#pragma HLS unroll factor=16

			l_edge_const[i+j].src=tmp[2*j].range(31,0).to_int();
			l_edge_const[i+j].Kcount=tmp[2*j].range(63,32).to_int();

			long data_l_K,data_l_G,data_l_B,data_l_hB,data_l_BIJ;
			long data_l_M_P_TLPF,data_l_M_Q_TLPF,data_l_Ri_eP,data_l_Ri_eQ;
			long data_l_M_P_TLPF_reverse,data_l_M_Q_TLPF_reverse,data_l_Ri_eP_reverse,data_l_Ri_eQ_reverse;
			data_l_K=tmp[2*j].range(127,64).to_long();
			data_l_G=tmp[2*j].range(191,128).to_long();
			data_l_B=tmp[2*j].range(255,192).to_long();
			data_l_hB=tmp[2*j].range(319,256).to_long();
			data_l_BIJ=tmp[2*j].range(383,320).to_long();
			data_l_M_P_TLPF=tmp[2*j].range(447,384).to_long();
			data_l_M_Q_TLPF=tmp[2*j].range(511,448).to_long();

			data_l_Ri_eP=tmp[2*j+1].range(63,0).to_long();
			data_l_Ri_eQ=tmp[2*j+1].range(127,64).to_long();
			data_l_M_P_TLPF_reverse=tmp[2*j+1].range(191,128).to_long();
			data_l_M_Q_TLPF_reverse=tmp[2*j+1].range(255,192).to_long();
			data_l_Ri_eP_reverse=tmp[2*j+1].range(319,256).to_long();
			data_l_Ri_eQ_reverse=tmp[2*j+1].range(383,320).to_long();

			union Int64_or_Double data_tmp_K,data_tmp_G,data_tmp_B,data_tmp_hB,data_tmp_BIJ;
			union Int64_or_Double data_tmp_M_P_TLPF,data_tmp_M_Q_TLPF,data_tmp_Ri_eP,data_tmp_Ri_eQ;
			union Int64_or_Double data_tmp_M_P_TLPF_reverse,data_tmp_M_Q_TLPF_reverse,data_tmp_Ri_eP_reverse,data_tmp_Ri_eQ_reverse;
			data_tmp_K.data_i=data_l_K;
			data_tmp_G.data_i=data_l_G;
			data_tmp_B.data_i=data_l_B;
			data_tmp_hB.data_i=data_l_hB;
			data_tmp_BIJ.data_i=data_l_BIJ;
			data_tmp_M_P_TLPF.data_i=data_l_M_P_TLPF;
			data_tmp_M_Q_TLPF.data_i=data_l_M_Q_TLPF;

			data_tmp_Ri_eP.data_i=data_l_Ri_eP;
			data_tmp_Ri_eQ.data_i=data_l_Ri_eQ;
			data_tmp_M_P_TLPF_reverse.data_i=data_l_M_P_TLPF_reverse;
			data_tmp_M_Q_TLPF_reverse.data_i=data_l_M_Q_TLPF_reverse;
			data_tmp_Ri_eP_reverse.data_i=data_l_Ri_eP_reverse;
			data_tmp_Ri_eQ_reverse.data_i=data_l_Ri_eQ_reverse;


			double data_K,data_G,data_B,data_hB,data_BIJ;
			double data_M_P_TLPF,data_M_Q_TLPF,data_Ri_eP,data_Ri_eQ;
			double data_M_P_TLPF_reverse,data_M_Q_TLPF_reverse,data_Ri_eP_reverse,data_Ri_eQ_reverse;

			data_K=data_tmp_K.data_d;
			data_G=data_tmp_G.data_d;
			data_B=data_tmp_B.data_d;
			data_hB=data_tmp_hB.data_d;
			data_BIJ=data_tmp_BIJ.data_d;
			data_M_P_TLPF=data_tmp_M_P_TLPF.data_d;
			data_M_Q_TLPF=data_tmp_M_Q_TLPF.data_d;

			data_Ri_eP=data_tmp_Ri_eP.data_d;
			data_Ri_eQ=data_tmp_Ri_eQ.data_d;
			data_M_P_TLPF_reverse=data_tmp_M_P_TLPF_reverse.data_d;
			data_M_Q_TLPF_reverse=data_tmp_M_Q_TLPF_reverse.data_d;
			data_Ri_eP_reverse=data_tmp_Ri_eP_reverse.data_d;
			data_Ri_eQ_reverse=data_tmp_Ri_eQ_reverse.data_d;


			l_edge_const[i+j].K=data_K;
			l_edge_const[i+j].G=data_G;
			l_edge_const[i+j].B=data_B;
			l_edge_const[i+j].hB=data_hB;
			l_edge_const[i+j].BIJ=data_BIJ;
			l_edge_const[i+j].M_P_TLPF=data_M_P_TLPF;
			l_edge_const[i+j].M_Q_TLPF=data_M_Q_TLPF;

			l_edge_const[i+j].Ri_eP=data_Ri_eP;
			l_edge_const[i+j].Ri_eQ=data_Ri_eQ;
			l_edge_const[i+j].M_P_TLPF_reverse=data_M_P_TLPF_reverse;
			l_edge_const[i+j].M_Q_TLPF_reverse=data_M_Q_TLPF_reverse;
			l_edge_const[i+j].Ri_eP_reverse=data_Ri_eP_reverse;
			l_edge_const[i+j].Ri_eQ_reverse=data_Ri_eQ_reverse;

		}
	}
}

void V_bank_memcpy(const TYPE* const  node_V,double*  l_Vm_dst,double  l_Vm_src[4][2][MAX_NODE_NUM/4],
				double*  l_Va_dst,double  l_Va_src[4][2][MAX_NODE_NUM/4],
				const int begin,
				const int end
				){
	//#pragma HLS INLINE
#pragma HLS allocation instances=V_bank_memcpy limit=4 function
	const int node_num=end-begin;

	for(int i=0;i<node_num;i+=4){
		#pragma HLS PIPELINE II=1

		TYPE tmp=node_V[i/4];
		for(int j=0;j<4;j++){
			#pragma HLS unroll factor=4

			long data_l,data_l2;
			data_l=tmp.range(128*j+63,128*j).to_long();
			data_l2=tmp.range(128*j+127,128*j+64).to_long();

			union Int64_or_Double data_tmp,data_tmp2;
			data_tmp.data_i=data_l;
			data_tmp2.data_i=data_l2;

			double data,data2;
			data=data_tmp.data_d;
			data2=data_tmp2.data_d;
			l_Vm_dst[i+j]=data;
			l_Va_dst[i+j]=data2;
			for(int k=0;k<4;k++){
				#pragma HLS unroll factor=4
				for(int kk=0;kk<2;kk++){
					#pragma HLS unroll factor=2
					l_Vm_src[k][kk][i+j]=data;
					l_Va_src[k][kk][i+j]=data2;
				}
			}
		}
	}
}

/*
void Vm_bank_memcpy(const TYPE* const  node_Vm,double*  l_Vm_dst,double  l_Vm_src[4][2][MAX_NODE_NUM/4],
				const int begin,
				const int end
				){
	//#pragma HLS INLINE
#pragma HLS allocation instances=Vm_bank_memcpy limit=4 function
	const int node_num=end-begin;

	for(int i=0;i<node_num;i+=8){
		#pragma HLS PIPELINE II=1

		TYPE tmp=node_Vm[i/8];
		for(int j=0;j<8;j++){
			#pragma HLS unroll factor=8

			long data_l;
			data_l=tmp.range(64*j+63,64*j).to_long();

			union Int64_or_Double data_tmp;
			data_tmp.data_i=data_l;

			double data;
			data=data_tmp.data_d;
			l_Vm_dst[i+j]=data;
			for(int k=0;k<4;k++){
				#pragma HLS unroll factor=4
				for(int kk=0;kk<2;kk++){
					#pragma HLS unroll factor=2
					l_Vm_src[k][kk][i+j]=data;
				}
			}
		}
	}
}
void Va_bank_memcpy(const TYPE* const  node_Va,double*  l_Va_dst,double  l_Va_src[4][2][MAX_NODE_NUM/4],
				const int begin,
				const int end
				){
	//#pragma HLS INLINE
#pragma HLS allocation instances=Va_bank_memcpy limit=4 function
	const int node_num=end-begin;

	for(int i=0;i<node_num;i+=8){
		#pragma HLS PIPELINE II=1

		TYPE tmp=node_Va[i/8];
		for(int j=0;j<8;j++){
			#pragma HLS unroll factor=8

			long data_l;
			data_l=tmp.range(64*j+63,64*j).to_long();

			union Int64_or_Double data_tmp;
			data_tmp.data_i=data_l;

			double data;
			data=data_tmp.data_d;
			l_Va_dst[i+j]=data;
			for(int k=0;k<4;k++){
				#pragma HLS unroll factor=4
				for(int kk=0;kk<2;kk++){
					#pragma HLS unroll factor=2
					l_Va_src[k][kk][i+j]=data;
				}
			}
		}
	}
}
*/
/*void H_r_PQ_bank_memcpy(TYPE*  node_H_r_PQ,const double* const  l_H_r_P,const double* const  l_H_r_Q,
				const int begin,
				const int end
				){
	//#pragma HLS INLINE
#pragma HLS allocation instances=H_r_PQ_bank_memcpy limit=4 function
	const int node_num=end-begin;

	for(int i=0;i<node_num;i+=4){
		#pragma HLS PIPELINE II=1

		TYPE tmp;
		for(int j=0;j<4;j++){
			#pragma HLS unroll factor=4

			double data=l_H_r_P[i+j];
			double data2=l_H_r_Q[i+j];
			union Int64_or_Double data_tmp,data_tmp2;
			data_tmp.data_d=data;
			data_tmp2.data_d=data2;
			long data_l,data_l2;
			data_l=data_tmp.data_i;
			data_l2=data_tmp2.data_i;
			tmp.range(128*j+63,128*j)=data_l;
			tmp.range(128*j+127,128*j+64)=data_l2;
		}
		node_H_r_PQ[i/4]=tmp;
	}
}*/
/*
void H_r_P_bank_memcpy(TYPE*  node_H_r_P,const double* const  l_H_r_P,
				const int begin,
				const int end
				){
	//#pragma HLS INLINE
#pragma HLS allocation instances=H_r_P_bank_memcpy limit=4 function
	const int node_num=end-begin;

	for(int i=0;i<node_num;i+=8){
		#pragma HLS PIPELINE II=1

		TYPE tmp;
		for(int j=0;j<8;j++){
			#pragma HLS unroll factor=8

			double data=l_H_r_P[i+j];
			union Int64_or_Double data_tmp;
			data_tmp.data_d=data;
			long data_l;
			data_l=data_tmp.data_i;
			tmp.range(64*j+63,64*j)=data_l;
		}
		node_H_r_P[i/8]=tmp;
	}
}
void H_r_Q_bank_memcpy(TYPE*  node_H_r_Q,const double* const  l_H_r_Q,
				const int begin,
				const int end
				){
	//#pragma HLS INLINE
#pragma HLS allocation instances=H_r_Q_bank_memcpy limit=4 function
	const int node_num=end-begin;

	for(int i=0;i<node_num;i+=8){
		#pragma HLS PIPELINE II=1

		TYPE tmp;
		for(int j=0;j<8;j++){
			#pragma HLS unroll factor=8

			double data=l_H_r_Q[i+j];
			union Int64_or_Double data_tmp;
			data_tmp.data_d=data;
			long data_l;
			data_l=data_tmp.data_i;
			tmp.range(64*j+63,64*j)=data_l;
		}
		node_H_r_Q[i/8]=tmp;
	}
}
*/

#define F2L_PRECISION (1024*1024*256)

void krnl_se_pipe(
		const double* const  l_Vm_dst,
		const double* const  l_Va_dst,
		const double  l_Vm_src_0[2][MAX_NODE_NUM/4],const double  l_Vm_src_1[2][MAX_NODE_NUM/4],const double  l_Vm_src_2[2][MAX_NODE_NUM/4],const double  l_Vm_src_3[2][MAX_NODE_NUM/4],
		const double  l_Va_src_0[2][MAX_NODE_NUM/4],const double  l_Va_src_1[2][MAX_NODE_NUM/4],const double  l_Va_src_2[2][MAX_NODE_NUM/4],const double  l_Va_src_3[2][MAX_NODE_NUM/4],
		const int l_src_offset_0,const int l_src_offset_1,const int l_src_offset_2,const int l_src_offset_3,
		const struct Node_SE_CONST* const  l_node_const,
		const int* const  l_off_0,const int* const  l_off_1,
		const struct Edge_SE_CONST* const  l_edge_const,
		/*double*  l_H_r_P,*/long*  ll_H_r_P,
		/*double*  l_H_r_Q,*/long*  ll_H_r_Q,
		TYPE*  node_H_r_PQ,
		const int begin,
		const int end
		){
//#pragma HLS INLINE
#pragma HLS allocation instances=krnl_se_pipe limit=4 function

	const int l_node_num=end-begin;

	for(int i=0;i<l_node_num;i+=8){
		#pragma HLS PIPELINE II=1
		for(int j=0;j<8;j++){
			#pragma HLS UNROLL factor=8
			ll_H_r_P[i+j]=0;
			ll_H_r_Q[i+j]=0;
		}
	}

	const int edge_offset=l_off_0[0];

	const int edge_num=l_off_0[l_node_num]-edge_offset;
	int edge_count=0;
	int dst_id=0;
	while(edge_count<edge_num&&dst_id<l_node_num){
		#pragma HLS PIPELINE II=1

		int buf_off_0,buf_off_1,buf_off_2,buf_off_3;
		int buf_off_4,buf_off_5,buf_off_6,buf_off_7;
		int buf_off_8;

		if(dst_id+8<l_node_num){
			buf_off_8=l_off_1[dst_id+8]-edge_offset;
		}
		else{
			buf_off_8=edge_num;
		}


		const int old_dst_id=dst_id;
		const int old_edge_count=edge_count;

		if(buf_off_8>edge_count+4){
			dst_id+=0;edge_count+=4;
		}
		else if(buf_off_8==edge_count+4){
			dst_id+=8;edge_count+=4;
		}
		else{
			dst_id+=8;
		}

		for(int i=0;i<8;i++){
			#pragma HLS UNROLL factor=8
			int tmp_off;
			if(old_dst_id+i<l_node_num){
				tmp_off=l_off_0[old_dst_id/8*8+i]-edge_offset;
			}
			else{
				tmp_off=edge_num;
			}
			switch(i){
				case 0:
					buf_off_0=tmp_off;break;
				case 1:
					buf_off_1=tmp_off;break;
				case 2:
					buf_off_2=tmp_off;break;
				case 3:
					buf_off_3=tmp_off;break;
				case 4:
					buf_off_4=tmp_off;break;
				case 5:
					buf_off_5=tmp_off;break;
				case 6:
					buf_off_6=tmp_off;break;
				case 7:
					buf_off_7=tmp_off;break;
			}
		}


		double dst_Vm[8];
		double dst_Va[8];
		struct Node_SE_CONST dst_node_const[8];
		for(int j=0;j<8;j++){
			#pragma HLS UNROLL factor=8
			dst_Vm[j]=l_Vm_dst[old_dst_id/8*8+j];
			dst_Va[j]=l_Va_dst[old_dst_id/8*8+j];
			dst_node_const[j]=l_node_const[old_dst_id/8*8+j];
		}

		long buf_ll_H_r_P_1_0,buf_ll_H_r_P_1_1,buf_ll_H_r_P_1_2,buf_ll_H_r_P_1_3;
		//long buf_ll_H_r_P_1_4,buf_ll_H_r_P_1_5,buf_ll_H_r_P_1_6,buf_ll_H_r_P_1_7;
		//long buf_ll_H_r_P_1_8,buf_ll_H_r_P_1_9,buf_ll_H_r_P_1_10,buf_ll_H_r_P_1_11;
		//long buf_ll_H_r_P_1_12,buf_ll_H_r_P_1_13,buf_ll_H_r_P_1_14,buf_ll_H_r_P_1_15;
		long buf_ll_H_r_Q_1_0,buf_ll_H_r_Q_1_1,buf_ll_H_r_Q_1_2,buf_ll_H_r_Q_1_3;
		//long buf_ll_H_r_Q_1_4,buf_ll_H_r_Q_1_5,buf_ll_H_r_Q_1_6,buf_ll_H_r_Q_1_7;
		//long buf_ll_H_r_Q_1_8,buf_ll_H_r_Q_1_9,buf_ll_H_r_Q_1_10,buf_ll_H_r_Q_1_11;
		//long buf_ll_H_r_Q_1_12,buf_ll_H_r_Q_1_13,buf_ll_H_r_Q_1_14,buf_ll_H_r_Q_1_15;

		/*long buf_ll_H_r_P_2_0,buf_ll_H_r_P_2_1,buf_ll_H_r_P_2_2,buf_ll_H_r_P_2_3;
		long buf_ll_H_r_P_2_4,buf_ll_H_r_P_2_5,buf_ll_H_r_P_2_6,buf_ll_H_r_P_2_7;
		long buf_ll_H_r_P_2_8,buf_ll_H_r_P_2_9,buf_ll_H_r_P_2_10,buf_ll_H_r_P_2_11;
		long buf_ll_H_r_P_2_12,buf_ll_H_r_P_2_13,buf_ll_H_r_P_2_14,buf_ll_H_r_P_2_15;
		long buf_ll_H_r_Q_2_0,buf_ll_H_r_Q_2_1,buf_ll_H_r_Q_2_2,buf_ll_H_r_Q_2_3;
		long buf_ll_H_r_Q_2_4,buf_ll_H_r_Q_2_5,buf_ll_H_r_Q_2_6,buf_ll_H_r_Q_2_7;
		long buf_ll_H_r_Q_2_8,buf_ll_H_r_Q_2_9,buf_ll_H_r_Q_2_10,buf_ll_H_r_Q_2_11;
		long buf_ll_H_r_Q_2_12,buf_ll_H_r_Q_2_13,buf_ll_H_r_Q_2_14,buf_ll_H_r_Q_2_15;

		long buf_ll_H_r_P_3_0,buf_ll_H_r_P_3_1,buf_ll_H_r_P_3_2,buf_ll_H_r_P_3_3;
		long buf_ll_H_r_P_3_4,buf_ll_H_r_P_3_5,buf_ll_H_r_P_3_6,buf_ll_H_r_P_3_7;
		long buf_ll_H_r_P_3_8,buf_ll_H_r_P_3_9,buf_ll_H_r_P_3_10,buf_ll_H_r_P_3_11;
		long buf_ll_H_r_P_3_12,buf_ll_H_r_P_3_13,buf_ll_H_r_P_3_14,buf_ll_H_r_P_3_15;
		long buf_ll_H_r_Q_3_0,buf_ll_H_r_Q_3_1,buf_ll_H_r_Q_3_2,buf_ll_H_r_Q_3_3;
		long buf_ll_H_r_Q_3_4,buf_ll_H_r_Q_3_5,buf_ll_H_r_Q_3_6,buf_ll_H_r_Q_3_7;
		long buf_ll_H_r_Q_3_8,buf_ll_H_r_Q_3_9,buf_ll_H_r_Q_3_10,buf_ll_H_r_Q_3_11;
		long buf_ll_H_r_Q_3_12,buf_ll_H_r_Q_3_13,buf_ll_H_r_Q_3_14,buf_ll_H_r_Q_3_15;

		long buf_ll_H_r_P_4_0,buf_ll_H_r_P_4_1,buf_ll_H_r_P_4_2,buf_ll_H_r_P_4_3;
		long buf_ll_H_r_P_4_4,buf_ll_H_r_P_4_5,buf_ll_H_r_P_4_6,buf_ll_H_r_P_4_7;
		long buf_ll_H_r_P_4_8,buf_ll_H_r_P_4_9,buf_ll_H_r_P_4_10,buf_ll_H_r_P_4_11;
		long buf_ll_H_r_P_4_12,buf_ll_H_r_P_4_13,buf_ll_H_r_P_4_14,buf_ll_H_r_P_4_15;
		long buf_ll_H_r_Q_4_0,buf_ll_H_r_Q_4_1,buf_ll_H_r_Q_4_2,buf_ll_H_r_Q_4_3;
		long buf_ll_H_r_Q_4_4,buf_ll_H_r_Q_4_5,buf_ll_H_r_Q_4_6,buf_ll_H_r_Q_4_7;
		long buf_ll_H_r_Q_4_8,buf_ll_H_r_Q_4_9,buf_ll_H_r_Q_4_10,buf_ll_H_r_Q_4_11;
		long buf_ll_H_r_Q_4_12,buf_ll_H_r_Q_4_13,buf_ll_H_r_Q_4_14,buf_ll_H_r_Q_4_15;*/

		long buf_ll_H_r_P_5_0,buf_ll_H_r_P_5_1,buf_ll_H_r_P_5_2,buf_ll_H_r_P_5_3;
		//long buf_ll_H_r_P_5_4,buf_ll_H_r_P_5_5,buf_ll_H_r_P_5_6,buf_ll_H_r_P_5_7;
		//long buf_ll_H_r_P_5_8,buf_ll_H_r_P_5_9,buf_ll_H_r_P_5_10,buf_ll_H_r_P_5_11;
		//long buf_ll_H_r_P_5_12,buf_ll_H_r_P_5_13,buf_ll_H_r_P_5_14,buf_ll_H_r_P_5_15;
		long buf_ll_H_r_Q_5_0,buf_ll_H_r_Q_5_1,buf_ll_H_r_Q_5_2,buf_ll_H_r_Q_5_3;
		//long buf_ll_H_r_Q_5_4,buf_ll_H_r_Q_5_5,buf_ll_H_r_Q_5_6,buf_ll_H_r_Q_5_7;
		//long buf_ll_H_r_Q_5_8,buf_ll_H_r_Q_5_9,buf_ll_H_r_Q_5_10,buf_ll_H_r_Q_5_11;
		//long buf_ll_H_r_Q_5_12,buf_ll_H_r_Q_5_13,buf_ll_H_r_Q_5_14,buf_ll_H_r_Q_5_15;

		for(int j=0;j<4;j++){
			#pragma HLS UNROLL factor=4

			struct Edge_SE_CONST e=l_edge_const[old_edge_count/4*4+j];
			double src_Vm=0;
			double src_Va=0;


			int src_id=e.src;
			if(src_id>=l_src_offset_3){
				src_Vm=l_Vm_src_3[j%2][src_id-l_src_offset_3];
				src_Va=l_Va_src_3[j%2][src_id-l_src_offset_3];
			}
			else if(src_id>=l_src_offset_2){
				src_Vm=l_Vm_src_2[j%2][src_id-l_src_offset_2];
				src_Va=l_Va_src_2[j%2][src_id-l_src_offset_2];
			}
			else if(src_id>=l_src_offset_1){
				src_Vm=l_Vm_src_1[j%2][src_id-l_src_offset_1];
				src_Va=l_Va_src_1[j%2][src_id-l_src_offset_1];
			}
			else{
				src_Vm=l_Vm_src_0[j%2][src_id-l_src_offset_0];
				src_Va=l_Va_src_0[j%2][src_id-l_src_offset_0];
			}
			double delta_Va;
			double edge_dst_Vm,edge_dst_Va;
			if(old_edge_count+j<buf_off_1) {edge_dst_Va=dst_Va[0];edge_dst_Vm=dst_Vm[0];}
			else if(old_edge_count+j<buf_off_2) {edge_dst_Va=dst_Va[1];edge_dst_Vm=dst_Vm[1];}
			else if(old_edge_count+j<buf_off_3) {edge_dst_Va=dst_Va[2];edge_dst_Vm=dst_Vm[2];}
			else if(old_edge_count+j<buf_off_4) {edge_dst_Va=dst_Va[3];edge_dst_Vm=dst_Vm[3];}
			else if(old_edge_count+j<buf_off_5) {edge_dst_Va=dst_Va[4];edge_dst_Vm=dst_Vm[4];}
			else if(old_edge_count+j<buf_off_6) {edge_dst_Va=dst_Va[5];edge_dst_Vm=dst_Vm[5];}
			else if(old_edge_count+j<buf_off_7) {edge_dst_Va=dst_Va[6];edge_dst_Vm=dst_Vm[6];}
			else {edge_dst_Va=dst_Va[7];edge_dst_Vm=dst_Vm[7];}
			delta_Va=edge_dst_Va-src_Va;
			double cosine_array=cos(delta_Va);//hls::cordic::
			double sine_array=sin(delta_Va);

			double d_H_r_P;//=- src_Vm*(src_G*cosine_array +src_B*sine_array);
			double d_H_r_Q;// = - src_Vm*(src_G*sine_array - src_B*cosine_array);


			double tap_ratio = fabs(e.K/e.Kcount);
			double tap_ratio_square = (e.K/e.Kcount)*(e.K/e.Kcount);

			if(e.K == 0){
				d_H_r_P = e.BIJ * (e.M_P_TLPF - (edge_dst_Vm * edge_dst_Vm * e.G - edge_dst_Vm * src_Vm * (e.G*cosine_array + (-e.B)*(sine_array)))) * e.Ri_eP
						+(-1) * e.BIJ * (e.M_P_TLPF_reverse - (src_Vm * src_Vm * e.G - src_Vm * edge_dst_Vm * (e.G*cosine_array + (-e.B)*(-sine_array)))) * e.Ri_eP_reverse;
				d_H_r_Q = (e.B - e.hB) * (e.M_Q_TLPF - (- edge_dst_Vm * edge_dst_Vm * (-e.B + 0.5*e.hB) - edge_dst_Vm * src_Vm * (e.G*sine_array - (-e.B)*cosine_array))) * e.Ri_eQ
			    		+(-1) * e.B * (e.M_Q_TLPF_reverse - (- src_Vm * src_Vm * (-e.B + 0.5*e.hB) - src_Vm * edge_dst_Vm * (e.G*(-sine_array) - (-e.B)*cosine_array))) * e.Ri_eQ_reverse;
			//if(i==0&&iter==0) printf("(%d<-(%d,%lf))\n",i,p,e.BIJ * (e.M_P_TLPF - (node[i].Vm * node[i].Vm * e.G - node[i].Vm * node[p].Vm * (e.G*cos(node[i].Va-node[p].Va) + (-e.B)*sin(node[i].Va-node[p].Va)))) * e.Ri_eP);
			}
			else if(e.K > 0){
				d_H_r_P = e.BIJ * (e.M_P_TLPF - (edge_dst_Vm * edge_dst_Vm * (e.G/tap_ratio_square) - edge_dst_Vm * src_Vm * ((e.G/tap_ratio)*cosine_array + (-e.B/tap_ratio)*(sine_array)))) * e.Ri_eP
			    		+(-1) * e.BIJ * (e.M_P_TLPF_reverse - (src_Vm * src_Vm * e.G - src_Vm * edge_dst_Vm * ((e.G/tap_ratio)*cosine_array + (-e.B/tap_ratio)*(-sine_array)))) * e.Ri_eP_reverse;
				d_H_r_Q = (e.B/tap_ratio - e.hB) * (e.M_Q_TLPF - (- edge_dst_Vm * edge_dst_Vm * (-e.B + 0.5*e.hB) / tap_ratio_square - edge_dst_Vm * src_Vm * ((e.G/tap_ratio)*(sine_array) - (-e.B/tap_ratio)*cosine_array))) * e.Ri_eQ
			    		+(-1) * (e.B/tap_ratio) * (e.M_Q_TLPF_reverse - (- src_Vm * src_Vm * (-e.B + 0.5*e.hB) - src_Vm * edge_dst_Vm * ((e.G/tap_ratio)*(-sine_array) - (-e.B/tap_ratio)*cosine_array))) * e.Ri_eQ_reverse;
			    //if(i==0&&iter==0) printf("(%d<-(%d,%lf))\n",i,p,e.BIJ * (e.M_P_TLPF - (node[i].Vm * node[i].Vm * (e.G/tap_ratio_square) - node[i].Vm * node[p].Vm * ((e.G/tap_ratio)*cos(node[i].Va-node[p].Va) + (-e.B/tap_ratio)*sin(node[i].Va-node[p].Va)))) * e.Ri_eP);
			}
			else{
				d_H_r_P = e.BIJ * (e.M_P_TLPF - (edge_dst_Vm * edge_dst_Vm * e.G - edge_dst_Vm * src_Vm * ((e.G/tap_ratio)*cosine_array + (-e.B/tap_ratio)*(sine_array)))) * e.Ri_eP
			    		+(-1) * e.BIJ * (e.M_P_TLPF_reverse - (src_Vm * src_Vm * (e.G/tap_ratio_square) - src_Vm * edge_dst_Vm * ((e.G/tap_ratio)*cosine_array + (-e.B/tap_ratio)*(-sine_array)))) * e.Ri_eP_reverse;
				d_H_r_Q = (e.B/tap_ratio - e.hB) * (e.M_Q_TLPF - (- edge_dst_Vm * edge_dst_Vm * (-e.B + 0.5*e.hB) - edge_dst_Vm * src_Vm * ((e.G/tap_ratio)*(sine_array) - (-e.B/tap_ratio)*cosine_array))) * e.Ri_eQ
			    		+(-1) * (e.B/tap_ratio) * (e.M_Q_TLPF_reverse - (- src_Vm * src_Vm * (-e.B + 0.5*e.hB) / tap_ratio_square - src_Vm * edge_dst_Vm * ((e.G/tap_ratio)*(-sine_array) - (-e.B/tap_ratio)*cosine_array))) * e.Ri_eQ_reverse;
			    //if(i==0&&iter==0) printf("(%d<-(%d,%lf))\n",i,p,e.BIJ * (e.M_P_TLPF - (node[i].Vm * node[i].Vm * e.G - node[i].Vm * node[p].Vm * ((e.G/tap_ratio)*cos(node[i].Va-node[p].Va) + (-e.B/tap_ratio)*sin(node[i].Va-node[p].Va)))) * e.Ri_eP);
			}


			//buf_ll_H_r_P_1[j]=(long)(d_H_r_P*F2L_PRECISION);
			//buf_ll_H_r_Q_1[j]=(long)(d_H_r_Q*F2L_PRECISION);
			long tmp_ll_H_r_P_1=(long)(d_H_r_P*F2L_PRECISION);
			long tmp_ll_H_r_Q_1=(long)(d_H_r_Q*F2L_PRECISION);
			switch(j){
			case 0:
				buf_ll_H_r_P_1_0=tmp_ll_H_r_P_1;buf_ll_H_r_Q_1_0=tmp_ll_H_r_Q_1;break;
			case 1:
				buf_ll_H_r_P_1_1=tmp_ll_H_r_P_1;buf_ll_H_r_Q_1_1=tmp_ll_H_r_Q_1;break;
			case 2:
				buf_ll_H_r_P_1_2=tmp_ll_H_r_P_1;buf_ll_H_r_Q_1_2=tmp_ll_H_r_Q_1;break;
			case 3:
				buf_ll_H_r_P_1_3=tmp_ll_H_r_P_1;buf_ll_H_r_Q_1_3=tmp_ll_H_r_Q_1;break;
			/*case 4:
				buf_ll_H_r_P_1_4=tmp_ll_H_r_P_1;buf_ll_H_r_Q_1_4=tmp_ll_H_r_Q_1;break;
			case 5:
				buf_ll_H_r_P_1_5=tmp_ll_H_r_P_1;buf_ll_H_r_Q_1_5=tmp_ll_H_r_Q_1;break;
			case 6:
				buf_ll_H_r_P_1_6=tmp_ll_H_r_P_1;buf_ll_H_r_Q_1_6=tmp_ll_H_r_Q_1;break;
			case 7:
				buf_ll_H_r_P_1_7=tmp_ll_H_r_P_1;buf_ll_H_r_Q_1_7=tmp_ll_H_r_Q_1;break;
			case 8:
				buf_ll_H_r_P_1_8=tmp_ll_H_r_P_1;buf_ll_H_r_Q_1_8=tmp_ll_H_r_Q_1;break;
			case 9:
				buf_ll_H_r_P_1_9=tmp_ll_H_r_P_1;buf_ll_H_r_Q_1_9=tmp_ll_H_r_Q_1;break;
			case 10:
				buf_ll_H_r_P_1_10=tmp_ll_H_r_P_1;buf_ll_H_r_Q_1_10=tmp_ll_H_r_Q_1;break;
			case 11:
				buf_ll_H_r_P_1_11=tmp_ll_H_r_P_1;buf_ll_H_r_Q_1_11=tmp_ll_H_r_Q_1;break;
			case 12:
				buf_ll_H_r_P_1_12=tmp_ll_H_r_P_1;buf_ll_H_r_Q_1_12=tmp_ll_H_r_Q_1;break;
			case 13:
				buf_ll_H_r_P_1_13=tmp_ll_H_r_P_1;buf_ll_H_r_Q_1_13=tmp_ll_H_r_Q_1;break;
			case 14:
				buf_ll_H_r_P_1_14=tmp_ll_H_r_P_1;buf_ll_H_r_Q_1_14=tmp_ll_H_r_Q_1;break;
			case 15:
				buf_ll_H_r_P_1_15=tmp_ll_H_r_P_1;buf_ll_H_r_Q_1_15=tmp_ll_H_r_Q_1;break;*/
			}
		}
/*
//level 1
		buf_ll_H_r_P_2_0=buf_ll_H_r_P_1_0;
		buf_ll_H_r_Q_2_0=buf_ll_H_r_Q_1_0;
		buf_ll_H_r_P_2_2=buf_ll_H_r_P_1_2;
		buf_ll_H_r_Q_2_2=buf_ll_H_r_Q_1_2;
		buf_ll_H_r_P_2_4=buf_ll_H_r_P_1_4;
		buf_ll_H_r_Q_2_4=buf_ll_H_r_Q_1_4;
		buf_ll_H_r_P_2_6=buf_ll_H_r_P_1_6;
		buf_ll_H_r_Q_2_6=buf_ll_H_r_Q_1_6;

		buf_ll_H_r_P_2_1=buf_ll_H_r_P_1_1+buf_ll_H_r_P_1_0;
		buf_ll_H_r_Q_2_1=buf_ll_H_r_Q_1_1+buf_ll_H_r_Q_1_0;
		buf_ll_H_r_P_2_3=buf_ll_H_r_P_1_3+buf_ll_H_r_P_1_2;
		buf_ll_H_r_Q_2_3=buf_ll_H_r_Q_1_3+buf_ll_H_r_Q_1_2;
		buf_ll_H_r_P_2_5=buf_ll_H_r_P_1_5+buf_ll_H_r_P_1_4;
		buf_ll_H_r_Q_2_5=buf_ll_H_r_Q_1_5+buf_ll_H_r_Q_1_4;
		buf_ll_H_r_P_2_7=buf_ll_H_r_P_1_7+buf_ll_H_r_P_1_6;
		buf_ll_H_r_Q_2_7=buf_ll_H_r_Q_1_7+buf_ll_H_r_Q_1_6;

		buf_ll_H_r_P_2_8=buf_ll_H_r_P_1_8;
		buf_ll_H_r_Q_2_8=buf_ll_H_r_Q_1_8;
		buf_ll_H_r_P_2_10=buf_ll_H_r_P_1_10;
		buf_ll_H_r_Q_2_10=buf_ll_H_r_Q_1_10;
		buf_ll_H_r_P_2_12=buf_ll_H_r_P_1_12;
		buf_ll_H_r_Q_2_12=buf_ll_H_r_Q_1_12;
		buf_ll_H_r_P_2_14=buf_ll_H_r_P_1_14;
		buf_ll_H_r_Q_2_14=buf_ll_H_r_Q_1_14;

		buf_ll_H_r_P_2_9=buf_ll_H_r_P_1_9+buf_ll_H_r_P_1_8;
		buf_ll_H_r_Q_2_9=buf_ll_H_r_Q_1_9+buf_ll_H_r_Q_1_8;
		buf_ll_H_r_P_2_11=buf_ll_H_r_P_1_11+buf_ll_H_r_P_1_10;
		buf_ll_H_r_Q_2_11=buf_ll_H_r_Q_1_11+buf_ll_H_r_Q_1_10;
		buf_ll_H_r_P_2_13=buf_ll_H_r_P_1_13+buf_ll_H_r_P_1_12;
		buf_ll_H_r_Q_2_13=buf_ll_H_r_Q_1_13+buf_ll_H_r_Q_1_12;
		buf_ll_H_r_P_2_15=buf_ll_H_r_P_1_15+buf_ll_H_r_P_1_14;
		buf_ll_H_r_Q_2_15=buf_ll_H_r_Q_1_15+buf_ll_H_r_Q_1_14;


//level 2
		buf_ll_H_r_P_3_0=buf_ll_H_r_P_2_0;
		buf_ll_H_r_Q_3_0=buf_ll_H_r_Q_2_0;
		buf_ll_H_r_P_3_1=buf_ll_H_r_P_2_1;
		buf_ll_H_r_Q_3_1=buf_ll_H_r_Q_2_1;
		buf_ll_H_r_P_3_4=buf_ll_H_r_P_2_4;
		buf_ll_H_r_Q_3_4=buf_ll_H_r_Q_2_4;
		buf_ll_H_r_P_3_5=buf_ll_H_r_P_2_5;
		buf_ll_H_r_Q_3_5=buf_ll_H_r_Q_2_5;

		buf_ll_H_r_P_3_2=buf_ll_H_r_P_2_2+buf_ll_H_r_P_2_1;
		buf_ll_H_r_Q_3_2=buf_ll_H_r_Q_2_2+buf_ll_H_r_Q_2_1;
		buf_ll_H_r_P_3_6=buf_ll_H_r_P_2_6+buf_ll_H_r_P_2_5;
		buf_ll_H_r_Q_3_6=buf_ll_H_r_Q_2_6+buf_ll_H_r_Q_2_5;

		buf_ll_H_r_P_3_3=buf_ll_H_r_P_2_3+buf_ll_H_r_P_2_1;
		buf_ll_H_r_Q_3_3=buf_ll_H_r_Q_2_3+buf_ll_H_r_Q_2_1;
		buf_ll_H_r_P_3_7=buf_ll_H_r_P_2_7+buf_ll_H_r_P_2_5;
		buf_ll_H_r_Q_3_7=buf_ll_H_r_Q_2_7+buf_ll_H_r_Q_2_5;

		buf_ll_H_r_P_3_8=buf_ll_H_r_P_2_8;
		buf_ll_H_r_Q_3_8=buf_ll_H_r_Q_2_8;
		buf_ll_H_r_P_3_9=buf_ll_H_r_P_2_9;
		buf_ll_H_r_Q_3_9=buf_ll_H_r_Q_2_9;
		buf_ll_H_r_P_3_12=buf_ll_H_r_P_2_12;
		buf_ll_H_r_Q_3_12=buf_ll_H_r_Q_2_12;
		buf_ll_H_r_P_3_13=buf_ll_H_r_P_2_13;
		buf_ll_H_r_Q_3_13=buf_ll_H_r_Q_2_13;

		buf_ll_H_r_P_3_10=buf_ll_H_r_P_2_10+buf_ll_H_r_P_2_9;
		buf_ll_H_r_Q_3_10=buf_ll_H_r_Q_2_10+buf_ll_H_r_Q_2_9;
		buf_ll_H_r_P_3_14=buf_ll_H_r_P_2_14+buf_ll_H_r_P_2_13;
		buf_ll_H_r_Q_3_14=buf_ll_H_r_Q_2_14+buf_ll_H_r_Q_2_13;

		buf_ll_H_r_P_3_11=buf_ll_H_r_P_2_11+buf_ll_H_r_P_2_9;
		buf_ll_H_r_Q_3_11=buf_ll_H_r_Q_2_11+buf_ll_H_r_Q_2_9;
		buf_ll_H_r_P_3_15=buf_ll_H_r_P_2_15+buf_ll_H_r_P_2_13;
		buf_ll_H_r_Q_3_15=buf_ll_H_r_Q_2_15+buf_ll_H_r_Q_2_13;


//level 3
		buf_ll_H_r_P_4_0=buf_ll_H_r_P_3_0;
		buf_ll_H_r_Q_4_0=buf_ll_H_r_Q_3_0;
		buf_ll_H_r_P_4_1=buf_ll_H_r_P_3_1;
		buf_ll_H_r_Q_4_1=buf_ll_H_r_Q_3_1;
		buf_ll_H_r_P_4_2=buf_ll_H_r_P_3_2;
		buf_ll_H_r_Q_4_2=buf_ll_H_r_Q_3_2;
		buf_ll_H_r_P_4_3=buf_ll_H_r_P_3_3;
		buf_ll_H_r_Q_4_3=buf_ll_H_r_Q_3_3;

		buf_ll_H_r_P_4_4=buf_ll_H_r_P_3_4+buf_ll_H_r_P_3_3;
		buf_ll_H_r_Q_4_4=buf_ll_H_r_Q_3_4+buf_ll_H_r_Q_3_3;
		buf_ll_H_r_P_4_5=buf_ll_H_r_P_3_5+buf_ll_H_r_P_3_3;
		buf_ll_H_r_Q_4_5=buf_ll_H_r_Q_3_5+buf_ll_H_r_Q_3_3;
		buf_ll_H_r_P_4_6=buf_ll_H_r_P_3_6+buf_ll_H_r_P_3_3;
		buf_ll_H_r_Q_4_6=buf_ll_H_r_Q_3_6+buf_ll_H_r_Q_3_3;
		buf_ll_H_r_P_4_7=buf_ll_H_r_P_3_7+buf_ll_H_r_P_3_3;
		buf_ll_H_r_Q_4_7=buf_ll_H_r_Q_3_7+buf_ll_H_r_Q_3_3;

		buf_ll_H_r_P_4_8=buf_ll_H_r_P_3_8;
		buf_ll_H_r_Q_4_8=buf_ll_H_r_Q_3_8;
		buf_ll_H_r_P_4_9=buf_ll_H_r_P_3_9;
		buf_ll_H_r_Q_4_9=buf_ll_H_r_Q_3_9;
		buf_ll_H_r_P_4_10=buf_ll_H_r_P_3_10;
		buf_ll_H_r_Q_4_10=buf_ll_H_r_Q_3_10;
		buf_ll_H_r_P_4_11=buf_ll_H_r_P_3_11;
		buf_ll_H_r_Q_4_11=buf_ll_H_r_Q_3_11;

		buf_ll_H_r_P_4_12=buf_ll_H_r_P_3_12+buf_ll_H_r_P_3_11;
		buf_ll_H_r_Q_4_12=buf_ll_H_r_Q_3_12+buf_ll_H_r_Q_3_11;
		buf_ll_H_r_P_4_13=buf_ll_H_r_P_3_13+buf_ll_H_r_P_3_11;
		buf_ll_H_r_Q_4_13=buf_ll_H_r_Q_3_13+buf_ll_H_r_Q_3_11;
		buf_ll_H_r_P_4_14=buf_ll_H_r_P_3_14+buf_ll_H_r_P_3_11;
		buf_ll_H_r_Q_4_14=buf_ll_H_r_Q_3_14+buf_ll_H_r_Q_3_11;
		buf_ll_H_r_P_4_15=buf_ll_H_r_P_3_15+buf_ll_H_r_P_3_11;
		buf_ll_H_r_Q_4_15=buf_ll_H_r_Q_3_15+buf_ll_H_r_Q_3_11;


//level 4
		buf_ll_H_r_P_5_0=buf_ll_H_r_P_4_0;
		buf_ll_H_r_Q_5_0=buf_ll_H_r_Q_4_0;
		buf_ll_H_r_P_5_1=buf_ll_H_r_P_4_1;
		buf_ll_H_r_Q_5_1=buf_ll_H_r_Q_4_1;
		buf_ll_H_r_P_5_2=buf_ll_H_r_P_4_2;
		buf_ll_H_r_Q_5_2=buf_ll_H_r_Q_4_2;
		buf_ll_H_r_P_5_3=buf_ll_H_r_P_4_3;
		buf_ll_H_r_Q_5_3=buf_ll_H_r_Q_4_3;
		buf_ll_H_r_P_5_4=buf_ll_H_r_P_4_4;
		buf_ll_H_r_Q_5_4=buf_ll_H_r_Q_4_4;
		buf_ll_H_r_P_5_5=buf_ll_H_r_P_4_5;
		buf_ll_H_r_Q_5_5=buf_ll_H_r_Q_4_5;
		buf_ll_H_r_P_5_6=buf_ll_H_r_P_4_6;
		buf_ll_H_r_Q_5_6=buf_ll_H_r_Q_4_6;
		buf_ll_H_r_P_5_7=buf_ll_H_r_P_4_7;
		buf_ll_H_r_Q_5_7=buf_ll_H_r_Q_4_7;

		buf_ll_H_r_P_5_8=buf_ll_H_r_P_4_8+buf_ll_H_r_P_4_7;
		buf_ll_H_r_Q_5_8=buf_ll_H_r_Q_4_8+buf_ll_H_r_Q_4_7;
		buf_ll_H_r_P_5_9=buf_ll_H_r_P_4_9+buf_ll_H_r_P_4_7;
		buf_ll_H_r_Q_5_9=buf_ll_H_r_Q_4_9+buf_ll_H_r_Q_4_7;
		buf_ll_H_r_P_5_10=buf_ll_H_r_P_4_10+buf_ll_H_r_P_4_7;
		buf_ll_H_r_Q_5_10=buf_ll_H_r_Q_4_10+buf_ll_H_r_Q_4_7;
		buf_ll_H_r_P_5_11=buf_ll_H_r_P_4_11+buf_ll_H_r_P_4_7;
		buf_ll_H_r_Q_5_11=buf_ll_H_r_Q_4_11+buf_ll_H_r_Q_4_7;
		buf_ll_H_r_P_5_12=buf_ll_H_r_P_4_12+buf_ll_H_r_P_4_7;
		buf_ll_H_r_Q_5_12=buf_ll_H_r_Q_4_12+buf_ll_H_r_Q_4_7;
		buf_ll_H_r_P_5_13=buf_ll_H_r_P_4_13+buf_ll_H_r_P_4_7;
		buf_ll_H_r_Q_5_13=buf_ll_H_r_Q_4_13+buf_ll_H_r_Q_4_7;
		buf_ll_H_r_P_5_14=buf_ll_H_r_P_4_14+buf_ll_H_r_P_4_7;
		buf_ll_H_r_Q_5_14=buf_ll_H_r_Q_4_14+buf_ll_H_r_Q_4_7;
		buf_ll_H_r_P_5_15=buf_ll_H_r_P_4_15+buf_ll_H_r_P_4_7;
		buf_ll_H_r_Q_5_15=buf_ll_H_r_Q_4_15+buf_ll_H_r_Q_4_7;
*/
		buf_ll_H_r_P_5_0=buf_ll_H_r_P_1_0;
		buf_ll_H_r_P_5_1=buf_ll_H_r_P_1_0+buf_ll_H_r_P_1_1;
		buf_ll_H_r_P_5_2=buf_ll_H_r_P_1_0+buf_ll_H_r_P_1_1+buf_ll_H_r_P_1_2;
		buf_ll_H_r_P_5_3=buf_ll_H_r_P_1_0+buf_ll_H_r_P_1_1+buf_ll_H_r_P_1_2+buf_ll_H_r_P_1_3;

		buf_ll_H_r_Q_5_0=buf_ll_H_r_Q_1_0;
		buf_ll_H_r_Q_5_1=buf_ll_H_r_Q_1_0+buf_ll_H_r_Q_1_1;
		buf_ll_H_r_Q_5_2=buf_ll_H_r_Q_1_0+buf_ll_H_r_Q_1_1+buf_ll_H_r_Q_1_2;
		buf_ll_H_r_Q_5_3=buf_ll_H_r_Q_1_0+buf_ll_H_r_Q_1_1+buf_ll_H_r_Q_1_2+buf_ll_H_r_Q_1_3;


		int locat_0,locat_1,locat_2,locat_3,locat_4,locat_5,locat_6,locat_7,locat_8;
		long H_r_P_value_0,H_r_P_value_1,H_r_P_value_2,H_r_P_value_3;
		long H_r_P_value_4,H_r_P_value_5,H_r_P_value_6,H_r_P_value_7;
		long H_r_P_value_8;
		long H_r_Q_value_0,H_r_Q_value_1,H_r_Q_value_2,H_r_Q_value_3;
		long H_r_Q_value_4,H_r_Q_value_5,H_r_Q_value_6,H_r_Q_value_7;
		long H_r_Q_value_8;


#define get_predix_sum(locat_0,H_r_P_value_0,H_r_Q_value_0)\
		switch(locat_0){\
			case 0:\
				H_r_P_value_0=buf_ll_H_r_P_5_0;H_r_Q_value_0=buf_ll_H_r_Q_5_0;break;\
			case 1:\
				H_r_P_value_0=buf_ll_H_r_P_5_1;H_r_Q_value_0=buf_ll_H_r_Q_5_1;break;\
			case 2:\
				H_r_P_value_0=buf_ll_H_r_P_5_2;H_r_Q_value_0=buf_ll_H_r_Q_5_2;break;\
			case 3:\
				H_r_P_value_0=buf_ll_H_r_P_5_3;H_r_Q_value_0=buf_ll_H_r_Q_5_3;break;\
			/*case 4:\
				H_r_P_value_0=buf_ll_H_r_P_5_4;H_r_Q_value_0=buf_ll_H_r_Q_5_4;break;\
			case 5:\
				H_r_P_value_0=buf_ll_H_r_P_5_5;H_r_Q_value_0=buf_ll_H_r_Q_5_5;break;\
			case 6:\
				H_r_P_value_0=buf_ll_H_r_P_5_6;H_r_Q_value_0=buf_ll_H_r_Q_5_6;break;\
			case 7:\
				H_r_P_value_0=buf_ll_H_r_P_5_7;H_r_Q_value_0=buf_ll_H_r_Q_5_7;break;\
			case 8:\
				H_r_P_value_0=buf_ll_H_r_P_5_8;H_r_Q_value_0=buf_ll_H_r_Q_5_8;break;\
			case 9:\
				H_r_P_value_0=buf_ll_H_r_P_5_9;H_r_Q_value_0=buf_ll_H_r_Q_5_9;break;\
			case 10:\
				H_r_P_value_0=buf_ll_H_r_P_5_10;H_r_Q_value_0=buf_ll_H_r_Q_5_10;break;\
			case 11:\
				H_r_P_value_0=buf_ll_H_r_P_5_11;H_r_Q_value_0=buf_ll_H_r_Q_5_11;break;\
			case 12:\
				H_r_P_value_0=buf_ll_H_r_P_5_12;H_r_Q_value_0=buf_ll_H_r_Q_5_12;break;\
			case 13:\
				H_r_P_value_0=buf_ll_H_r_P_5_13;H_r_Q_value_0=buf_ll_H_r_Q_5_13;break;\
			case 14:\
				H_r_P_value_0=buf_ll_H_r_P_5_14;H_r_Q_value_0=buf_ll_H_r_Q_5_14;break;\
			case 15:\
				H_r_P_value_0=buf_ll_H_r_P_5_15;H_r_Q_value_0=buf_ll_H_r_Q_5_15;break;*/\
		}

#define get_predix_sum_2(buf_off_0,locat_0,H_r_P_value_0,H_r_Q_value_0) \
		if(buf_off_0-old_edge_count-1<0){\
			H_r_P_value_0=0;H_r_Q_value_0=0;\
		}\
		else{\
			if(buf_off_0-old_edge_count-1>3){\
				locat_0=3;\
			}\
			else{\
				locat_0=buf_off_0-old_edge_count-1;\
			}\
			get_predix_sum(locat_0,H_r_P_value_0,H_r_Q_value_0)\
		}

		get_predix_sum_2(buf_off_0,locat_0,H_r_P_value_0,H_r_Q_value_0)
		get_predix_sum_2(buf_off_1,locat_1,H_r_P_value_1,H_r_Q_value_1)
		get_predix_sum_2(buf_off_2,locat_2,H_r_P_value_2,H_r_Q_value_2)
		get_predix_sum_2(buf_off_3,locat_3,H_r_P_value_3,H_r_Q_value_3)
		get_predix_sum_2(buf_off_4,locat_4,H_r_P_value_4,H_r_Q_value_4)
		get_predix_sum_2(buf_off_5,locat_5,H_r_P_value_5,H_r_Q_value_5)
		get_predix_sum_2(buf_off_6,locat_6,H_r_P_value_6,H_r_Q_value_6)
		get_predix_sum_2(buf_off_7,locat_7,H_r_P_value_7,H_r_Q_value_7)
		get_predix_sum_2(buf_off_8,locat_8,H_r_P_value_8,H_r_Q_value_8)

		for(int j=0;j<8;j++){
			#pragma HLS UNROLL factor=8
			long tmpP,tmpQ;
			switch(j){
			case 0:
				tmpP=H_r_P_value_1-H_r_P_value_0;tmpQ=H_r_Q_value_1-H_r_Q_value_0;break;
			case 1:
				tmpP=H_r_P_value_2-H_r_P_value_1;tmpQ=H_r_Q_value_2-H_r_Q_value_1;break;
			case 2:
				tmpP=H_r_P_value_3-H_r_P_value_2;tmpQ=H_r_Q_value_3-H_r_Q_value_2;break;
			case 3:
				tmpP=H_r_P_value_4-H_r_P_value_3;tmpQ=H_r_Q_value_4-H_r_Q_value_3;break;
			case 4:
				tmpP=H_r_P_value_5-H_r_P_value_4;tmpQ=H_r_Q_value_5-H_r_Q_value_4;break;
			case 5:
				tmpP=H_r_P_value_6-H_r_P_value_5;tmpQ=H_r_Q_value_6-H_r_Q_value_5;break;
			case 6:
				tmpP=H_r_P_value_7-H_r_P_value_6;tmpQ=H_r_Q_value_7-H_r_Q_value_6;break;
			case 7:
				tmpP=H_r_P_value_8-H_r_P_value_7;tmpQ=H_r_Q_value_8-H_r_Q_value_7;break;
			}

			ll_H_r_P[old_dst_id/8*8+j]+=tmpP;
			ll_H_r_Q[old_dst_id/8*8+j]+=tmpQ;
		}
	}
	for(int i=0;i<l_node_num;i+=4){
			#pragma HLS PIPELINE II=1

			TYPE tmp;
			for(int j=0;j<4;j++){
				#pragma HLS UNROLL factor=4

				double d_l_H_r_P=(double)ll_H_r_P[i+j]/(double)F2L_PRECISION;
				double d_l_H_r_Q=(double)ll_H_r_Q[i+j]/(double)F2L_PRECISION;

				struct Node_SE_CONST dst_const=l_node_const[i+j];
				double dst_Vm=l_Vm_dst[i+j];

				double data=d_l_H_r_P;
				double data2=d_l_H_r_Q+(dst_const.M_Vm-dst_Vm)*dst_const.Ri_V;

				union Int64_or_Double data_tmp,data_tmp2;
				data_tmp.data_d=data;
				data_tmp2.data_d=data2;
				long data_l,data_l2;
				data_l=data_tmp.data_i;
				data_l2=data_tmp2.data_i;
				tmp.range(128*j+63,128*j)=data_l;
				tmp.range(128*j+127,128*j+64)=data_l2;
			}
			node_H_r_PQ[i/4]=tmp;
	}
}


extern "C" {
void krnl_se(
				const TYPE* const  node_V_0,
				const TYPE* const  node_V_1,
				const TYPE* const  node_V_2,
				const TYPE* const  node_V_3,
				//const TYPE* const  node_Va_0,
				//const TYPE* const  node_Va_1,
				//const TYPE* const  node_Va_2,
				//const TYPE* const  node_Va_3,
				const TYPE* const  node_const_0,
				const TYPE* const  node_const_1,
				const TYPE* const  node_const_2,
				const TYPE* const  node_const_3,
				const TYPE* const  edge_off_0,
				const TYPE* const  edge_off_1,
				const TYPE* const  edge_off_2,
				const TYPE* const  edge_off_3,
                const TYPE* const  edge_const_0,
				const TYPE* const  edge_const_1,
				const TYPE* const  edge_const_2,
				const TYPE* const  edge_const_3,
				TYPE*  node_H_r_PQ_0,
				TYPE*  node_H_r_PQ_1,
				TYPE*  node_H_r_PQ_2,
				TYPE*  node_H_r_PQ_3,
				//TYPE*  node_H_r_Q_0,
				//TYPE*  node_H_r_Q_1,
				//TYPE*  node_H_r_Q_2,
				//TYPE*  node_H_r_Q_3,
				const int begin_0,
				const int begin_1,
				const int begin_2,
				const int begin_3,
				const int end_0,
				const int end_1,
				const int end_2,
				const int end_3,
				const int node_num,
				const bool init_flag
				//const bool Vm_flag,
				//const bool Va_flag,
				//const bool H_r_P_flag,
				//const bool H_r_Q_flag,
				//const bool compute_flag
				)
{
#pragma HLS INTERFACE m_axi port = node_V_0 offset = slave bundle = gmem0 //depth=512
#pragma HLS INTERFACE m_axi port = node_V_1 offset = slave bundle = gmem1 //depth=512
#pragma HLS INTERFACE m_axi port = node_V_2 offset = slave bundle = gmem2 //depth=512
#pragma HLS INTERFACE m_axi port = node_V_3 offset = slave bundle = gmem3 //depth=512

//#pragma HLS INTERFACE m_axi port = node_Va_0 offset = slave bundle = gmem0 //depth=512
//#pragma HLS INTERFACE m_axi port = node_Va_1 offset = slave bundle = gmem1 //depth=512
//#pragma HLS INTERFACE m_axi port = node_Va_2 offset = slave bundle = gmem2 //depth=512
//#pragma HLS INTERFACE m_axi port = node_Va_3 offset = slave bundle = gmem3 //depth=512

#pragma HLS INTERFACE m_axi port = node_const_0 offset = slave bundle = gmem0 //depth=512
#pragma HLS INTERFACE m_axi port = node_const_1 offset = slave bundle = gmem1 //depth=512
#pragma HLS INTERFACE m_axi port = node_const_2 offset = slave bundle = gmem2 //depth=512
#pragma HLS INTERFACE m_axi port = node_const_3 offset = slave bundle = gmem3 //depth=512

#pragma HLS INTERFACE m_axi port = edge_off_0 offset = slave bundle = gmem0 //depth=512
#pragma HLS INTERFACE m_axi port = edge_off_1 offset = slave bundle = gmem1 //depth=512
#pragma HLS INTERFACE m_axi port = edge_off_2 offset = slave bundle = gmem2 //depth=512
#pragma HLS INTERFACE m_axi port = edge_off_3 offset = slave bundle = gmem3 //depth=512

#pragma HLS INTERFACE m_axi port = edge_const_0 offset = slave bundle = gmem0 //depth=512
#pragma HLS INTERFACE m_axi port = edge_const_1 offset = slave bundle = gmem1 //depth=512
#pragma HLS INTERFACE m_axi port = edge_const_2 offset = slave bundle = gmem2 //depth=512
#pragma HLS INTERFACE m_axi port = edge_const_3 offset = slave bundle = gmem3 //depth=512

#pragma HLS INTERFACE m_axi port = node_H_r_PQ_0 offset = slave bundle = gmem0 //depth=512
#pragma HLS INTERFACE m_axi port = node_H_r_PQ_1 offset = slave bundle = gmem1 //depth=512
#pragma HLS INTERFACE m_axi port = node_H_r_PQ_2 offset = slave bundle = gmem2 //depth=512
#pragma HLS INTERFACE m_axi port = node_H_r_PQ_3 offset = slave bundle = gmem3 //depth=512

//#pragma HLS INTERFACE m_axi port = node_H_r_Q_0 offset = slave bundle = gmem0 //depth=512
//#pragma HLS INTERFACE m_axi port = node_H_r_Q_1 offset = slave bundle = gmem1 //depth=512
//#pragma HLS INTERFACE m_axi port = node_H_r_Q_2 offset = slave bundle = gmem2 //depth=512
//#pragma HLS INTERFACE m_axi port = node_H_r_Q_3 offset = slave bundle = gmem3 //depth=512


#pragma HLS INTERFACE s_axilite port=node_V_0 bundle=control
#pragma HLS INTERFACE s_axilite port=node_V_1 bundle=control
#pragma HLS INTERFACE s_axilite port=node_V_2 bundle=control
#pragma HLS INTERFACE s_axilite port=node_V_3 bundle=control

//#pragma HLS INTERFACE s_axilite port=node_Va_0 bundle=control
//#pragma HLS INTERFACE s_axilite port=node_Va_1 bundle=control
//#pragma HLS INTERFACE s_axilite port=node_Va_2 bundle=control
//#pragma HLS INTERFACE s_axilite port=node_Va_3 bundle=control

#pragma HLS INTERFACE s_axilite port=node_const_0 bundle=control
#pragma HLS INTERFACE s_axilite port=node_const_1 bundle=control
#pragma HLS INTERFACE s_axilite port=node_const_2 bundle=control
#pragma HLS INTERFACE s_axilite port=node_const_3 bundle=control

#pragma HLS INTERFACE s_axilite port=edge_off_0 bundle=control
#pragma HLS INTERFACE s_axilite port=edge_off_1 bundle=control
#pragma HLS INTERFACE s_axilite port=edge_off_2 bundle=control
#pragma HLS INTERFACE s_axilite port=edge_off_3 bundle=control

#pragma HLS INTERFACE s_axilite port=edge_const_0 bundle=control
#pragma HLS INTERFACE s_axilite port=edge_const_1 bundle=control
#pragma HLS INTERFACE s_axilite port=edge_const_2 bundle=control
#pragma HLS INTERFACE s_axilite port=edge_const_3 bundle=control

#pragma HLS INTERFACE s_axilite port=node_H_r_PQ_0 bundle=control
#pragma HLS INTERFACE s_axilite port=node_H_r_PQ_1 bundle=control
#pragma HLS INTERFACE s_axilite port=node_H_r_PQ_2 bundle=control
#pragma HLS INTERFACE s_axilite port=node_H_r_PQ_3 bundle=control

//#pragma HLS INTERFACE s_axilite port=node_H_r_Q_0 bundle=control
//#pragma HLS INTERFACE s_axilite port=node_H_r_Q_1 bundle=control
//#pragma HLS INTERFACE s_axilite port=node_H_r_Q_2 bundle=control
//#pragma HLS INTERFACE s_axilite port=node_H_r_Q_3 bundle=control

#pragma HLS INTERFACE s_axilite port=begin_0 bundle=control
#pragma HLS INTERFACE s_axilite port=begin_1 bundle=control
#pragma HLS INTERFACE s_axilite port=begin_2 bundle=control
#pragma HLS INTERFACE s_axilite port=begin_3 bundle=control

#pragma HLS INTERFACE s_axilite port=end_0 bundle=control
#pragma HLS INTERFACE s_axilite port=end_1 bundle=control
#pragma HLS INTERFACE s_axilite port=end_2 bundle=control
#pragma HLS INTERFACE s_axilite port=end_3 bundle=control

#pragma HLS INTERFACE s_axilite port=node_num bundle=control
#pragma HLS INTERFACE s_axilite port=init_flag bundle=control
//#pragma HLS INTERFACE s_axilite port=Vm_flag bundle=control
//#pragma HLS INTERFACE s_axilite port=Va_flag bundle=control
//#pragma HLS INTERFACE s_axilite port=H_r_P_flag bundle=control
//#pragma HLS INTERFACE s_axilite port=H_r_Q_flag bundle=control
//#pragma HLS INTERFACE s_axilite port=compute_flag bundle=control

#pragma HLS INTERFACE s_axilite port=return bundle=control



	static struct Node_SE_CONST l_node_const[4][MAX_NODE_NUM/4];
#pragma HLS ARRAY_PARTITION variable=l_node_const complete dim=1
#pragma HLS ARRAY_PARTITION variable=l_node_const dim=2 cyclic factor=4

	static int l_off[4][2][MAX_NODE_NUM/4];
#pragma HLS ARRAY_PARTITION variable=l_off complete dim=1
#pragma HLS ARRAY_PARTITION variable=l_off complete dim=2
#pragma HLS ARRAY_PARTITION variable=l_off dim=3 cyclic factor=4
	static struct Edge_SE_CONST l_edge_const[4][MAX_EDGE_NUM/4];
#pragma HLS ARRAY_PARTITION variable=l_edge_const complete dim=1
#pragma HLS ARRAY_PARTITION variable=l_edge_const dim=2 cyclic factor=4


	static double l_Vm_dst[4][MAX_NODE_NUM/4];
#pragma HLS ARRAY_PARTITION variable=l_Vm_dst complete dim=1
#pragma HLS ARRAY_PARTITION variable=l_Vm_dst dim=2 cyclic factor=4

	static double l_Va_dst[4][MAX_NODE_NUM/4];
#pragma HLS ARRAY_PARTITION variable=l_Va_dst complete dim=1
#pragma HLS ARRAY_PARTITION variable=l_Va_dst dim=2 cyclic factor=4

	static double l_Vm_src[4][4][2][MAX_NODE_NUM/4];
#pragma HLS ARRAY_PARTITION variable=l_Vm_src complete dim=1
#pragma HLS ARRAY_PARTITION variable=l_Vm_src complete dim=2
#pragma HLS ARRAY_PARTITION variable=l_Vm_src complete dim=3
#pragma HLS ARRAY_PARTITION variable=l_Vm_src dim=4 cyclic factor=4
	static double l_Va_src[4][4][2][MAX_NODE_NUM/4];
#pragma HLS ARRAY_PARTITION variable=l_Va_src complete dim=1
#pragma HLS ARRAY_PARTITION variable=l_Va_src complete dim=2
#pragma HLS ARRAY_PARTITION variable=l_Va_src complete dim=3
#pragma HLS ARRAY_PARTITION variable=l_Va_src dim=4 cyclic factor=4

	/*double  l_H_r_P[4][MAX_NODE_NUM/4];
#pragma HLS ARRAY_PARTITION variable=l_H_r_P complete dim=1
#pragma HLS ARRAY_PARTITION variable=l_H_r_P dim=2 cyclic factor=4
	double  l_H_r_Q[4][MAX_NODE_NUM/4];
#pragma HLS ARRAY_PARTITION variable=l_H_r_Q complete dim=1
#pragma HLS ARRAY_PARTITION variable=l_H_r_Q dim=2 cyclic factor=4
*/
	long  ll_H_r_P[4][MAX_NODE_NUM/4];
	#pragma HLS ARRAY_PARTITION variable=ll_H_r_P complete dim=1
	#pragma HLS ARRAY_PARTITION variable=ll_H_r_P dim=2 cyclic factor=4
	long  ll_H_r_Q[4][MAX_NODE_NUM/4];
	#pragma HLS ARRAY_PARTITION variable=ll_H_r_Q complete dim=1
	#pragma HLS ARRAY_PARTITION variable=ll_H_r_Q dim=2 cyclic factor=4


if(init_flag){
	init_bank_memcpy(node_const_0,l_node_const[0],
						edge_off_0,l_off[0][0],l_off[0][1],
		                edge_const_0,l_edge_const[0],
						begin_0,
						end_0
						);
		init_bank_memcpy(node_const_1,l_node_const[1],
							edge_off_1,l_off[1][0],l_off[1][1],
			                edge_const_1,l_edge_const[1],
							begin_1,
							end_1
							);
		init_bank_memcpy(node_const_2,l_node_const[2],
							edge_off_2,l_off[2][0],l_off[2][1],
			                edge_const_2,l_edge_const[2],
							begin_2,
							end_2
							);
		init_bank_memcpy(node_const_3,l_node_const[3],
							edge_off_3,l_off[3][0],l_off[3][1],
			                edge_const_3,l_edge_const[3],
							begin_3,
							end_3
							);
}

	V_bank_memcpy(node_V_0,l_Vm_dst[0],l_Vm_src[0],l_Va_dst[0],l_Va_src[0],begin_0,end_0);
	V_bank_memcpy(node_V_1,l_Vm_dst[1],l_Vm_src[1],l_Va_dst[1],l_Va_src[1],begin_1,end_1);
	V_bank_memcpy(node_V_2,l_Vm_dst[2],l_Vm_src[2],l_Va_dst[2],l_Va_src[2],begin_2,end_2);
	V_bank_memcpy(node_V_3,l_Vm_dst[3],l_Vm_src[3],l_Va_dst[3],l_Va_src[3],begin_3,end_3);

/*if(Vm_flag){
	Vm_bank_memcpy(node_Vm_0,l_Vm_dst[0],l_Vm_src[0],begin_0,end_0);
	Vm_bank_memcpy(node_Vm_1,l_Vm_dst[1],l_Vm_src[1],begin_1,end_1);
	Vm_bank_memcpy(node_Vm_2,l_Vm_dst[2],l_Vm_src[2],begin_2,end_2);
	Vm_bank_memcpy(node_Vm_3,l_Vm_dst[3],l_Vm_src[3],begin_3,end_3);
}
if(Va_flag){
	Va_bank_memcpy(node_Va_0,l_Va_dst[0],l_Va_src[0],begin_0,end_0);
	Va_bank_memcpy(node_Va_1,l_Va_dst[1],l_Va_src[1],begin_1,end_1);
	Va_bank_memcpy(node_Va_2,l_Va_dst[2],l_Va_src[2],begin_2,end_2);
	Va_bank_memcpy(node_Va_3,l_Va_dst[3],l_Va_src[3],begin_3,end_3);
}*/
	int begin[4],end[4];
	begin[0]=begin_0;end[0]=end_0;
	begin[1]=begin_1;end[1]=end_1;
	begin[2]=begin_2;end[2]=end_2;
	begin[3]=begin_3;end[3]=end_3;


	krnl_se_pipe(
					l_Vm_dst[0],
					l_Va_dst[0],
					l_Vm_src[0][0],l_Vm_src[1][0],l_Vm_src[2][0],l_Vm_src[3][0],
					l_Va_src[0][0],l_Va_src[1][0],l_Va_src[2][0],l_Va_src[3][0],
					begin[0],begin[1],begin[2],begin[3],
					l_node_const[0],
					l_off[0][0],l_off[0][1],
					l_edge_const[0],
					/*l_H_r_P[0],*/ll_H_r_P[0],
					/*l_H_r_Q[0],*/ll_H_r_Q[0],
					node_H_r_PQ_0,
					begin[0],
					end[0]
			);
	krnl_se_pipe(
					l_Vm_dst[1],
					l_Va_dst[1],
					l_Vm_src[0][1],l_Vm_src[1][1],l_Vm_src[2][1],l_Vm_src[3][1],
					l_Va_src[0][1],l_Va_src[1][1],l_Va_src[2][1],l_Va_src[3][1],
					begin[0],begin[1],begin[2],begin[3],
					l_node_const[1],
					l_off[1][0],l_off[1][1],
					l_edge_const[1],
					/*l_H_r_P[1],*/ll_H_r_P[1],
					/*l_H_r_Q[1],*/ll_H_r_Q[1],
					node_H_r_PQ_1,
					begin[1],
					end[1]
			);
	krnl_se_pipe(
					l_Vm_dst[2],
					l_Va_dst[2],
					l_Vm_src[0][2],l_Vm_src[1][2],l_Vm_src[2][2],l_Vm_src[3][2],
					l_Va_src[0][2],l_Va_src[1][2],l_Va_src[2][2],l_Va_src[3][2],
					begin[0],begin[1],begin[2],begin[3],
					l_node_const[2],
					l_off[2][0],l_off[2][1],
					l_edge_const[2],
					/*l_H_r_P[2],*/ll_H_r_P[2],
					/*l_H_r_Q[2],*/ll_H_r_Q[2],
					node_H_r_PQ_2,
					begin[2],
					end[2]
			);
	krnl_se_pipe(
					l_Vm_dst[3],
					l_Va_dst[3],
					l_Vm_src[0][3],l_Vm_src[1][3],l_Vm_src[2][3],l_Vm_src[3][3],
					l_Va_src[0][3],l_Va_src[1][3],l_Va_src[2][3],l_Va_src[3][3],
					begin[0],begin[1],begin[2],begin[3],
					l_node_const[3],
					l_off[3][0],l_off[3][1],
					l_edge_const[3],
					/*l_H_r_P[3],*/ll_H_r_P[3],
					/*l_H_r_Q[3],*/ll_H_r_Q[3],
					node_H_r_PQ_3,
					begin[3],
					end[3]
			);



/*if(H_r_P_flag){
	H_r_P_bank_memcpy(node_H_r_P_0,l_H_r_P[0],begin_0,end_0);
	H_r_P_bank_memcpy(node_H_r_P_1,l_H_r_P[1],begin_1,end_1);
	H_r_P_bank_memcpy(node_H_r_P_2,l_H_r_P[2],begin_2,end_2);
	H_r_P_bank_memcpy(node_H_r_P_3,l_H_r_P[3],begin_3,end_3);
}
if(H_r_Q_flag){
	H_r_Q_bank_memcpy(node_H_r_Q_0,l_H_r_Q[0],begin_0,end_0);
	H_r_Q_bank_memcpy(node_H_r_Q_1,l_H_r_Q[1],begin_1,end_1);
	H_r_Q_bank_memcpy(node_H_r_Q_2,l_H_r_Q[2],begin_2,end_2);
	H_r_Q_bank_memcpy(node_H_r_Q_3,l_H_r_Q[3],begin_3,end_3);
}*/
	//H_r_PQ_bank_memcpy(node_H_r_PQ_0,ll_H_r_P[0],ll_H_r_Q[0],begin_0,end_0);
	//H_r_PQ_bank_memcpy(node_H_r_PQ_1,ll_H_r_P[1],ll_H_r_Q[1],begin_1,end_1);
	//H_r_PQ_bank_memcpy(node_H_r_PQ_2,ll_H_r_P[2],ll_H_r_Q[2],begin_2,end_2);
	//H_r_PQ_bank_memcpy(node_H_r_PQ_3,ll_H_r_P[3],ll_H_r_Q[3],begin_3,end_3);

}
}

