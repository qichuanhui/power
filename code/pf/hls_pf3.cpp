#include "hls_math.h"
#include <ap_int.h>
#include <string.h>

struct Node_PF_PQ{
	double P,Q;
};

struct Node_PF_V{
	double Vm,Va;
};
struct Edge_PF{
	double G,B;
    int src;
};
struct PF_PQ{
	double deltaP;
	double deltaQ;
};


#define EDGE_BLOCK_SIZE 4
struct Edge_PF_Block{
	int src[EDGE_BLOCK_SIZE];
	double G[EDGE_BLOCK_SIZE];
	double B[EDGE_BLOCK_SIZE];
};

#define MAX_NODE_NUM (12*1024)
#define MAX_EDGE_NUM (48*1024)


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

void init_bank_memcpy(const TYPE* const  node_P,double*  l_P,
				const TYPE* const  node_Q,double*  l_Q,
				const TYPE* const  node_type,int*  l_type,
				const TYPE* const  edge_off,int*  l_off_0,int*  l_off_1,
                const TYPE* const  edge_src,int*  l_src,
				const TYPE* const  edge_G, double*  l_G,
				const TYPE* const  edge_B,double*  l_B,
				const int begin,
				const int end
				){
	//#pragma HLS INLINE
	#pragma HLS allocation instances=init_bank_memcpy limit=4 function
	const int node_num=end-begin;

	for(int i=0;i<node_num;i+=8){
		#pragma HLS PIPELINE II=1

		TYPE tmp=node_P[i/8];
		for(int j=0;j<8;j++){
			#pragma HLS unroll factor=8

			long data_l;
			data_l=tmp.range(64*j+63,64*j).to_long();

			union Int64_or_Double data_tmp;
			data_tmp.data_i=data_l;

			double data;
			data=data_tmp.data_d;
			l_P[i+j]=data;
		}
	}

	for(int i=0;i<node_num;i+=8){
		#pragma HLS PIPELINE II=1

		TYPE tmp=node_Q[i/8];
		for(int j=0;j<8;j++){
			#pragma HLS unroll factor=8

			long data_l;
			data_l=tmp.range(64*j+63,64*j).to_long();

			union Int64_or_Double data_tmp;
			data_tmp.data_i=data_l;

			double data;
			data=data_tmp.data_d;
			l_Q[i+j]=data;
		}
	}


	for(int i=0;i<node_num;i+=16){
		#pragma HLS PIPELINE II=1

		TYPE tmp=node_type[i/16];
		for(int j=0;j<16;j++){
			#pragma HLS unroll factor=16

			int data;
			data=tmp.range(32*j+31,32*j).to_int();

			l_type[i+j]=data;
		}
	}
	const int off_length=node_num+1;
	for(int i=0;i<off_length;i+=16){
		#pragma HLS PIPELINE II=2

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

		TYPE tmp=edge_src[i/16];
		for(int j=0;j<16;j++){
			#pragma HLS unroll factor=16

			int data;
			data=tmp.range(32*j+31,32*j).to_int();

			l_src[i+j]=data;
		}
	}
	for(int i=0;i<edge_num;i+=16){
		#pragma HLS PIPELINE II=2

		TYPE tmp=edge_G[i/8];
		TYPE tmp2=edge_G[i/8+1];
		double data[16];
		for(int j=0;j<8;j++){
			#pragma HLS unroll factor=8

			long data_l;
			data_l=tmp.range(64*j+63,64*j).to_long();

			union Int64_or_Double data_tmp;
			data_tmp.data_i=data_l;

			long data_l2;
			data_l2=tmp2.range(64*j+63,64*j).to_long();

			union Int64_or_Double data_tmp2;
			data_tmp2.data_i=data_l2;

			data[j]=data_tmp.data_d;
			data[8+j]=data_tmp2.data_d;
		}
		for(int j=0;j<16;j++){
			#pragma HLS unroll factor=16
			l_G[i+j]=data[j];
		}
		/*for(int j=0;j<8;j++){
			#pragma HLS unroll factor=8

			long data_l;
			data_l=tmp.range(64*j+63,64*j).to_long();

			union Int64_or_Double data_tmp;
			data_tmp.data_i=data_l;

			double data;
			data=data_tmp.data_d;
			l_G[i+j]=data;
		}*/
	}

	for(int i=0;i<edge_num;i+=16){
		#pragma HLS PIPELINE II=2

		TYPE tmp=edge_B[i/8];
		TYPE tmp2=edge_B[i/8+1];
		double data[16];
		for(int j=0;j<8;j++){
			#pragma HLS unroll factor=8

			long data_l;
			data_l=tmp.range(64*j+63,64*j).to_long();

			union Int64_or_Double data_tmp;
			data_tmp.data_i=data_l;

			long data_l2;
			data_l2=tmp2.range(64*j+63,64*j).to_long();

			union Int64_or_Double data_tmp2;
			data_tmp2.data_i=data_l2;

			data[j]=data_tmp.data_d;
			data[8+j]=data_tmp2.data_d;
		}

		for(int j=0;j<16;j++){
			#pragma HLS unroll factor=16
			l_B[i+j]=data[j];
		}
		/*for(int j=0;j<8;j++){
			#pragma HLS unroll factor=8

			long data_l;
			data_l=tmp.range(64*j+63,64*j).to_long();

			union Int64_or_Double data_tmp;
			data_tmp.data_i=data_l;

			double data;
			data=data_tmp.data_d;
			l_B[i+j]=data;
		}*/
	}
}

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
void deltaP_bank_memcpy(TYPE*  node_deltaP,const double* const  l_deltaP,
				const int begin,
				const int end
				){
	//#pragma HLS INLINE
#pragma HLS allocation instances=deltaP_bank_memcpy limit=4 function
	const int node_num=end-begin;

	for(int i=0;i<node_num;i+=8){
		#pragma HLS PIPELINE II=1

		TYPE tmp;
		for(int j=0;j<8;j++){
			#pragma HLS unroll factor=8

			double data=l_deltaP[i+j];
			union Int64_or_Double data_tmp;
			data_tmp.data_d=data;
			long data_l;
			data_l=data_tmp.data_i;
			tmp.range(64*j+63,64*j)=data_l;
		}
		node_deltaP[i/8]=tmp;
	}
}
void deltaQ_bank_memcpy(TYPE*  node_deltaQ,const double* const  l_deltaQ,
				const int begin,
				const int end
				){
	//#pragma HLS INLINE
#pragma HLS allocation instances=deltaQ_bank_memcpy limit=4 function
	const int node_num=end-begin;

	for(int i=0;i<node_num;i+=8){
		#pragma HLS PIPELINE II=1

		TYPE tmp;
		for(int j=0;j<8;j++){
			#pragma HLS unroll factor=8

			double data=l_deltaQ[i+j];
			union Int64_or_Double data_tmp;
			data_tmp.data_d=data;
			long data_l;
			data_l=data_tmp.data_i;
			tmp.range(64*j+63,64*j)=data_l;
		}
		node_deltaQ[i/8]=tmp;
	}
}

#define F2L_PRECISION (1024*1024*256)

void krnl_pf_pipe(
		const double* const  l_P,
		const double* const  l_Q,
		const double* const  l_Vm_dst,
		const double* const  l_Va_dst,
		const double  l_Vm_src_0[2][MAX_NODE_NUM/4],const double  l_Vm_src_1[2][MAX_NODE_NUM/4],const double  l_Vm_src_2[2][MAX_NODE_NUM/4],const double  l_Vm_src_3[2][MAX_NODE_NUM/4],
		const double  l_Va_src_0[2][MAX_NODE_NUM/4],const double  l_Va_src_1[2][MAX_NODE_NUM/4],const double  l_Va_src_2[2][MAX_NODE_NUM/4],const double  l_Va_src_3[2][MAX_NODE_NUM/4],
		const int l_src_offset_0,const int l_src_offset_1,const int l_src_offset_2,const int l_src_offset_3,
		const int* const  l_type,
		const int* const  l_off_0,const int* const  l_off_1,
		const int* const  l_src,
		const double* const  l_G,
		const double* const  l_B,
		double*  l_deltaP,long*  ll_deltaP,
		double*  l_deltaQ,long*  ll_deltaQ,
		double*  l_max_deltaP,
		double*  l_max_deltaQ,
		const int begin,
		const int end
		){
//#pragma HLS INLINE
#pragma HLS allocation instances=krnl_pf_pipe limit=4 function

	const int l_node_num=end-begin;

	for(int i=0;i<l_node_num;i+=8){
		#pragma HLS PIPELINE II=1
		for(int j=0;j<8;j++){
			#pragma HLS UNROLL factor=8
			ll_deltaP[i+j]=0;
			ll_deltaQ[i+j]=0;
		}
	}

	double p_max_deltaP=0;
	double p_max_deltaQ=0;
	const int edge_offset=l_off_0[0];

	const int edge_num=l_off_0[l_node_num]-edge_offset;
	int edge_count=0;
	int dst_id=0;
	while(edge_count<edge_num&&dst_id<l_node_num){
		#pragma HLS PIPELINE II=1
/*if((dst_id<20||dst_id>2500)&&begin==0)
	printf("dst_id=%d,dst_id+1=%d,edge_count=%d,l_off_0[%d]=%d,l_off_0[%d]=%d,l_off_0[%d]=%d\n",
			dst_id,dst_id+1,edge_count,dst_id,l_off_0[dst_id],
			dst_id+1,l_off_0[dst_id+1],dst_id+2,l_off_0[dst_id+2]);*/

		int buf_off_0,buf_off_1,buf_off_2,buf_off_3;
		int buf_off_4,buf_off_5,buf_off_6,buf_off_7;
		int buf_off_8;

		if(dst_id+8<l_node_num){
			buf_off_8=l_off_1[dst_id+8]-edge_offset;
		}
		else{
			buf_off_8=edge_num;
		}

		/*int buf_off_0=l_off_0[dst_id+0]-edge_offset;
		int buf_off_1=l_off_1[dst_id+1]-edge_offset;
		int buf_off_2=l_off_2[dst_id+2]-edge_offset;*/

		const int old_dst_id=dst_id;
		const int old_edge_count=edge_count;

		if(buf_off_8>edge_count+16){
			dst_id+=0;edge_count+=16;
		}
		else if(buf_off_8==edge_count+16){
			dst_id+=8;edge_count+=16;
		}
		else{
			dst_id+=8;
		}

		for(int i=0;i<8;i++){
			#pragma HLS UNROLL factor=8
			int tmp_off;//=l_off_0[dst_id*8/8+i]-edge_offset;
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
		int dst_type[8];
		for(int j=0;j<8;j++){
			#pragma HLS UNROLL factor=8
			dst_Vm[j]=l_Vm_dst[old_dst_id/8*8+j];
			dst_Va[j]=l_Va_dst[old_dst_id/8*8+j];
			dst_type[j]=l_type[old_dst_id/8*8+j];
		}

		
    	 // 在顶层作用域定义数组（假设类型为 int）
		long buf_ll_deltaP_1[16];  // 索引0~15
		long buf_ll_deltaQ_1[16];  // 索引0~15
		#pragma HLS ARRAY_PARTITION variable=buf_ll_deltaP_1 complete dim=1
		#pragma HLS ARRAY_PARTITION variable=buf_ll_deltaQ_1 complete dim=1

		long buf_ll_deltaP_5[16];  // 索引0~15
		long buf_ll_deltaQ_5[16];  // 索引0~15
		#pragma HLS ARRAY_PARTITION variable=buf_ll_deltaP_5 complete dim=1
		#pragma HLS ARRAY_PARTITION variable=buf_ll_deltaQ_5 complete dim=1
		for(int j=0;j<16;j++){
			#pragma HLS UNROLL factor=16

			int src_id=l_src[old_edge_count/16*16+j];
			double src_G=l_G[old_edge_count/16*16+j];
			double src_B=l_B[old_edge_count/16*16+j];
			double src_Vm=0;
			double src_Va=0;
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
			if(old_edge_count+j<buf_off_1) delta_Va=dst_Va[0]-src_Va;
			else if(old_edge_count+j<buf_off_2) delta_Va=dst_Va[1]-src_Va;
			else if(old_edge_count+j<buf_off_3) delta_Va=dst_Va[2]-src_Va;
			else if(old_edge_count+j<buf_off_4) delta_Va=dst_Va[3]-src_Va;
			else if(old_edge_count+j<buf_off_5) delta_Va=dst_Va[4]-src_Va;
			else if(old_edge_count+j<buf_off_6) delta_Va=dst_Va[5]-src_Va;
			else if(old_edge_count+j<buf_off_7) delta_Va=dst_Va[6]-src_Va;
			else delta_Va=dst_Va[7]-src_Va;
			double cosine_array=cos(delta_Va);//hls::cordic::
			double sine_array=sin(delta_Va);

			double d_deltaP	=- src_Vm*(src_G*cosine_array +src_B*sine_array);
			double d_deltaQ =- src_Vm*(src_G*sine_array - src_B*cosine_array);
			//buf_ll_deltaP_1[j]=(long)(d_deltaP*F2L_PRECISION);
			//buf_ll_deltaQ_1[j]=(long)(d_deltaQ*F2L_PRECISION);
			long tmp_ll_deltaP_1=(long)(d_deltaP*F2L_PRECISION);
			long tmp_ll_deltaQ_1=(long)(d_deltaQ*F2L_PRECISION);
			switch(j){
			case 0:
				buf_ll_deltaP_1[0]=	tmp_ll_deltaP_1;buf_ll_deltaQ_1[1]=tmp_ll_deltaQ_1;break;
			case 1:
				buf_ll_deltaP_1[1]=	tmp_ll_deltaP_1;buf_ll_deltaQ_1[1]=tmp_ll_deltaQ_1;break;
			case 2:
				buf_ll_deltaP_1[2]=	tmp_ll_deltaP_1;buf_ll_deltaQ_1[1]=tmp_ll_deltaQ_1;break;
			case 3:
				buf_ll_deltaP_1[3]=	tmp_ll_deltaP_1;buf_ll_deltaQ_1[1]=tmp_ll_deltaQ_1;break;
			case 4:
				buf_ll_deltaP_1[4]=	tmp_ll_deltaP_1;buf_ll_deltaQ_1[1]=tmp_ll_deltaQ_1;break;
			case 5:
				buf_ll_deltaP_1[5]=	tmp_ll_deltaP_1;buf_ll_deltaQ_1[1]=tmp_ll_deltaQ_1;break;
			case 6:
				buf_ll_deltaP_1[6]=	tmp_ll_deltaP_1;buf_ll_deltaQ_1[1]=tmp_ll_deltaQ_1;break;
			case 7:
				buf_ll_deltaP_1[7]=	tmp_ll_deltaP_1;buf_ll_deltaQ_1[1]=tmp_ll_deltaQ_1;break;
			case 8:
				buf_ll_deltaP_1[8]=	tmp_ll_deltaP_1;buf_ll_deltaQ_1[1]=tmp_ll_deltaQ_1;break;
			case 9:
				buf_ll_deltaP_1[9]=	tmp_ll_deltaP_1;buf_ll_deltaQ_1[1]=tmp_ll_deltaQ_1;break;
			case 10:
				buf_ll_deltaP_1[10]=tmp_ll_deltaP_1;buf_ll_deltaQ_1[1]=tmp_ll_deltaQ_1;break;
			case 11:
				buf_ll_deltaP_1[11]=tmp_ll_deltaP_1;buf_ll_deltaQ_1[1]=tmp_ll_deltaQ_1;break;
			case 12:
				buf_ll_deltaP_1[12]=tmp_ll_deltaP_1;buf_ll_deltaQ_1[1]=tmp_ll_deltaQ_1;break;
			case 13:
				buf_ll_deltaP_1[13]=tmp_ll_deltaP_1;buf_ll_deltaQ_1[1]=tmp_ll_deltaQ_1;break;
			case 14:
				buf_ll_deltaP_1[14]=tmp_ll_deltaP_1;buf_ll_deltaQ_1[1]=tmp_ll_deltaQ_1;break;
			case 15:
				buf_ll_deltaP_1[15]=tmp_ll_deltaP_1;buf_ll_deltaQ_1[1]=tmp_ll_deltaQ_1;break;
			}
		}

// 原始输入寄存器级
long in_p[16], in_q[16];
#pragma HLS ARRAY_PARTITION variable=in_p complete
#pragma HLS ARRAY_PARTITION variable=in_q complete
#pragma HLS BIND_STORAGE variable=in_p type=ram_1p impl=lutram
#pragma HLS BIND_STORAGE variable=in_q type=ram_1p impl=lutram
in_p[0] = buf_ll_deltaP_1[0];
in_p[1] = buf_ll_deltaP_1[1];
in_p[2] = buf_ll_deltaP_1[2];
in_p[3] = buf_ll_deltaP_1[3];
in_p[4] = buf_ll_deltaP_1[4];
in_p[5] = buf_ll_deltaP_1[5];
in_p[6] = buf_ll_deltaP_1[6];
in_p[7] = buf_ll_deltaP_1[7];
in_p[8] = buf_ll_deltaP_1[8];
in_p[9] = buf_ll_deltaP_1[9];
in_p[10] = buf_ll_deltaP_1[10];
in_p[11] = buf_ll_deltaP_1[11];
in_p[12] = buf_ll_deltaP_1[12];
in_p[13] = buf_ll_deltaP_1[13];
in_p[14] = buf_ll_deltaP_1[14];
in_p[15] = buf_ll_deltaP_1[15];

in_q[0] = buf_ll_deltaQ_1[0];
in_q[1] = buf_ll_deltaQ_1[1];
in_q[2] = buf_ll_deltaQ_1[2];
in_q[3] = buf_ll_deltaQ_1[3];
in_q[4] = buf_ll_deltaQ_1[4];
in_q[5] = buf_ll_deltaQ_1[5];
in_q[6] = buf_ll_deltaQ_1[6];
in_q[7] = buf_ll_deltaQ_1[7];
in_q[8] = buf_ll_deltaQ_1[8];
in_q[9] = buf_ll_deltaQ_1[9];
in_q[10] = buf_ll_deltaQ_1[10];
in_q[11] = buf_ll_deltaQ_1[11];
in_q[12] = buf_ll_deltaQ_1[12];
in_q[13] = buf_ll_deltaQ_1[13];
in_q[14] = buf_ll_deltaQ_1[14];
in_q[15] = buf_ll_deltaQ_1[15];

// 第一级计算
long stage1_p[16], stage1_q[16];
#pragma HLS ARRAY_PARTITION variable=stage1_p complete
#pragma HLS ARRAY_PARTITION variable=stage1_q complete
#pragma HLS BIND_STORAGE variable=stage1_p type=ram_1p impl=lutram
#pragma HLS BIND_STORAGE variable=stage1_q type=ram_1p impl=lutram
stage1_p[0] = in_p[0];
stage1_q[0] = in_q[0];
stage1_p[2] = in_p[2];
stage1_q[2] = in_q[2];
stage1_p[4] = in_p[4];
stage1_q[4] = in_q[4];
stage1_p[6] = in_p[6];
stage1_q[6] = in_q[6];
stage1_p[8] = in_p[8];
stage1_q[8] = in_q[8];
stage1_p[10] = in_p[10];
stage1_q[10] = in_q[10];
stage1_p[12] = in_p[12];
stage1_q[12] = in_q[12];
stage1_p[14] = in_p[14];
stage1_q[14] = in_q[14];

stage1_p[1] = in_p[1] + in_p[0];
stage1_q[1] = in_q[1] + in_q[0];
stage1_p[3] = in_p[3] + in_p[2];
stage1_q[3] = in_q[3] + in_q[2];
stage1_p[5] = in_p[5] + in_p[4];
stage1_q[5] = in_q[5] + in_q[4];
stage1_p[7] = in_p[7] + in_p[6];
stage1_q[7] = in_q[7] + in_q[6];
stage1_p[9] = in_p[9] + in_p[8];
stage1_q[9] = in_q[9] + in_q[8];
stage1_p[11] = in_p[11] + in_p[10];
stage1_q[11] = in_q[11] + in_q[10];
stage1_p[13] = in_p[13] + in_p[12];
stage1_q[13] = in_q[13] + in_q[12];
stage1_p[15] = in_p[15] + in_p[14];
stage1_q[15] = in_q[15] + in_q[14];

// 第二级计算
long stage2_p[16], stage2_q[16];
#pragma HLS ARRAY_PARTITION variable=stage2_p complete
#pragma HLS ARRAY_PARTITION variable=stage2_q complete
#pragma HLS BIND_STORAGE variable=stage2_p type=ram_1p impl=lutram
#pragma HLS BIND_STORAGE variable=stage2_q type=ram_1p impl=lutram
stage2_p[0] = stage1_p[0];
stage2_q[0] = stage1_q[0];
stage2_p[1] = stage1_p[1];
stage2_q[1] = stage1_q[1];
stage2_p[4] = stage1_p[4];
stage2_q[4] = stage1_q[4];
stage2_p[5] = stage1_p[5];
stage2_q[5] = stage1_q[5];
stage2_p[8] = stage1_p[8];
stage2_q[8] = stage1_q[8];
stage2_p[9] = stage1_p[9];
stage2_q[9] = stage1_q[9];
stage2_p[12] = stage1_p[12];
stage2_q[12] = stage1_q[12];
stage2_p[13] = stage1_p[13];
stage2_q[13] = stage1_q[13];

stage2_p[2] = stage1_p[2] + stage1_p[1];
stage2_q[2] = stage1_q[2] + stage1_q[1];
stage2_p[6] = stage1_p[6] + stage1_p[5];
stage2_q[6] = stage1_q[6] + stage1_q[5];
stage2_p[10] = stage1_p[10] + stage1_p[9];
stage2_q[10] = stage1_q[10] + stage1_q[9];
stage2_p[14] = stage1_p[14] + stage1_p[13];
stage2_q[14] = stage1_q[14] + stage1_q[13];

stage2_p[3] = stage1_p[3] + stage1_p[1];
stage2_q[3] = stage1_q[3] + stage1_q[1];
stage2_p[7] = stage1_p[7] + stage1_p[5];
stage2_q[7] = stage1_q[7] + stage1_q[5];
stage2_p[11] = stage1_p[11] + stage1_p[9];
stage2_q[11] = stage1_q[11] + stage1_q[9];
stage2_p[15] = stage1_p[15] + stage1_p[13];
stage2_q[15] = stage1_q[15] + stage1_q[13];

// 第三级计算
long stage3_p[16], stage3_q[16];
#pragma HLS ARRAY_PARTITION variable=stage3_p complete
#pragma HLS ARRAY_PARTITION variable=stage3_q complete
#pragma HLS BIND_STORAGE variable=stage3_p type=ram_1p impl=lutram
#pragma HLS BIND_STORAGE variable=stage3_q type=ram_1p impl=lutram
stage3_p[0] = stage2_p[0];
stage3_q[0] = stage2_q[0];
stage3_p[1] = stage2_p[1];
stage3_q[1] = stage2_q[1];
stage3_p[2] = stage2_p[2];
stage3_q[2] = stage2_q[2];
stage3_p[3] = stage2_p[3];
stage3_q[3] = stage2_q[3];
stage3_p[8] = stage2_p[8];
stage3_q[8] = stage2_q[8];
stage3_p[9] = stage2_p[9];
stage3_q[9] = stage2_q[9];
stage3_p[10] = stage2_p[10];
stage3_q[10] = stage2_q[10];
stage3_p[11] = stage2_p[11];
stage3_q[11] = stage2_q[11];

stage3_p[4] = stage2_p[4] + stage2_p[3];
stage3_q[4] = stage2_q[4] + stage2_q[3];
stage3_p[5] = stage2_p[5] + stage2_p[3];
stage3_q[5] = stage2_q[5] + stage2_q[3];
stage3_p[6] = stage2_p[6] + stage2_p[3];
stage3_q[6] = stage2_q[6] + stage2_q[3];
stage3_p[7] = stage2_p[7] + stage2_p[3];
stage3_q[7] = stage2_q[7] + stage2_q[3];

stage3_p[12] = stage2_p[12] + stage2_p[11];
stage3_q[12] = stage2_q[12] + stage2_q[11];
stage3_p[13] = stage2_p[13] + stage2_p[11];
stage3_q[13] = stage2_q[13] + stage2_q[11];
stage3_p[14] = stage2_p[14] + stage2_p[11];
stage3_q[14] = stage2_q[14] + stage2_q[11];
stage3_p[15] = stage2_p[15] + stage2_p[11];
stage3_q[15] = stage2_q[15] + stage2_q[11];

// 第四级计算
// 第四级计算 - 优化版本
long stage4_p[16], stage4_q[16];
#pragma HLS ARRAY_PARTITION variable=stage4_p complete
#pragma HLS ARRAY_PARTITION variable=stage4_q complete
#pragma HLS BIND_STORAGE variable=stage4_p type=ram_1p impl=lutram
#pragma HLS BIND_STORAGE variable=stage4_q type=ram_1p impl=lutram
// 步骤1: 0-7位置直接复制
for(int i = 0; i < 8; i++) {
    #pragma HLS UNROLL
    stage4_p[i] = stage3_p[i];
    stage4_q[i] = stage3_q[i];
}

// 步骤2: 创建渐进式加法链（消除stage3_p[7]高扇出）
stage4_p[8] = stage3_p[8] + stage3_p[7];  // 第一个使用stage3_p[7]
stage4_q[8] = stage3_q[8] + stage3_q[7];

// 后续位置依次基于前一个结果
stage4_p[9] = stage3_p[9] + stage4_p[8];  // 用stage4_p[8]而非stage3_p[7]
stage4_q[9] = stage3_q[9] + stage4_q[8];

stage4_p[10] = stage3_p[10] + stage4_p[9];
stage4_q[10] = stage3_q[10] + stage4_q[9];

stage4_p[11] = stage3_p[11] + stage4_p[10];
stage4_q[11] = stage3_q[11] + stage4_q[10];

stage4_p[12] = stage3_p[12] + stage4_p[11];
stage4_q[12] = stage3_q[12] + stage4_q[11];

stage4_p[13] = stage3_p[13] + stage4_p[12];
stage4_q[13] = stage3_q[13] + stage4_q[12];

stage4_p[14] = stage3_p[14] + stage4_p[13];
stage4_q[14] = stage3_q[14] + stage4_q[13];

stage4_p[15] = stage3_p[15] + stage4_p[14];
stage4_q[15] = stage3_q[15] + stage4_q[14];

// 第五级计算 (最终输出)
buf_ll_deltaP_5[0] = stage4_p[0];
buf_ll_deltaP_5[1] = stage4_p[1];
buf_ll_deltaP_5[2] = stage4_p[2];
buf_ll_deltaP_5[3] = stage4_p[3];
buf_ll_deltaP_5[4] = stage4_p[4];
buf_ll_deltaP_5[5] = stage4_p[5];
buf_ll_deltaP_5[6] = stage4_p[6];
buf_ll_deltaP_5[7] = stage4_p[7];
buf_ll_deltaP_5[8] = stage4_p[8];
buf_ll_deltaP_5[9] = stage4_p[9];
buf_ll_deltaP_5[10] = stage4_p[10];
buf_ll_deltaP_5[11] = stage4_p[11];
buf_ll_deltaP_5[12] = stage4_p[12];
buf_ll_deltaP_5[13] = stage4_p[13];
buf_ll_deltaP_5[14] = stage4_p[14];
buf_ll_deltaP_5[15] = stage4_p[15];

buf_ll_deltaQ_5[0]   = stage4_q[0];
buf_ll_deltaQ_5[1]   = stage4_q[1];
buf_ll_deltaQ_5[2]   = stage4_q[2];
buf_ll_deltaQ_5[3]   = stage4_q[3];
buf_ll_deltaQ_5[4]   = stage4_q[4];
buf_ll_deltaQ_5[5]   = stage4_q[5];
buf_ll_deltaQ_5[6]   = stage4_q[6];
buf_ll_deltaQ_5[7]   = stage4_q[7];
buf_ll_deltaQ_5[8]   = stage4_q[8];
buf_ll_deltaQ_5[9]   = stage4_q[9];
buf_ll_deltaQ_5[10]   = stage4_q[10];
buf_ll_deltaQ_5[11]   = stage4_q[11];
buf_ll_deltaQ_5[12]   = stage4_q[12];
buf_ll_deltaQ_5[13]   = stage4_q[13];
buf_ll_deltaQ_5[14]   = stage4_q[14];
buf_ll_deltaQ_5[15]   = stage4_q[15];


		/*if(begin==0&&old_edge_count<50){
			for(int j=0;j<8;j++){
				printf("buf_ll_deltaP_1[%d]=%lf\tbuf_ll_deltaP_4[%d]=%lf\n",
					old_edge_count+j,(double)buf_ll_deltaP_1[j]/(double)F2L_PRECISION,
					old_edge_count+j,(double)buf_ll_deltaP_4[j]/(double)F2L_PRECISION);
			}
			printf("ll_deltaP[0]=%lf\tll_deltaP[1]=%lf\tll_deltaP[2]=%lf\n",
					(double)ll_deltaP[0]/(double)F2L_PRECISION,
					(double)ll_deltaP[1]/(double)F2L_PRECISION,
					(double)ll_deltaP[2]/(double)F2L_PRECISION);
		}*/
		int locat_0,locat_1,locat_2,locat_3,locat_4,locat_5,locat_6,locat_7,locat_8;
		long deltaP_value_0,deltaP_value_1,deltaP_value_2,deltaP_value_3;
		long deltaP_value_4,deltaP_value_5,deltaP_value_6,deltaP_value_7;
		long deltaP_value_8;
		long deltaQ_value_0,deltaQ_value_1,deltaQ_value_2,deltaQ_value_3;
		long deltaQ_value_4,deltaQ_value_5,deltaQ_value_6,deltaQ_value_7;
		long deltaQ_value_8;


#define get_predix_sum(locat_0,deltaP_value_0,deltaQ_value_0)\
		switch(locat_0){\
			case 0:\
				deltaP_value_0=buf_ll_deltaP_5[0];	deltaQ_value_0=buf_ll_deltaQ_5[0];break;\
			case 1:\
				deltaP_value_0=buf_ll_deltaP_5[1];	deltaQ_value_0=buf_ll_deltaQ_5[1];break;\
			case 2:\
				deltaP_value_0=buf_ll_deltaP_5[2];	deltaQ_value_0=buf_ll_deltaQ_5[2];break;\
			case 3:\
				deltaP_value_0=buf_ll_deltaP_5[3];	deltaQ_value_0=buf_ll_deltaQ_5[3];break;\
			case 4:\
				deltaP_value_0=buf_ll_deltaP_5[4];	deltaQ_value_0=buf_ll_deltaQ_5[4];break;\
			case 5:\
				deltaP_value_0=buf_ll_deltaP_5[5];	deltaQ_value_0=buf_ll_deltaQ_5[5];break;\
			case 6:\
				deltaP_value_0=buf_ll_deltaP_5[6];	deltaQ_value_0=buf_ll_deltaQ_5[6];break;\
			case 7:\
				deltaP_value_0=buf_ll_deltaP_5[7];	deltaQ_value_0=buf_ll_deltaQ_5[7];break;\
			case 8:\
				deltaP_value_0=buf_ll_deltaP_5[8];	deltaQ_value_0=buf_ll_deltaQ_5[8];break;\
			case 9:\
				deltaP_value_0=buf_ll_deltaP_5[9];	deltaQ_value_0=buf_ll_deltaQ_5[9];break;\
			case 10:\
				deltaP_value_0=buf_ll_deltaP_5[10];	deltaQ_value_0=buf_ll_deltaQ_5[10];break;\
			case 11:\
				deltaP_value_0=buf_ll_deltaP_5[11];	deltaQ_value_0=buf_ll_deltaQ_5[11];break;\
			case 12:\
				deltaP_value_0=buf_ll_deltaP_5[12];	deltaQ_value_0=buf_ll_deltaQ_5[12];break;\
			case 13:\
				deltaP_value_0=buf_ll_deltaP_5[13];	deltaQ_value_0=buf_ll_deltaQ_5[13];break;\
			case 14:\
				deltaP_value_0=buf_ll_deltaP_5[14];	deltaQ_value_0=buf_ll_deltaQ_5[14];break;\
			case 15:\
				deltaP_value_0=buf_ll_deltaP_5[15];	deltaQ_value_0=buf_ll_deltaQ_5[15];break;\
		}

#define get_predix_sum_2(buf_off_0,locat_0,deltaP_value_0,deltaQ_value_0) \
		if(buf_off_0-old_edge_count-1<0){\
			deltaP_value_0=0;deltaQ_value_0=0;\
		}\
		else{\
			if(buf_off_0-old_edge_count-1>15){\
				locat_0=15;\
			}\
			else{\
				locat_0=buf_off_0-old_edge_count-1;\
			}\
			get_predix_sum(locat_0,deltaP_value_0,deltaQ_value_0)\
		}

		get_predix_sum_2(buf_off_0,locat_0,deltaP_value_0,deltaQ_value_0)
		get_predix_sum_2(buf_off_1,locat_1,deltaP_value_1,deltaQ_value_1)
		get_predix_sum_2(buf_off_2,locat_2,deltaP_value_2,deltaQ_value_2)
		get_predix_sum_2(buf_off_3,locat_3,deltaP_value_3,deltaQ_value_3)
		get_predix_sum_2(buf_off_4,locat_4,deltaP_value_4,deltaQ_value_4)
		get_predix_sum_2(buf_off_5,locat_5,deltaP_value_5,deltaQ_value_5)
		get_predix_sum_2(buf_off_6,locat_6,deltaP_value_6,deltaQ_value_6)
		get_predix_sum_2(buf_off_7,locat_7,deltaP_value_7,deltaQ_value_7)
		get_predix_sum_2(buf_off_8,locat_8,deltaP_value_8,deltaQ_value_8)

		for(int j=0;j<8;j++){
			#pragma HLS UNROLL factor=8
			long tmpP,tmpQ;
			switch(j){
			case 0:
				tmpP=deltaP_value_1-deltaP_value_0;tmpQ=deltaQ_value_1-deltaQ_value_0;break;
			case 1:
				tmpP=deltaP_value_2-deltaP_value_1;tmpQ=deltaQ_value_2-deltaQ_value_1;break;
			case 2:
				tmpP=deltaP_value_3-deltaP_value_2;tmpQ=deltaQ_value_3-deltaQ_value_2;break;
			case 3:
				tmpP=deltaP_value_4-deltaP_value_3;tmpQ=deltaQ_value_4-deltaQ_value_3;break;
			case 4:
				tmpP=deltaP_value_5-deltaP_value_4;tmpQ=deltaQ_value_5-deltaQ_value_4;break;
			case 5:
				tmpP=deltaP_value_6-deltaP_value_5;tmpQ=deltaQ_value_6-deltaQ_value_5;break;
			case 6:
				tmpP=deltaP_value_7-deltaP_value_6;tmpQ=deltaQ_value_7-deltaQ_value_6;break;
			case 7:
				tmpP=deltaP_value_8-deltaP_value_7;tmpQ=deltaQ_value_8-deltaQ_value_7;break;
			}

			ll_deltaP[old_dst_id/8*8+j]+=tmpP;
			ll_deltaQ[old_dst_id/8*8+j]+=tmpQ;
		}
	}
	for(int i=0;i<l_node_num;i+=16){
			#pragma HLS PIPELINE II=1

			double buf_deltaP_1[16],buf_deltaQ_1[16];
			double buf_deltaP_2[8],buf_deltaQ_2[8];
			double buf_deltaP_3[4],buf_deltaQ_3[4];
			double buf_deltaP_4[2],buf_deltaQ_4[2];
			double buf_deltaP_5[1],buf_deltaQ_5[1];
			for(int j=0;j<16;j++){
				#pragma HLS UNROLL factor=16

				double d_l_deltaP=(double)ll_deltaP[i+j]/(double)F2L_PRECISION;
				double d_l_deltaQ=(double)ll_deltaQ[i+j]/(double)F2L_PRECISION;

				double dst_Vm=l_Vm_dst[i+j];
				int dst_type=l_type[i+j];
				d_l_deltaP+=l_P[i+j]/dst_Vm;
				d_l_deltaQ+=l_Q[i+j]/dst_Vm;
				if(dst_type==3||i+j>=l_node_num) d_l_deltaP=0;
				if(dst_type>=2||i+j>=l_node_num) d_l_deltaQ=0;
				l_deltaP[i+j]=d_l_deltaP;
				l_deltaQ[i+j]=d_l_deltaQ;

				buf_deltaP_1[j]=fabs(d_l_deltaP);
				buf_deltaQ_1[j]=fabs(d_l_deltaQ);
			}
			for(int j=0;j<8;j++){
				#pragma HLS UNROLL factor=8
				buf_deltaP_2[j]=fmax(buf_deltaP_1[j],buf_deltaP_1[8+j]);
				buf_deltaQ_2[j]=fmax(buf_deltaQ_1[j],buf_deltaQ_1[8+j]);
			}
			for(int j=0;j<4;j++){
				#pragma HLS UNROLL factor=4
				buf_deltaP_3[j]=fmax(buf_deltaP_2[j],buf_deltaP_2[4+j]);
				buf_deltaQ_3[j]=fmax(buf_deltaQ_2[j],buf_deltaQ_2[4+j]);
			}
			for(int j=0;j<2;j++){
				#pragma HLS UNROLL factor=2
				buf_deltaP_4[j]=fmax(buf_deltaP_3[j],buf_deltaP_3[2+j]);
				buf_deltaQ_4[j]=fmax(buf_deltaQ_3[j],buf_deltaQ_3[2+j]);
			}

			buf_deltaP_5[0]=fmax(buf_deltaP_4[0],buf_deltaP_4[1]);
			buf_deltaQ_5[0]=fmax(buf_deltaQ_4[0],buf_deltaQ_4[1]);

			p_max_deltaP=fmax(p_max_deltaP,buf_deltaP_5[0]);
			p_max_deltaQ=fmax(p_max_deltaQ,buf_deltaQ_5[0]);
	}

	l_max_deltaP[0]=p_max_deltaP;
	l_max_deltaQ[0]=p_max_deltaQ;
}

extern "C" {
void krnl_pf(
				const TYPE* const  node_P_0,
				const TYPE* const  node_P_1,
				const TYPE* const  node_P_2,
				const TYPE* const  node_P_3,
				const TYPE* const  node_Q_0,
				const TYPE* const  node_Q_1,
				const TYPE* const  node_Q_2,
				const TYPE* const  node_Q_3,
				const TYPE* const  node_Vm_0,
				const TYPE* const  node_Vm_1,
				const TYPE* const  node_Vm_2,
				const TYPE* const  node_Vm_3,
				const TYPE* const  node_Va_0,
				const TYPE* const  node_Va_1,
				const TYPE* const  node_Va_2,
				const TYPE* const  node_Va_3,
				const TYPE* const  node_type_0,
				const TYPE* const  node_type_1,
				const TYPE* const  node_type_2,
				const TYPE* const  node_type_3,
				const TYPE* const  edge_off_0,
				const TYPE* const  edge_off_1,
				const TYPE* const  edge_off_2,
				const TYPE* const  edge_off_3,
                const TYPE* const  edge_src_0,
				const TYPE* const  edge_src_1,
				const TYPE* const  edge_src_2,
				const TYPE* const  edge_src_3,
				const TYPE* const  edge_G_0,
				const TYPE* const  edge_G_1,
				const TYPE* const  edge_G_2,
				const TYPE* const  edge_G_3,
				const TYPE* const  edge_B_0,
				const TYPE* const  edge_B_1,
				const TYPE* const  edge_B_2,
				const TYPE* const  edge_B_3,
				TYPE*  node_deltaP_0,
				TYPE*  node_deltaP_1,
				TYPE*  node_deltaP_2,
				TYPE*  node_deltaP_3,
				TYPE*  node_deltaQ_0,
				TYPE*  node_deltaQ_1,
				TYPE*  node_deltaQ_2,
				TYPE*  node_deltaQ_3,
				TYPE*  max_delta,
				const int begin_0,
				const int begin_1,
				const int begin_2,
				const int begin_3,
				const int end_0,
				const int end_1,
				const int end_2,
				const int end_3,
				const int node_num,
				const bool init_flag,
				const bool Vm_flag,
				const bool Va_flag,
				const bool deltaP_flag,
				const bool deltaQ_flag,
				const bool compute_flag
				)
{
#pragma HLS INTERFACE m_axi port = node_P_0 offset = slave bundle = gmem0 //depth=512
#pragma HLS INTERFACE m_axi port = node_P_1 offset = slave bundle = gmem1 //depth=512
#pragma HLS INTERFACE m_axi port = node_P_2 offset = slave bundle = gmem2 //depth=512
#pragma HLS INTERFACE m_axi port = node_P_3 offset = slave bundle = gmem3 //depth=512

#pragma HLS INTERFACE m_axi port = node_Q_0 offset = slave bundle = gmem0 //depth=512
#pragma HLS INTERFACE m_axi port = node_Q_1 offset = slave bundle = gmem1 //depth=512
#pragma HLS INTERFACE m_axi port = node_Q_2 offset = slave bundle = gmem2 //depth=512
#pragma HLS INTERFACE m_axi port = node_Q_3 offset = slave bundle = gmem3 //depth=512

#pragma HLS INTERFACE m_axi port = node_Vm_0 offset = slave bundle = gmem0 //depth=512
#pragma HLS INTERFACE m_axi port = node_Vm_1 offset = slave bundle = gmem1 //depth=512
#pragma HLS INTERFACE m_axi port = node_Vm_2 offset = slave bundle = gmem2 //depth=512
#pragma HLS INTERFACE m_axi port = node_Vm_3 offset = slave bundle = gmem3 //depth=512

#pragma HLS INTERFACE m_axi port = node_Va_0 offset = slave bundle = gmem0 //depth=512
#pragma HLS INTERFACE m_axi port = node_Va_1 offset = slave bundle = gmem1 //depth=512
#pragma HLS INTERFACE m_axi port = node_Va_2 offset = slave bundle = gmem2 //depth=512
#pragma HLS INTERFACE m_axi port = node_Va_3 offset = slave bundle = gmem3 //depth=512

#pragma HLS INTERFACE m_axi port = node_type_0 offset = slave bundle = gmem0 //depth=512
#pragma HLS INTERFACE m_axi port = node_type_1 offset = slave bundle = gmem1 //depth=512
#pragma HLS INTERFACE m_axi port = node_type_2 offset = slave bundle = gmem2 //depth=512
#pragma HLS INTERFACE m_axi port = node_type_3 offset = slave bundle = gmem3 //depth=512

#pragma HLS INTERFACE m_axi port = edge_off_0 offset = slave bundle = gmem0 //depth=512
#pragma HLS INTERFACE m_axi port = edge_off_1 offset = slave bundle = gmem1 //depth=512
#pragma HLS INTERFACE m_axi port = edge_off_2 offset = slave bundle = gmem2 //depth=512
#pragma HLS INTERFACE m_axi port = edge_off_3 offset = slave bundle = gmem3 //depth=512

#pragma HLS INTERFACE m_axi port = edge_src_0 offset = slave bundle = gmem0 //depth=512
#pragma HLS INTERFACE m_axi port = edge_src_1 offset = slave bundle = gmem1 //depth=512
#pragma HLS INTERFACE m_axi port = edge_src_2 offset = slave bundle = gmem2 //depth=512
#pragma HLS INTERFACE m_axi port = edge_src_3 offset = slave bundle = gmem3 //depth=512

#pragma HLS INTERFACE m_axi port = edge_G_0 offset = slave bundle = gmem0 //depth=512
#pragma HLS INTERFACE m_axi port = edge_G_1 offset = slave bundle = gmem1 //depth=512
#pragma HLS INTERFACE m_axi port = edge_G_2 offset = slave bundle = gmem2 //depth=512
#pragma HLS INTERFACE m_axi port = edge_G_3 offset = slave bundle = gmem3 //depth=512

#pragma HLS INTERFACE m_axi port = edge_B_0 offset = slave bundle = gmem0 //depth=512
#pragma HLS INTERFACE m_axi port = edge_B_1 offset = slave bundle = gmem1 //depth=512
#pragma HLS INTERFACE m_axi port = edge_B_2 offset = slave bundle = gmem2 //depth=512
#pragma HLS INTERFACE m_axi port = edge_B_3 offset = slave bundle = gmem3 //depth=512

#pragma HLS INTERFACE m_axi port = node_deltaP_0 offset = slave bundle = gmem0 //depth=512
#pragma HLS INTERFACE m_axi port = node_deltaP_1 offset = slave bundle = gmem1 //depth=512
#pragma HLS INTERFACE m_axi port = node_deltaP_2 offset = slave bundle = gmem2 //depth=512
#pragma HLS INTERFACE m_axi port = node_deltaP_3 offset = slave bundle = gmem3 //depth=512

#pragma HLS INTERFACE m_axi port = node_deltaQ_0 offset = slave bundle = gmem0 //depth=512
#pragma HLS INTERFACE m_axi port = node_deltaQ_1 offset = slave bundle = gmem1 //depth=512
#pragma HLS INTERFACE m_axi port = node_deltaQ_2 offset = slave bundle = gmem2 //depth=512
#pragma HLS INTERFACE m_axi port = node_deltaQ_3 offset = slave bundle = gmem3 //depth=512

#pragma HLS INTERFACE m_axi port = max_delta offset = slave bundle = gmem0 //depth=512


#pragma HLS INTERFACE s_axilite port=node_P_0 bundle=control
#pragma HLS INTERFACE s_axilite port=node_P_1 bundle=control
#pragma HLS INTERFACE s_axilite port=node_P_2 bundle=control
#pragma HLS INTERFACE s_axilite port=node_P_3 bundle=control

#pragma HLS INTERFACE s_axilite port=node_Q_0 bundle=control
#pragma HLS INTERFACE s_axilite port=node_Q_1 bundle=control
#pragma HLS INTERFACE s_axilite port=node_Q_2 bundle=control
#pragma HLS INTERFACE s_axilite port=node_Q_3 bundle=control

#pragma HLS INTERFACE s_axilite port=node_Vm_0 bundle=control
#pragma HLS INTERFACE s_axilite port=node_Vm_1 bundle=control
#pragma HLS INTERFACE s_axilite port=node_Vm_2 bundle=control
#pragma HLS INTERFACE s_axilite port=node_Vm_3 bundle=control

#pragma HLS INTERFACE s_axilite port=node_Va_0 bundle=control
#pragma HLS INTERFACE s_axilite port=node_Va_1 bundle=control
#pragma HLS INTERFACE s_axilite port=node_Va_2 bundle=control
#pragma HLS INTERFACE s_axilite port=node_Va_3 bundle=control

#pragma HLS INTERFACE s_axilite port=node_type_0 bundle=control
#pragma HLS INTERFACE s_axilite port=node_type_1 bundle=control
#pragma HLS INTERFACE s_axilite port=node_type_2 bundle=control
#pragma HLS INTERFACE s_axilite port=node_type_3 bundle=control

#pragma HLS INTERFACE s_axilite port=edge_off_0 bundle=control
#pragma HLS INTERFACE s_axilite port=edge_off_1 bundle=control
#pragma HLS INTERFACE s_axilite port=edge_off_2 bundle=control
#pragma HLS INTERFACE s_axilite port=edge_off_3 bundle=control

#pragma HLS INTERFACE s_axilite port=edge_src_0 bundle=control
#pragma HLS INTERFACE s_axilite port=edge_src_1 bundle=control
#pragma HLS INTERFACE s_axilite port=edge_src_2 bundle=control
#pragma HLS INTERFACE s_axilite port=edge_src_3 bundle=control

#pragma HLS INTERFACE s_axilite port=edge_G_0 bundle=control
#pragma HLS INTERFACE s_axilite port=edge_G_1 bundle=control
#pragma HLS INTERFACE s_axilite port=edge_G_2 bundle=control
#pragma HLS INTERFACE s_axilite port=edge_G_3 bundle=control

#pragma HLS INTERFACE s_axilite port=edge_B_0 bundle=control
#pragma HLS INTERFACE s_axilite port=edge_B_1 bundle=control
#pragma HLS INTERFACE s_axilite port=edge_B_2 bundle=control
#pragma HLS INTERFACE s_axilite port=edge_B_3 bundle=control

#pragma HLS INTERFACE s_axilite port=node_deltaP_0 bundle=control
#pragma HLS INTERFACE s_axilite port=node_deltaP_1 bundle=control
#pragma HLS INTERFACE s_axilite port=node_deltaP_2 bundle=control
#pragma HLS INTERFACE s_axilite port=node_deltaP_3 bundle=control

#pragma HLS INTERFACE s_axilite port=node_deltaQ_0 bundle=control
#pragma HLS INTERFACE s_axilite port=node_deltaQ_1 bundle=control
#pragma HLS INTERFACE s_axilite port=node_deltaQ_2 bundle=control
#pragma HLS INTERFACE s_axilite port=node_deltaQ_3 bundle=control

#pragma HLS INTERFACE s_axilite port=max_delta bundle=control

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
#pragma HLS INTERFACE s_axilite port=Vm_flag bundle=control
#pragma HLS INTERFACE s_axilite port=Va_flag bundle=control
#pragma HLS INTERFACE s_axilite port=deltaP_flag bundle=control
#pragma HLS INTERFACE s_axilite port=deltaQ_flag bundle=control
#pragma HLS INTERFACE s_axilite port=compute_flag bundle=control

#pragma HLS INTERFACE s_axilite port=return bundle=control



	static double l_P[4][MAX_NODE_NUM/4];
#pragma HLS ARRAY_PARTITION variable=l_P complete dim=1
#pragma HLS ARRAY_PARTITION variable=l_P dim=2 cyclic factor=8
	static double l_Q[4][MAX_NODE_NUM/4];
#pragma HLS ARRAY_PARTITION variable=l_Q complete dim=1
#pragma HLS ARRAY_PARTITION variable=l_Q dim=2 cyclic factor=8
	static int l_type[4][MAX_NODE_NUM/4];
#pragma HLS ARRAY_PARTITION variable=l_type complete dim=1
#pragma HLS ARRAY_PARTITION variable=l_type dim=2 cyclic factor=8

	static int l_off[4][2][MAX_NODE_NUM/4];
#pragma HLS ARRAY_PARTITION variable=l_off complete dim=1
#pragma HLS ARRAY_PARTITION variable=l_off complete dim=2
#pragma HLS ARRAY_PARTITION variable=l_off dim=3 cyclic factor=8
	static int l_src[4][MAX_EDGE_NUM/4];
#pragma HLS ARRAY_PARTITION variable=l_src complete dim=1
#pragma HLS ARRAY_PARTITION variable=l_src dim=2 cyclic factor=16
	static double l_G[4][MAX_EDGE_NUM/4];
#pragma HLS ARRAY_PARTITION variable=l_G complete dim=1
#pragma HLS ARRAY_PARTITION variable=l_G dim=2 cyclic factor=16
	static double l_B[4][MAX_EDGE_NUM/4];
#pragma HLS ARRAY_PARTITION variable=l_B complete dim=1
#pragma HLS ARRAY_PARTITION variable=l_B dim=2 cyclic factor=16


	static double l_Vm_dst[4][MAX_NODE_NUM/4];
#pragma HLS ARRAY_PARTITION variable=l_Vm_dst complete dim=1
#pragma HLS ARRAY_PARTITION variable=l_Vm_dst dim=2 cyclic factor=8

	static double l_Va_dst[4][MAX_NODE_NUM/4];
#pragma HLS ARRAY_PARTITION variable=l_Va_dst complete dim=1
#pragma HLS ARRAY_PARTITION variable=l_Va_dst dim=2 cyclic factor=8

	static double l_Vm_src[4][4][2][MAX_NODE_NUM/4];

#pragma HLS ARRAY_PARTITION variable=l_Vm_src complete dim=1
#pragma HLS ARRAY_PARTITION variable=l_Vm_src dim=4 cyclic factor=8
	static double l_Va_src[4][4][2][MAX_NODE_NUM/4];
#pragma HLS ARRAY_PARTITION variable=l_Va_src complete dim=1
#pragma HLS ARRAY_PARTITION variable=l_Va_src dim=4 cyclic factor=8

	double  l_deltaP[4][MAX_NODE_NUM/4];
#pragma HLS ARRAY_PARTITION variable=l_deltaP complete dim=1
#pragma HLS ARRAY_PARTITION variable=l_deltaP dim=2 cyclic factor=8
	double  l_deltaQ[4][MAX_NODE_NUM/4];
#pragma HLS ARRAY_PARTITION variable=l_deltaQ complete dim=1
#pragma HLS ARRAY_PARTITION variable=l_deltaQ dim=2 cyclic factor=8

	long  ll_deltaP[4][MAX_NODE_NUM/4];
	#pragma HLS ARRAY_PARTITION variable=ll_deltaP complete dim=1
	#pragma HLS ARRAY_PARTITION variable=ll_deltaP dim=2 cyclic factor=8
	long  ll_deltaQ[4][MAX_NODE_NUM/4];
	#pragma HLS ARRAY_PARTITION variable=ll_deltaQ complete dim=1
	#pragma HLS ARRAY_PARTITION variable=ll_deltaQ dim=2 cyclic factor=8


if(init_flag){
		init_bank_memcpy(node_P_0,l_P[0],
							node_Q_0,l_Q[0],
							node_type_0,l_type[0],
							edge_off_0,l_off[0][0],l_off[0][1],
							edge_src_0,l_src[0],
							edge_G_0,l_G[0],
							edge_B_0,l_B[0],
							begin_0,
							end_0
							);
		init_bank_memcpy(node_P_1,l_P[1],
							node_Q_1,l_Q[1],
							node_type_1,l_type[1],
							edge_off_1,l_off[1][0],l_off[1][1],
			                edge_src_1,l_src[1],
							edge_G_1,l_G[1],
							edge_B_1,l_B[1],
							begin_1,
							end_1
							);
		init_bank_memcpy(node_P_2,l_P[2],
							node_Q_2,l_Q[2],
							node_type_2,l_type[2],
							edge_off_2,l_off[2][0],l_off[2][1],
			                edge_src_2,l_src[2],
							edge_G_2,l_G[2],
							edge_B_2,l_B[2],
							begin_2,
							end_2
							);
		init_bank_memcpy(node_P_3,l_P[3],
							node_Q_3,l_Q[3],
							node_type_3,l_type[3],
							edge_off_3,l_off[3][0],l_off[3][1],
			                edge_src_3,l_src[3],
							edge_G_3,l_G[3],
							edge_B_3,l_B[3],
							begin_3,
							end_3
							);
}

if(Vm_flag){
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
}
	int begin[4],end[4];
	begin[0]=begin_0;end[0]=end_0;
	begin[1]=begin_1;end[1]=end_1;
	begin[2]=begin_2;end[2]=end_2;
	begin[3]=begin_3;end[3]=end_3;

	double l_max_deltaP[4]={0,0,0,0};
	double l_max_deltaQ[4]={0,0,0,0};
if(compute_flag){
	//for(int i=0;i<4;i++){
		//#pragma HLS UNROLL factor=4

		/*krnl_pf_pipe(
				l_P[i],
				l_Q[i],
				l_Vm_dst[i],
				l_Va_dst[i],
				l_Vm_src[0][i],l_Vm_src[1][i],l_Vm_src[2][i],l_Vm_src[3][i],
				l_Va_src[0][i],l_Va_src[1][i],l_Va_src[2][i],l_Va_src[3][i],
				begin[0],begin[1],begin[2],begin[3],
				l_type[i],
				l_off[i][0],l_off[i][1],l_off[i][2],
				l_src[i],
				l_G[i],
				l_B[i],
				l_deltaP[i],ll_deltaP[i],
				l_deltaQ[i],ll_deltaQ[i],
				&l_max_deltaP[i],
				&l_max_deltaQ[i],
				begin[i],
				end[i]
		);*/
	krnl_pf_pipe(
					l_P[0],
					l_Q[0],
					l_Vm_dst[0],
					l_Va_dst[0],
					l_Vm_src[0][0],l_Vm_src[1][0],l_Vm_src[2][0],l_Vm_src[3][0],
					l_Va_src[0][0],l_Va_src[1][0],l_Va_src[2][0],l_Va_src[3][0],
					begin[0],begin[1],begin[2],begin[3],
					l_type[0],
					l_off[0][0],l_off[0][1],
					l_src[0],
					l_G[0],
					l_B[0],
					l_deltaP[0],ll_deltaP[0],
					l_deltaQ[0],ll_deltaQ[0],
					&l_max_deltaP[0],
					&l_max_deltaQ[0],
					begin[0],
					end[0]
			);
	krnl_pf_pipe(
					l_P[1],
					l_Q[1],
					l_Vm_dst[1],
					l_Va_dst[1],
					l_Vm_src[0][1],l_Vm_src[1][1],l_Vm_src[2][1],l_Vm_src[3][1],
					l_Va_src[0][1],l_Va_src[1][1],l_Va_src[2][1],l_Va_src[3][1],
					begin[0],begin[1],begin[2],begin[3],
					l_type[1],
					l_off[1][0],l_off[1][1],
					l_src[1],
					l_G[1],
					l_B[1],
					l_deltaP[1],ll_deltaP[1],
					l_deltaQ[1],ll_deltaQ[1],
					&l_max_deltaP[1],
					&l_max_deltaQ[1],
					begin[1],
					end[1]
			);
	krnl_pf_pipe(
					l_P[2],
					l_Q[2],
					l_Vm_dst[2],
					l_Va_dst[2],
					l_Vm_src[0][2],l_Vm_src[1][2],l_Vm_src[2][2],l_Vm_src[3][2],
					l_Va_src[0][2],l_Va_src[1][2],l_Va_src[2][2],l_Va_src[3][2],
					begin[0],begin[1],begin[2],begin[3],
					l_type[2],
					l_off[2][0],l_off[2][1],
					l_src[2],
					l_G[2],
					l_B[2],
					l_deltaP[2],ll_deltaP[2],
					l_deltaQ[2],ll_deltaQ[2],
					&l_max_deltaP[2],
					&l_max_deltaQ[2],
					begin[2],
					end[2]
			);
	krnl_pf_pipe(
					l_P[3],
					l_Q[3],
					l_Vm_dst[3],
					l_Va_dst[3],
					l_Vm_src[0][3],l_Vm_src[1][3],l_Vm_src[2][3],l_Vm_src[3][3],
					l_Va_src[0][3],l_Va_src[1][3],l_Va_src[2][3],l_Va_src[3][3],
					begin[0],begin[1],begin[2],begin[3],
					l_type[3],
					l_off[3][0],l_off[3][1],
					l_src[3],
					l_G[3],
					l_B[3],
					l_deltaP[3],ll_deltaP[3],
					l_deltaQ[3],ll_deltaQ[3],
					&l_max_deltaP[3],
					&l_max_deltaQ[3],
					begin[3],
					end[3]
			);

	//}
}

	double max_delta_tmp[8];
	max_delta_tmp[0]=fmax(fmax(l_max_deltaP[0],l_max_deltaP[1]),fmax(l_max_deltaP[2],l_max_deltaP[3]));
	max_delta_tmp[1]=fmax(fmax(l_max_deltaQ[0],l_max_deltaQ[1]),fmax(l_max_deltaQ[2],l_max_deltaQ[3]));
	max_delta_tmp[2]=0;
	max_delta_tmp[3]=0;
	max_delta_tmp[4]=0;
	max_delta_tmp[5]=0;
	max_delta_tmp[6]=0;
	max_delta_tmp[7]=0;

	TYPE tmp;
	for(int j=0;j<8;j++){
		#pragma HLS unroll factor=8

		double data=max_delta_tmp[j];
		union Int64_or_Double data_tmp;
		data_tmp.data_d=data;
		long data_l;
		data_l=data_tmp.data_i;
		tmp.range(64*j+63,64*j)=data_l;
	}
	max_delta[0]=tmp;


if(deltaP_flag){
	deltaP_bank_memcpy(node_deltaP_0,l_deltaP[0],begin_0,end_0);
	deltaP_bank_memcpy(node_deltaP_1,l_deltaP[1],begin_1,end_1);
	deltaP_bank_memcpy(node_deltaP_2,l_deltaP[2],begin_2,end_2);
	deltaP_bank_memcpy(node_deltaP_3,l_deltaP[3],begin_3,end_3);
}
if(deltaQ_flag){
	deltaQ_bank_memcpy(node_deltaQ_0,l_deltaQ[0],begin_0,end_0);
	deltaQ_bank_memcpy(node_deltaQ_1,l_deltaQ[1],begin_1,end_1);
	deltaQ_bank_memcpy(node_deltaQ_2,l_deltaQ[2],begin_2,end_2);
	deltaQ_bank_memcpy(node_deltaQ_3,l_deltaQ[3],begin_3,end_3);
}

}
}

