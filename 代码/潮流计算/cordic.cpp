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

//angles = atan(2.^-(0:27));
static double angles[36] =  {
	0.78539816339745,   0.46364760900081,   0.24497866312686,   0.12435499454676,
	0.06241880999596,   0.03123983343027,   0.01562372862048,   0.00781234106010,
	0.00390623013197,   0.00195312251648,   0.00097656218956,   0.00048828121119,
	0.00024414062015,   0.00012207031189,   0.00006103515617,   0.00003051757812,
	0.00001525878906,   0.00000762939453,   0.00000381469727,   0.00000190734863,
	0.00000095367432,   0.00000047683716,  	0.00000023841858,   0.00000011920929,
	0.00000005960464,   0.00000002980232,   0.00000001490116,   0.00000000745058,

	0.00000000372529,	0.00000000186265,	0.00000000093132,	0.00000000046566,
	0.00000000023283,	0.00000000011642,	0.00000000005821,	0.00000000002910

};

static double Kvalues[24] = {
	0.70710678118655,   0.63245553203368,   0.61357199107790,   0.60883391251775,
	0.60764825625617,   0.60735177014130,   0.60727764409353,   0.60725911229889,
	0.60725447933256,   0.60725332108988,   0.60725303152913,   0.60725295913894,
	0.60725294104140,   0.60725293651701,   0.60725293538591,   0.60725293510314,
	0.60725293503245,   0.60725293501477,   0.60725293501035,   0.60725293500925,
	0.60725293500897,   0.60725293500890,  	0.60725293500889,   0.60725293500888
};

void cordic_code_generate(){
	for(int i=0;i<32;i++){
		printf("//Loop %d\n",i);
		printf("double factor_%d,beta_%d;\n",i,i+1);

		double left_i;
		if(i<31) left_i=1.0/(1<<i);
		else left_i=1.0/(1<<(i-4))/16;
		printf("if(beta_%d<0) {factor_%d=-%.14lf;beta_%d = beta_%d + %.14lf;}\n",
				i,i,left_i,i+1,i,angles[i]);
		printf("else {factor_%d=%.14lf;beta_%d = beta_%d - %.14lf;}\n",
				i,left_i,i+1,i,angles[i]);

		printf("double x_cos_%d = x_cos_%d - factor_%d * y_sin_%d;\n",i+1,i,i,i);
		printf("double y_sin_%d = factor_%d * x_cos_%d + y_sin_%d;\n",i+1,i,i,i);

		printf("\n");
	}
}

void cordic(double angle_para,double *cosine, double *sine){
	double beta;
	if(angle_para>=PI) beta=angle_para-2*PI;
	else if(angle_para<-PI) beta=angle_para+2*PI;
	else beta=angle_para;

	double beta_0;
	bool flag;
	if(beta>=PI/2) {beta_0=beta-PI;flag=true;}
	else if (beta<-PI/2) {beta_0=beta+PI;flag=true;}
	else {beta_0=beta;flag=false;}

	double x_cos_0=0.60725293500888,y_sin_0=0;//initial

	//Loop 0
	double factor_0,beta_1;
	if(beta_0<0) {factor_0=-1.00000000000000;beta_1 = beta_0 + 0.78539816339745;}
	else {factor_0=1.00000000000000;beta_1 = beta_0 - 0.78539816339745;}
	double x_cos_1 = x_cos_0 - factor_0 * y_sin_0;
	double y_sin_1 = factor_0 * x_cos_0 + y_sin_0;

	//Loop 1
	double factor_1,beta_2;
	if(beta_1<0) {factor_1=-0.50000000000000;beta_2 = beta_1 + 0.46364760900081;}
	else {factor_1=0.50000000000000;beta_2 = beta_1 - 0.46364760900081;}
	double x_cos_2 = x_cos_1 - factor_1 * y_sin_1;
	double y_sin_2 = factor_1 * x_cos_1 + y_sin_1;

	//Loop 2
	double factor_2,beta_3;
	if(beta_2<0) {factor_2=-0.25000000000000;beta_3 = beta_2 + 0.24497866312686;}
	else {factor_2=0.25000000000000;beta_3 = beta_2 - 0.24497866312686;}
	double x_cos_3 = x_cos_2 - factor_2 * y_sin_2;
	double y_sin_3 = factor_2 * x_cos_2 + y_sin_2;

	//Loop 3
	double factor_3,beta_4;
	if(beta_3<0) {factor_3=-0.12500000000000;beta_4 = beta_3 + 0.12435499454676;}
	else {factor_3=0.12500000000000;beta_4 = beta_3 - 0.12435499454676;}
	double x_cos_4 = x_cos_3 - factor_3 * y_sin_3;
	double y_sin_4 = factor_3 * x_cos_3 + y_sin_3;

	//Loop 4
	double factor_4,beta_5;
	if(beta_4<0) {factor_4=-0.06250000000000;beta_5 = beta_4 + 0.06241880999596;}
	else {factor_4=0.06250000000000;beta_5 = beta_4 - 0.06241880999596;}
	double x_cos_5 = x_cos_4 - factor_4 * y_sin_4;
	double y_sin_5 = factor_4 * x_cos_4 + y_sin_4;

	//Loop 5
	double factor_5,beta_6;
	if(beta_5<0) {factor_5=-0.03125000000000;beta_6 = beta_5 + 0.03123983343027;}
	else {factor_5=0.03125000000000;beta_6 = beta_5 - 0.03123983343027;}
	double x_cos_6 = x_cos_5 - factor_5 * y_sin_5;
	double y_sin_6 = factor_5 * x_cos_5 + y_sin_5;

	//Loop 6
	double factor_6,beta_7;
	if(beta_6<0) {factor_6=-0.01562500000000;beta_7 = beta_6 + 0.01562372862048;}
	else {factor_6=0.01562500000000;beta_7 = beta_6 - 0.01562372862048;}
	double x_cos_7 = x_cos_6 - factor_6 * y_sin_6;
	double y_sin_7 = factor_6 * x_cos_6 + y_sin_6;

	//Loop 7
	double factor_7,beta_8;
	if(beta_7<0) {factor_7=-0.00781250000000;beta_8 = beta_7 + 0.00781234106010;}
	else {factor_7=0.00781250000000;beta_8 = beta_7 - 0.00781234106010;}
	double x_cos_8 = x_cos_7 - factor_7 * y_sin_7;
	double y_sin_8 = factor_7 * x_cos_7 + y_sin_7;

	//Loop 8
	double factor_8,beta_9;
	if(beta_8<0) {factor_8=-0.00390625000000;beta_9 = beta_8 + 0.00390623013197;}
	else {factor_8=0.00390625000000;beta_9 = beta_8 - 0.00390623013197;}
	double x_cos_9 = x_cos_8 - factor_8 * y_sin_8;
	double y_sin_9 = factor_8 * x_cos_8 + y_sin_8;

	//Loop 9
	double factor_9,beta_10;
	if(beta_9<0) {factor_9=-0.00195312500000;beta_10 = beta_9 + 0.00195312251648;}
	else {factor_9=0.00195312500000;beta_10 = beta_9 - 0.00195312251648;}
	double x_cos_10 = x_cos_9 - factor_9 * y_sin_9;
	double y_sin_10 = factor_9 * x_cos_9 + y_sin_9;

	//Loop 10
	double factor_10,beta_11;
	if(beta_10<0) {factor_10=-0.00097656250000;beta_11 = beta_10 + 0.00097656218956;}
	else {factor_10=0.00097656250000;beta_11 = beta_10 - 0.00097656218956;}
	double x_cos_11 = x_cos_10 - factor_10 * y_sin_10;
	double y_sin_11 = factor_10 * x_cos_10 + y_sin_10;

	//Loop 11
	double factor_11,beta_12;
	if(beta_11<0) {factor_11=-0.00048828125000;beta_12 = beta_11 + 0.00048828121119;}
	else {factor_11=0.00048828125000;beta_12 = beta_11 - 0.00048828121119;}
	double x_cos_12 = x_cos_11 - factor_11 * y_sin_11;
	double y_sin_12 = factor_11 * x_cos_11 + y_sin_11;

	//Loop 12
	double factor_12,beta_13;
	if(beta_12<0) {factor_12=-0.00024414062500;beta_13 = beta_12 + 0.00024414062015;}
	else {factor_12=0.00024414062500;beta_13 = beta_12 - 0.00024414062015;}
	double x_cos_13 = x_cos_12 - factor_12 * y_sin_12;
	double y_sin_13 = factor_12 * x_cos_12 + y_sin_12;

	//Loop 13
	double factor_13,beta_14;
	if(beta_13<0) {factor_13=-0.00012207031250;beta_14 = beta_13 + 0.00012207031189;}
	else {factor_13=0.00012207031250;beta_14 = beta_13 - 0.00012207031189;}
	double x_cos_14 = x_cos_13 - factor_13 * y_sin_13;
	double y_sin_14 = factor_13 * x_cos_13 + y_sin_13;

	//Loop 14
	double factor_14,beta_15;
	if(beta_14<0) {factor_14=-0.00006103515625;beta_15 = beta_14 + 0.00006103515617;}
	else {factor_14=0.00006103515625;beta_15 = beta_14 - 0.00006103515617;}
	double x_cos_15 = x_cos_14 - factor_14 * y_sin_14;
	double y_sin_15 = factor_14 * x_cos_14 + y_sin_14;

	//Loop 15
	double factor_15,beta_16;
	if(beta_15<0) {factor_15=-0.00003051757812;beta_16 = beta_15 + 0.00003051757812;}
	else {factor_15=0.00003051757812;beta_16 = beta_15 - 0.00003051757812;}
	double x_cos_16 = x_cos_15 - factor_15 * y_sin_15;
	double y_sin_16 = factor_15 * x_cos_15 + y_sin_15;

	//Loop 16
	double factor_16,beta_17;
	if(beta_16<0) {factor_16=-0.00001525878906;beta_17 = beta_16 + 0.00001525878906;}
	else {factor_16=0.00001525878906;beta_17 = beta_16 - 0.00001525878906;}
	double x_cos_17 = x_cos_16 - factor_16 * y_sin_16;
	double y_sin_17 = factor_16 * x_cos_16 + y_sin_16;

	//Loop 17
	double factor_17,beta_18;
	if(beta_17<0) {factor_17=-0.00000762939453;beta_18 = beta_17 + 0.00000762939453;}
	else {factor_17=0.00000762939453;beta_18 = beta_17 - 0.00000762939453;}
	double x_cos_18 = x_cos_17 - factor_17 * y_sin_17;
	double y_sin_18 = factor_17 * x_cos_17 + y_sin_17;

	//Loop 18
	double factor_18,beta_19;
	if(beta_18<0) {factor_18=-0.00000381469727;beta_19 = beta_18 + 0.00000381469727;}
	else {factor_18=0.00000381469727;beta_19 = beta_18 - 0.00000381469727;}
	double x_cos_19 = x_cos_18 - factor_18 * y_sin_18;
	double y_sin_19 = factor_18 * x_cos_18 + y_sin_18;

	//Loop 19
	double factor_19,beta_20;
	if(beta_19<0) {factor_19=-0.00000190734863;beta_20 = beta_19 + 0.00000190734863;}
	else {factor_19=0.00000190734863;beta_20 = beta_19 - 0.00000190734863;}
	double x_cos_20 = x_cos_19 - factor_19 * y_sin_19;
	double y_sin_20 = factor_19 * x_cos_19 + y_sin_19;

	//Loop 20
	double factor_20,beta_21;
	if(beta_20<0) {factor_20=-0.00000095367432;beta_21 = beta_20 + 0.00000095367432;}
	else {factor_20=0.00000095367432;beta_21 = beta_20 - 0.00000095367432;}
	double x_cos_21 = x_cos_20 - factor_20 * y_sin_20;
	double y_sin_21 = factor_20 * x_cos_20 + y_sin_20;

	//Loop 21
	double factor_21,beta_22;
	if(beta_21<0) {factor_21=-0.00000047683716;beta_22 = beta_21 + 0.00000047683716;}
	else {factor_21=0.00000047683716;beta_22 = beta_21 - 0.00000047683716;}
	double x_cos_22 = x_cos_21 - factor_21 * y_sin_21;
	double y_sin_22 = factor_21 * x_cos_21 + y_sin_21;

	//Loop 22
	double factor_22,beta_23;
	if(beta_22<0) {factor_22=-0.00000023841858;beta_23 = beta_22 + 0.00000023841858;}
	else {factor_22=0.00000023841858;beta_23 = beta_22 - 0.00000023841858;}
	double x_cos_23 = x_cos_22 - factor_22 * y_sin_22;
	double y_sin_23 = factor_22 * x_cos_22 + y_sin_22;

	//Loop 23
	double factor_23,beta_24;
	if(beta_23<0) {factor_23=-0.00000011920929;beta_24 = beta_23 + 0.00000011920929;}
	else {factor_23=0.00000011920929;beta_24 = beta_23 - 0.00000011920929;}
	double x_cos_24 = x_cos_23 - factor_23 * y_sin_23;
	double y_sin_24 = factor_23 * x_cos_23 + y_sin_23;

	//Loop 24
	double factor_24,beta_25;
	if(beta_24<0) {factor_24=-0.00000005960464;beta_25 = beta_24 + 0.00000005960464;}
	else {factor_24=0.00000005960464;beta_25 = beta_24 - 0.00000005960464;}
	double x_cos_25 = x_cos_24 - factor_24 * y_sin_24;
	double y_sin_25 = factor_24 * x_cos_24 + y_sin_24;

	//Loop 25
	double factor_25,beta_26;
	if(beta_25<0) {factor_25=-0.00000002980232;beta_26 = beta_25 + 0.00000002980232;}
	else {factor_25=0.00000002980232;beta_26 = beta_25 - 0.00000002980232;}
	double x_cos_26 = x_cos_25 - factor_25 * y_sin_25;
	double y_sin_26 = factor_25 * x_cos_25 + y_sin_25;

	//Loop 26
	double factor_26,beta_27;
	if(beta_26<0) {factor_26=-0.00000001490116;beta_27 = beta_26 + 0.00000001490116;}
	else {factor_26=0.00000001490116;beta_27 = beta_26 - 0.00000001490116;}
	double x_cos_27 = x_cos_26 - factor_26 * y_sin_26;
	double y_sin_27 = factor_26 * x_cos_26 + y_sin_26;

	//Loop 27
	double factor_27,beta_28;
	if(beta_27<0) {factor_27=-0.00000000745058;beta_28 = beta_27 + 0.00000000745058;}
	else {factor_27=0.00000000745058;beta_28 = beta_27 - 0.00000000745058;}
	double x_cos_28 = x_cos_27 - factor_27 * y_sin_27;
	double y_sin_28 = factor_27 * x_cos_27 + y_sin_27;

/*	//Loop 28
	double factor_28,beta_29;
	if(beta_28<0) {factor_28=-0.00000000372529;beta_29 = beta_28 + 0.00000000372529;}
	else {factor_28=0.00000000372529;beta_29 = beta_28 - 0.00000000372529;}
	double x_cos_29 = x_cos_28 - factor_28 * y_sin_28;
	double y_sin_29 = factor_28 * x_cos_28 + y_sin_28;

	//Loop 29
	double factor_29,beta_30;
	if(beta_29<0) {factor_29=-0.00000000186265;beta_30 = beta_29 + 0.00000000186265;}
	else {factor_29=0.00000000186265;beta_30 = beta_29 - 0.00000000186265;}
	double x_cos_30 = x_cos_29 - factor_29 * y_sin_29;
	double y_sin_30 = factor_29 * x_cos_29 + y_sin_29;

	//Loop 30
	double factor_30,beta_31;
	if(beta_30<0) {factor_30=-0.00000000093132;beta_31 = beta_30 + 0.00000000093132;}
	else {factor_30=0.00000000093132;beta_31 = beta_30 - 0.00000000093132;}
	double x_cos_31 = x_cos_30 - factor_30 * y_sin_30;
	double y_sin_31 = factor_30 * x_cos_30 + y_sin_30;

	//Loop 31
	double factor_31,beta_32;
	if(beta_31<0) {factor_31=-0.00000000046566;beta_32 = beta_31 + 0.00000000046566;}
	else {factor_31=0.00000000046566;beta_32 = beta_31 - 0.00000000046566;}
	double x_cos_32 = x_cos_31 - factor_31 * y_sin_31;
	double y_sin_32 = factor_31 * x_cos_31 + y_sin_31;
*/


	double c,s;
	if(flag){
		c=-x_cos_28;
		s=-y_sin_28;
	}
	else{
		c=x_cos_28;
		s=y_sin_28;
	}
	*cosine=c;
	*sine=s;

}


/*
void cordic(double angle_para,double *cosine, double *sine){
	//assert(angle_para<=2*PI&&angle_para>=-2*PI);

	double beta;
	int flag;

	if(angle_para>=PI) beta=angle_para-2*PI;
	else if(angle_para<-PI) beta=angle_para+2*PI;
	else beta=angle_para;

	if(beta>=PI/2) {beta-=PI;flag=1;}
	else if (beta<-PI/2) {beta+=PI;flag=1;}
	else {flag=0;}

	double c=Kvalues[23],s=0;//initial

	double poweroftwo=1;
	double angle=angles[0];
	for(int i=0;i<28;i++){//<28-1
		double sigma;
		if(beta<0) sigma=-1;
		else sigma=1;

		double factor = sigma * poweroftwo;
		double temp_c=c;
		double temp_s=s;
		c=temp_c-factor*temp_s;
		s=factor*temp_c+temp_s;

		beta = beta - sigma * angle;// % update the remaining angle
		poweroftwo = poweroftwo / 2;
		angle = angles[i+1];//if i+2>28 then angle=angle/2
	}

	if(flag){
		c=-c;
		s=-s;
	}
	*cosine=c;
	*sine=s;
}
*/

/*void cordic(double angle_para,double *cosine, double *sine){
	const double tangent[ ] ={1.0, 1/2.0, 1/4.0, 1/8.0, 1/16.0, 1/32.0, 1/64.0, 1/128.0, 1/512.0};
	const double angle [ ] = { 45.0, 26.6, 14, 7.1, 3.6, 1.8, 0.9, 0.4, 0.2, 0.1};
	int i, signal;
	double x_cos, y_sin, x_temp, y_temp, z, z_next;
	x_cos =0; y_sin=0; z=angle_para*180/PI; z_next=0;
	x_temp =0.6073; y_temp=0;
	signal = 1;
	for (i=0; i<9; i++){
		x_cos=x_temp - signal * y_temp * tangent [ i ];
		y_sin=y_temp+ signal * x_temp * tangent [ i ];
		z_next = z - signal * angle [ i ];

		x_temp = x_cos;
		y_temp = y_sin;

		z = z_next;
		if(z_next>0) signal = +1;
		else	signal = -1;
	}
	*cosine=x_cos;
	*sine=y_sin;
}*/

