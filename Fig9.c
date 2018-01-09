#include <stdio.h>
#include <stdlib.h>
#include "math.h"
#include "mkj.h"
#ifndef PI
	#define PI 3.1415926535897932384626433832795028841971693993751
#endif

void main(){
	FILE *f,*g;
	//f=fopen("Gxx_Fig91.txt","w");
	double complex value1,value2;
	f=fopen("Gxx_Fig9.txt","w");
	g=fopen("Gzx_FIg9.txt","w");
	double depth[2]={10.0,0.0};
	Sys *S; S=createSys((parm){.N=3,.d=depth, .k0=2*PI/633.0});
	S->Eps[0][0]=1.0;S->Eps[0][1]=1.0;S->Eps[0][2]=1.0;
	S->Eps[1][0]=1.0;S->Eps[1][1]=-2.475+0.5054*I;S->Eps[1][2]=-2.475+0.5054*I;
	S->Eps[2][0]=1.0;S->Eps[2][1]=1.0;S->Eps[2][2]=1.0;
	double rho=2.0;
	double phi=0.0;
	for(int i=0;i<100;i++){
		//value1=0.25/PI*(integral_imag(S,rho,phi,-1000.0+20.0*i,-1000.0+20.0*i, fpxx)+integral_imag(S,rho,phi,-1000.0+20.0*i,-1000.0+20.0*i, fsxx));
		value1=0.25/PI*(integral_imag(S,rho,phi,-10.0+0.2*i,1.5, fsxx)+integral_imag(S,rho,phi,-10.0+0.2*i,1.5, fpxx));
		//value3=0.25/PI*(integral_imag(S,rho,phi+PI,750.0, -1000.0+20.0*i,fsxx)+integral_imag(S,rho,phi+PI,750.0, -1000.0+20.0*i,fpxx));
		
		value2=0.25/PI*(integral_imag(S,rho,phi,-10.0+0.2*i,1.5,fszx)+integral_imag(S,rho,phi,-10.0+0.2*i,1.5, fpzx));
		//value2=0.25/PI*(integral_imag(S,rho,phi,-1000.0+20.0*i,-1000.0+20.0*i, fpzz));
		//value4=0.25/PI*(integral_imag(S,rho,phi+PI,750.0,-1000.0+20.0*i,fsxz)+integral_imag(S,rho,phi+PI,750.0, -1000.0+20.0*i,fpxz));
		//fprintf(f,"%d,%.15f,%.15f,%.15f,%.15f\n",-1000+20*i,creal(value1),cimag(value1),creal(value3),cimag(value3));
		//fprintf(g,"%d,%.15f,%.15f,%.15f,%.15f\n",-1000+20*i,creal(value2),cimag(value2),creal(value4),cimag(value4));
		fprintf(f,"%f,%.15f\n",-10+0.2*i,cabs(value1));
		fprintf(g,"%f,%.15f\n",-10+0.2*i,cabs(value2));
	}
	/*
	ItgParm data;
	double depth[3]={5.0,0.0,-5.0};
	Sys *S; S=createSys((parm){.N=4,.d=depth, .k0=2*PI/633.0});
	S->Eps[0][0]=1.0;S->Eps[0][1]=1.0;S->Eps[0][2]=1.0;
	S->Eps[1][0]=1.0;S->Eps[1][1]=2.0;S->Eps[1][2]=2.0;
	S->Eps[2][0]=1.0;S->Eps[2][1]=10.0;S->Eps[2][2]=10.0;
	S->Eps[3][0]=1.0;S->Eps[3][1]=1.0;S->Eps[3][2]=1.0;
	//double rho=0.0;
	double rho=6.33;
	double phi=PI*0.25;
	//double phi=0.0;
	double complex value1,value2,value3,value4;
	data=(ItgParm){.func_num=1, .rho=sqrt(12.0), .phi=PI*0.25, .z=5.0, .zs=3.0, .S=S, .f=fsxx};
	value1=fpxx(1000.0,&data);
	printf("%f,%f",creal(value1),cimag(value1));*/
	/*
	for(int i=0;i<100;i++){
		//value1=0.25/PI*(integral_imag(S,rho,phi,-1000.0+20.0*i,-1000.0+20.0*i, fpxx)+integral_imag(S,rho,phi,-1000.0+20.0*i,-1000.0+20.0*i, fsxx));
		value1=0.25/PI*(integral_imag(S,rho,phi,-10.0+0.2*i,7.5, fsxx)+integral_imag(S,rho,phi,-10.0+0.2*i,7.5, fpxx));
		//value3=0.25/PI*(integral_imag(S,rho,phi+PI,750.0, -1000.0+20.0*i,fsxx)+integral_imag(S,rho,phi+PI,750.0, -1000.0+20.0*i,fpxx));
		
		value2=0.25/PI*(integral_imag(S,rho,phi,-10.0+0.2*i,7.5,fszx)+integral_imag(S,rho,phi,-10.0+0.2*i,7.5, fpzx));
		//value2=0.25/PI*(integral_imag(S,rho,phi,-1000.0+20.0*i,-1000.0+20.0*i, fpzz));
		//value4=0.25/PI*(integral_imag(S,rho,phi+PI,750.0,-1000.0+20.0*i,fsxz)+integral_imag(S,rho,phi+PI,750.0, -1000.0+20.0*i,fpxz));
		//fprintf(f,"%d,%.15f,%.15f,%.15f,%.15f\n",-1000+20*i,creal(value1),cimag(value1),creal(value3),cimag(value3));
		//fprintf(g,"%d,%.15f,%.15f,%.15f,%.15f\n",-1000+20*i,creal(value2),cimag(value2),creal(value4),cimag(value4));
		fprintf(f,"%f,%.15f\n",-10+0.2*i,cabs(value1));
		fprintf(g,"%f,%.15f\n",-10+0.2*i,cabs(value2));
	}*/
	fclose(f);
	fclose(g);
}