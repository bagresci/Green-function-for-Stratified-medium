#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mkj.h"
#define PI 3.1415926535897932384626433832795028841971693993751
/*double f(double x, void* data)
{
	return sin(x);
}*/

int main()
{
/*
	double complex z=besselJ(0,1+I);
	printf("please%.15g\n",creal(z));
	*/
	FILE *f;
	f=fopen("bessel_test.txt","w");
//	fprintf(f,"besselJ(0,x+yi)/n");
//	fprintf(f,"real,     imag,    Rvalue,                         Ivalue/n");
	double complex z;
/*	for(int i=0;i<10000;i++){
		z=besselJ(0,(i/100)*0.1+I*(i%100)*0.1);
		fprintf(f,"%g,%g,%.15f,%.15f\n",(i/100)*0.1,(i%100)*0.1,creal(z),cimag(z));
	}
	//fprintf(f,"/n/n besselJ(1,x+yi)/n");
	for(int i=0;i<10000;i++){
		z=besselJ(1,(i/100)*0.1+I*(i%100)*0.1);
		fprintf(f,"%g,%g,%.15f,%.15f\n",(i/100)*0.1,(i%100)*0.1,creal(z),cimag(z));
	}
	*/
	for(int j=0;j<100;j++){
		z=besselJ(0,j*0.1);
		fprintf(f,"%g,%d,%.15f,%.15f\n",0.1*j,0,creal(z),cimag(z));
	}
	for(int j=0;j<100;j++){
		z=besselJ(0,I*0.1*j);
		fprintf(f,"%d,%g,%.15f,%.15f\n",0,j*0.1,creal(z),cimag(z));
	}
	for(int j=0;j<100;j++){
		z=besselJ(0,(1+I)*0.1*j);
		fprintf(f,"%g,%g,%.15f,%.15f\n",0.1*j,j*0.1,creal(z),cimag(z));
	}
	fclose(f);
}

