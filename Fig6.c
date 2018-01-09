#include <stdio.h>
#include <stdlib.h>
#include "math.h"
#include "mkj.h"


inline double MPI(){ return 3.1415926535897932384626433832795028841971693993751;}


void main(){
	FILE *f;
	f=fopen("Gzx_test1.txt","w");
	//f=fopen("Gzx_test3.txt","w");

	//double depth[5]={1000.0,500.0,0.0,-500.0,-1000.0};
	//Sys *S; S=createSys((parm){.N=6, .d=depth,.k0=2*MPI()/633.0});
	Sys *S; S=createSys((parm){.N=1, .d=NULL,.k0=2*MPI()/633.0});
	S->Eps[0][0]=1.0;S->Eps[0][1]=1.0;S->Eps[0][2]=1.0;
	//S->Eps[1][0]=1;S->Eps[1][1]=1;S->Eps[1][2]=1;S->Eps[2][0]=1;S->Eps[2][1]=1;S->Eps[2][2]=1;S->Eps[3][0]=1;S->Eps[3][1]=1;S->Eps[3][2]=1;
	//S->Eps[4][0]=1;S->Eps[4][1]=1;S->Eps[4][2]=1;S->Eps[5][0]=1;S->Eps[5][1]=1;S->Eps[5][2]=1;
	double complex k1,k2;
	//double rho=0.0;
	double rho=633.0;
	double phi=0.0;
	double complex krho;
	//k1=0.25*I/MPI()*integral_imag(S,rho,phi,-1000.0,750.0,fsxx);
	//printf("%f,%f\n",creal(k1),cimag(k1));
	
	for(int i=0;i<1;i++){
		//S->k_0=2*MPI()/(500.0+5.0*i);

		//k1=0.25*I/MPI()*(integral_imag(S,rho,phi,0.0,0.0,fsxx)+integral_imag(S,rho,phi,0.0,0.0,fpxx));
		//k2=0.25*I/MPI()*integral_imag(S,rho,phi,0.0,0.0,fpzz);
		k1=0.25*I/MPI()*(integral_imag(S,rho,phi,-1000.0+20.0*i,750.0,fpxx)+integral_imag(S,rho,phi,-1000.0+20.0*i,750.0,fsxx));
		k2=0.25*I/MPI()*integral_imag(S,rho,phi,-1000.0+20.0*i,750.0,fpzx);
		//fprintf(f,"%d,%.20f,%.20f\n",500+5*i,cimag(k1),cimag(k2));
		fprintf(f,"%d,%.20f,%.20f\n",-1000+20*i,cabs(k1),cabs(k2));
	}
	
	fclose(f);

}