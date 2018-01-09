#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mkj.h"

inline double MPI(){
	return  3.1415926535897932384626433832795028841971693993751;
}

void main(){
	FILE *f;
	f=fopen("fxx.txt","w");
	Sys *S; S=createSys((parm){.N=1,.d=NULL, .k0=2*MPI()/500.0});
	S->Eps[0][0]=1;S->Eps[0][1]=1;S->Eps[0][2]=1;
	double complex value;
	ItgParm P;
	P.rho=10.0;P.phi=MPI()/2.0;P.z=1000.0;P.zs=0.0;P.S=S;P.f=fpxx;P.func_num=2;
//	void * data; data = (void *) &P;
	double kmaj=S->k_0;


	for(int i=0;i<1000;i++){
		value=fsxx(-i*I+2.0*kmaj,&P);
		fprintf(f,"%.15f,%.15f\n",creal(value),cimag(value));
	}
	fclose(f);

}