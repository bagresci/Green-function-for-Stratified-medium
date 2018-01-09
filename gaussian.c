#include <stdlib.h>
#include "mkj.h"
#define PI 3.1415926535897932384626433832795028841971693993751
double complex gammafunc(double t,void * data){
	double complex v;
	v=*(double complex *)data;
	return cexp(1/(t*t-1)+1.5+I*sin(PI*t))*cpow(1/(t*t-1)+1.5+I*sin(PI*t),v)*(-2*t/(t*t-1)/(t*t-1)+I*PI*cos(PI*t))*cexp(PI*I*(v+1))/(cexp(I*2*PI*v)-1);
}//branch cut:negative real axis
void main(){
	FILE *f;
	f=fopen("gamma_gaussian.txt","w");
	double complex z;
	double complex vv;
	void * data=&vv;
	for(int i=0;i<9;i++){
		vv=1.1+i*0.1;
		z=gauss_legendre_complex(256,gammafunc,data,-1,1);
		fprintf(f,"%f,%.15f\n",creal(vv),creal(z));
	}
	for(int i=0;i<9;i++){
		vv=2.1+i*0.1;
		z=gauss_legendre_complex(256,gammafunc,data,-1,1);
		fprintf(f,"%f,%.15f\n",creal(vv),creal(z));
	}
	for(int i=0;i<9;i++){
		vv=3.1+i*0.1;
		z=gauss_legendre_complex(256,gammafunc,data,-1,1);
		fprintf(f,"%f,%.15f\n",creal(vv),creal(z));
	}
	for(int i=0;i<9;i++){
		vv=4.1+i*0.1;
		z=gauss_legendre_complex(256,gammafunc,data,-1,1);
		fprintf(f,"%f,%.15f\n",creal(vv),creal(z));
	}		
	fclose(f);
	vv=0.5;
	z=gauss_legendre_complex(256,gammafunc,data,-1,1);
	printf("%.15f,%.15f\n",creal(z),cimag(z));
}
