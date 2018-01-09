//http://www.holoborodko.com/pavel/numerical-methods/numerical-integration/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "gauss_legendre.h"

#ifndef PI
	#define PI 3.1415926535897932384626433832795028841971693993751
#endif
#ifndef FABS
	#define FABS(a) ((a)>=0?(a):-(a))
#endif

/*
double f(double x, void* data)
{
	return sin(x);
}
*/
int main(int argc, char* argv[])
{
/*
	double approx;		
	
	double exact = 2.0; 

	double error;       

	int i;

	printf("Numerical Approximation of int(sin(x), x=0..Pi) by Gauss-Legendre Quadrature:\n");
	for (i=2;i<=128;i++)
	{
		approx = gauss_legendre(i,f,NULL,0,PI);
		error = approx-exact;
		printf("n = %4d: error = %.15g\n",i,FABS(error));
	}
	*/
	FILE *f;
	f=fopen("bessel_test.txt","w");
	fprintf(f,"besselJ(0,x+yi)/n");
	fprintf(f,"real,     imag,    Rvalue,                         Ivalue/n");
	for(int i=0;i<10000;i++){
		fprintf(f,"%g,%g,%.15g,%.15g/n",double(i/100)*0.1,double(i%100)*0.1,creal(besselJ(0,double(i/100)*0.1+I*double(i%100)*0.1)),cimag(besselJ(0,0.1*double(i/100)+I*0.1*double(i%100)));
	}
	fprintf(f,"/n/n besselJ(1,x+yi)/n");
	for(int i=0;i<10000;i++){
		fprintf(f,"%g,%g,%.15g,%.15g/n",double(i/100)*0.1,double(i%100)*0.1,creal(besselJ(1,double(i/100)*0.1+I*double(i%100)*0.1)),cimag(besselJ(1,0.1*double(i/100)+I*0.1*double(i%100)));
	}
	fclose(f);
}

