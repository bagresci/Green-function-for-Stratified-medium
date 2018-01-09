#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <string.h>
#include "gauss_legendre.h"
//#define PI 3.1415926535897932384626433832795028841971693993751



typedef struct {
	int N; double * d; double k0;
} parm;

typedef struct {
	int N; double * d; double k_0;
	//double complex (*Eps)[3];
	double complex ** Eps;
} Sys;


typedef struct {
	int func_num;
	double rho, phi, z, zs;
	Sys * S;
	double complex (*f) (double complex, void *);
} ItgParm;

typedef struct {
	//int func_num;
	double rho,phi,z,zs,delta;
	double complex h;
	Sys * S;
	double complex (*f) (double complex, void *);
} ScatParm;



//bessel.c
extern void zbesj_wrap(double,double,double,int,int,double*,double*,int*,int*);
extern double complex besselJ(int order, double complex z);
extern void zbesy_wrap(double,double,double,int,int,double*,double*,int*, double*,double*,int*);
extern double complex besselY(int order, double complex z);
extern void zbesh_wrap(double ,double ,double ,int,int,int ,double *,double *,int *,int *);
extern double complex besselH1(int order, double complex z);
extern double complex besselH2(int order, double complex z);
// Legendre.c
//extern double gauss_legendre(int n, double (*f)(double,void*), void* data, double a, double b);
extern double complex gauss_legendre_complex (int n, double complex (*f)(double ,void * ), void * data, double a, double  b);

//mem.c
extern double complex ** makefield_doublecomplex(int row, int col);

//fun.c
extern inline double MPI();

int index_finder(Sys *S,double z);
void mul_mat(double complex * M1, double complex * M2, double complex * Mout);
double complex ellipse_f(double t, void *);

extern inline double complex kz(Sys * S, const int i, const double complex krho);
extern inline double complex Fres(Sys * S,const int q, const int i, const int j, const double complex krho);
extern inline double complex Gamma(Sys * S,const int q,const int i, const int j, const double complex krho);
extern int index_finder(Sys *S,double z);
extern double complex ellipse_f(double t, void *);

void M1(Sys * S, const int q,const int i,const int j,const double complex krho, double complex * Mat);
//void M2(Sys * S, const int q,const int i,const int j,const double complex krho, double complex * Mat);
void Mlow(Sys *S,const double complex krho,const int q,double complex *Ml,int n);
void Mhigh(Sys *S,const double complex krho,const int q,double complex *Mh,int n);
extern void qspace(Sys *S,double complex *A,double complex *B,const double complex krho,const int q,const  int qw, const int l,const int ls,const double zs);

extern Sys * createSys(parm P);
extern double complex integral_imag(Sys *S, double rho,double phi,double z,double zs, double complex (*f)(double complex ,void * data));
extern double complex fsxx(double complex krho,void * p);
extern double complex fpxx(double complex krho, void * p);
extern double complex fpxx(double complex krho,void * p);
extern double complex fsxy(double complex krho,void * p);
extern double complex fpxy(double complex krho,void * p);
extern double complex fsyx(double complex krho,void * p);
extern double complex fpyx(double complex krho,void * p);
extern double complex fsxz(double complex krho,void * p);
extern double complex fpxz(double complex krho,void * p);
extern double complex fszx(double complex krho,void * p);
extern double complex fpzx(double complex krho,void * p);
extern double complex fsyy(double complex krho,void * p);
extern double complex fpyy(double complex krho,void * p);
extern double complex fsyz(double complex krho,void * p);
extern double complex fpyz(double complex krho,void * p);
extern double complex fszy(double complex krho,void * p);
extern double complex fpzy(double complex krho,void * p);
extern double complex fszz(double complex krho,void * p);
extern double complex fpzz(double complex krho,void * p);

//scattering.c
double complex infinite_f(double t, void *);
extern double complex integral_real(Sys *S, double rho, double phi ,double z,double zs, double delta,double h, double complex (*f)(double complex ,void * data));
/*
extern double complex Msxx(double complex krho,void * p);
extern double complex Mpxx(double complex krho, void * p);
extern double complex Mpxx(double complex krho,void * p);
extern double complex Msxy(double complex krho,void * p);
extern double complex Mpxy(double complex krho,void * p);
extern double complex Msyx(double complex krho,void * p);
extern double complex Mpyx(double complex krho,void * p);
extern double complex Msxz(double complex krho,void * p);
extern double complex Mpxz(double complex krho,void * p);
extern double complex Mszx(double complex krho,void * p);
extern double complex Mpzx(double complex krho,void * p);
extern double complex Msyy(double complex krho,void * p);
extern double complex Mpyy(double complex krho,void * p);
extern double complex Msyz(double complex krho,void * p);
extern double complex Mpyz(double complex krho,void * p);
extern double complex Mszy(double complex krho,void * p);
extern double complex Mpzy(double complex krho,void * p);
extern double complex Mszz(double complex krho,void * p);
extern double complex Mpzz(double complex krho,void * p);
*/

extern double complex Matsxx(double complex krho,void * p);
extern double complex Matpxx(double complex krho,void * p);
extern double complex Matpxx(double complex krho,void * p);
extern double complex Matsxy(double complex krho,void * p);
extern double complex Matpxy(double complex krho,void * p);
extern double complex Matsyx(double complex krho,void * p);
extern double complex Matpyx(double complex krho,void * p);
extern double complex Matsxz(double complex krho,void * p);
extern double complex Matpxz(double complex krho,void * p);
extern double complex Matszx(double complex krho,void * p);
extern double complex Matpzx(double complex krho,void * p);
extern double complex Matsyy(double complex krho,void * p);
extern double complex Matpyy(double complex krho,void * p);
extern double complex Matsyz(double complex krho,void * p);
extern double complex Matpyz(double complex krho,void * p);
extern double complex Matszy(double complex krho,void * p);
extern double complex Matpzy(double complex krho,void * p);
extern double complex Matszz(double complex krho,void * p);
extern double complex Matpzz(double complex krho,void * p);