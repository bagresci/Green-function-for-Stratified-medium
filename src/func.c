#include "mkj.h"
//#include <string.h>
//q=0, s-polarization, q=1, p-polarizaiton & alpha=x or y, q=2, p-polarization & alpha=z
// krho,d,phi,rho,Eps[],N,ls,k_0(w)

inline double MPI(){
	return  3.1415926535897932384626433832795028841971693993751;
}

Sys * createSys(parm P){
	Sys * S;
	S= calloc(1, sizeof(Sys));
	if(S==NULL) {printf("initilisation fail! [1]\n"); exit(1);}
	S->Eps = makefield_doublecomplex(P.N,3);
	S->d=P.d;
	S->k_0= P.k0;
	S->N=P.N;
	return S;
}


inline double complex kz(Sys * S, const int i, const double complex krho){
	return I*csqrt(krho*krho-(S->Eps[i][1]*S->Eps[i][0]*S->k_0*S->k_0));
	//return csqrt(S->Eps[i][1]*S->Eps[i][0]*S->k_0*S->k_0-krho*krho);
}

inline double complex Fres(Sys * S,const int q, const int i, const int j, const double complex krho){
	return (S->Eps[j][q]*kz(S,i,krho)-S->Eps[i][q]*kz(S,j,krho))/(S->Eps[j][q]*kz(S,i,krho)+S->Eps[i][q]*kz(S,j,krho));
}

inline double complex Gamma(Sys * S,const int q,const int i, const int j, const double complex krho){
	return (S->Eps[j][q]*kz(S,i,krho)+S->Eps[i][q]*kz(S,j,krho))/(2*S->Eps[j][q]*kz(S,j,krho))
	*((q-1)*(q-2)/2-kz(S,j,krho)*S->Eps[i][0]*q*(q-2)/(kz(S,i,krho)*S->Eps[j][0])+S->Eps[i][0]*q*(q-1)/(2*S->Eps[j][0]));
}


int index_finder(Sys *S,double z){
	if((S->N)<=1){return 1;}
	int i=0;
	while((S->d[i])>=z && i<(S->N)-1)i++;
	return i+1;
}


void mul_mat(double complex * M1, double complex * M2, double complex * Mout){
	//double complex M[2][2];
	Mout[0]=M1[0]*M2[0]+M1[1]*M2[2];
	Mout[1]=M1[0]*M2[1]+M1[1]*M2[3];
	Mout[2]=M1[2]*M2[0]+M1[3]*M2[2];
	Mout[3]=M1[2]*M2[1]+M1[3]*M2[3];
}

double complex ellipse_f(double t, void * data){
	ItgParm * P; P = (ItgParm *) data;/////doubt
	const double rho = P->rho;
	const double phi = P->phi;
	const double z = P->z;
	const double zs = P->zs;
	Sys * S; S=P->S;
	//f = P->f;
	double kmaj, kmax=S->k_0;
	int i=0;
	while(i<S->N){
		kmax=(kmax<(S->k_0)*creal(csqrt(S->Eps[i][1])))? (S->k_0)*creal(csqrt(S->Eps[i][1])):kmax;
		i++;
	}
	kmaj=(S->k_0+kmax)*0.5;
	return MPI()*kmaj*0.5*(cos(MPI()*t*0.5)+I*sin(MPI()*t*0.5)*0.01)*P->f(kmaj+kmaj*sin(MPI()*t*0.5)-I*cos(MPI()*t*0.5)*kmaj*0.01, P);

}


double complex imag_route1(double t,void * data){
	ItgParm * P; P = (ItgParm *) data;/////doubt
	Sys * S; S=P->S;
	double kmaj, kmax=S->k_0;
	int i=0;
	while(i<S->N){
		kmax=(kmax<(S->k_0)*creal(csqrt(S->Eps[i][1])))? (S->k_0)*creal(csqrt(S->Eps[i][1])):kmax;
		i++;
	}
	kmaj=(S->k_0+kmax);
	return 2.0*I/(1.0+t)/(1.0+t)*P->f(kmaj+I*(1.0-t)/(1.0+t),data);
}

double complex imag_route2(double t,void * data){
	ItgParm * P; P = (ItgParm *) data;/////doubt
	Sys * S; S=P->S;
	double kmaj, kmax=S->k_0;
	int i=0;
	while(i<S->N){
		kmax=(kmax<(S->k_0)*creal(csqrt(S->Eps[i][1])))? (S->k_0)*creal(csqrt(S->Eps[i][1])):kmax;
		i++;
	}
	kmaj=(S->k_0+kmax);
	return -2.0*I/(1.0-t)/(1.0-t)*P->f(kmaj+I*(1.0+t)/(t-1.0),data);
}

double complex integral_imag(Sys *S, double rho,double phi,double z,double zs, double complex (*f)(double complex ,void * data)){
	ItgParm P; P=(ItgParm){.rho = rho, .phi =phi, .z=z, .zs=zs, .S=S, .f=f};
	double complex a,b,c;
//	if(rho==0){
//		return gauss_legendre_complex(100,ellipse_f,&P,-1,1)+gauss_legendre_complex(100,infinite_f,&P,-1,1);
//	}
	P.func_num=1;a=gauss_legendre_complex(1024,ellipse_f, &P, -1.0, 1.0);//printf("%.15f,%.15f/",creal(a),cimag(a));
	P.func_num=2;b=gauss_legendre_complex(1024,imag_route1,&P,-1.0,1.0);//printf("%.15f,%.15f/",creal(b),cimag(b));
	P.func_num=3;c=gauss_legendre_complex(1024,imag_route2,&P,-1.0,1.0);//printf("%.15f,%.15f\n",creal(c),cimag(c));
	//return a+b+c;
	return a+creal(b+c)+I*cimag(b+c);
	//return a;
}


void M1(Sys * S, const int q,const int i,const int j,const double complex krho, double complex * Mat){
	int t;
	t=(i<j)?i:j;
	double a=S->d[t];
	double complex common = Gamma(S,q, i ,j,krho);
	Mat[0]= common * cexp(I*a*(kz(S,j,krho)-kz(S,i,krho)));
	Mat[1]= common * Fres(S,q,i,j,krho)*cexp(-I*a*(kz(S,j,krho)+kz(S,i,krho)));
	Mat[2]= common * Fres(S,q,i,j,krho)*cexp(I*a*(kz(S,j,krho)+kz(S,i,krho)));
	Mat[3]= common * cexp(-I*a*(kz(S,j,krho)-kz(S,i,krho)));
}

void Mlow(Sys *S,const double complex krho,const int q, double complex *Ml,int n){
	double complex A[4]={1,0,0,1};
	double complex B[4];
	double complex m[4];
	for(int i=0;i<n-1;i++){
		M1(S,q,i+1,i,krho,m);
		mul_mat(m,A,B);
		memcpy(A,B,sizeof(double complex)*4);
	}
	memcpy(Ml,A,sizeof(double complex)*4);
}

void Mhigh(Sys *S,const double complex krho,const int q, double complex *Mh,int n){
	double complex A[4]={1,0,0,1};
	double complex B[4];
	double complex m[4];
	for(int i=0;i<S->N-n;i++){
		M1(S,q,S->N-i-2,S->N-i-1,krho,m);
		mul_mat(m,A,B);
		memcpy(A,B,sizeof(double complex)*4);
	}
	memcpy(Mh,A,sizeof(double complex)*4);
}


// qw=1,s-pol or p-pol&beta=z, else,p-pol&beta=x,y
void qspace(Sys *S,double complex *A,double complex *B,const double complex krho,const int q,const  int qw,const int l,const int ls,const double zs){
	double complex A1,BN;
	double complex Ml[4], Mh[4],C[4];
	Mlow(S,krho,q,Ml,ls);
	Mhigh(S,krho,q,Mh,ls);
	if(qw==1){
		A1=(Mh[1]*cexp(I*kz(S,ls-1,krho)*zs)/Mh[3]+cexp(-I*kz(S,ls-1,krho)*zs))/(1-Ml[2]*Mh[1]/(Ml[0]*Mh[3]))/Ml[0];
		BN=(Ml[2]*cexp(-I*kz(S,ls-1,krho)*zs)/Ml[0]+cexp(I*kz(S,ls-1,krho)*zs))/(1-Ml[2]*Mh[1]/(Ml[0]*Mh[3]))/Mh[3];
	}
	else{
		A1=-(Mh[1]*cexp(I*kz(S,ls-1,krho)*zs)/Mh[3]-cexp(-I*kz(S,ls-1,krho)*zs))/(1-Ml[2]*Mh[1]/(Ml[0]*Mh[3]))/Ml[0];
		BN=-(Ml[2]*cexp(-I*kz(S,ls-1,krho)*zs)/Ml[0]-cexp(I*kz(S,ls-1,krho)*zs))/(1-Ml[2]*Mh[1]/(Ml[0]*Mh[3]))/Mh[3];
	}
	if(l<=ls){
		Mlow(S,krho,q,C,l);
		*A=C[0]*A1;
		*B=C[2]*A1;
	}
	else{
		Mhigh(S,krho,q,C,l);
		*A=C[1]*BN;
		*B=C[3]*BN;
	}
}

double complex fsxx(double complex krho,void * data){
	ItgParm * P; P = (ItgParm *) data;
	double rho=P->rho;
	double phi=P->phi;
	double z=P->z;
	double zs=P->zs;
	Sys *S; S=P->S;
	int num=P->func_num;
	double complex A,B;
	int l=index_finder(S,z);
	int ls=index_finder(S,zs);
	qspace(S,&A,&B,krho,0,1,l,ls,zs);
	if(l==ls && z<zs){
		A-=cexp(-I*kz(S,ls-1,krho)*zs);
		B+=cexp(I*kz(S,ls-1,krho)*zs);
	}
	//printf("%f,%f,%f,%f\n",creal(A),cimag(A),creal(B),cimag(B));
	switch(num){
		case 1:
			if(rho==0.0){
				return krho*0.5*(A*cexp(I*kz(S,l-1,krho)*z)+B*cexp(-I*kz(S,l-1,krho)*z))/kz(S,l-1,krho);
			}
			return (krho*besselJ(0,krho*rho)*sin(phi)*sin(phi)+besselJ(1,krho*rho)*cos(2*phi)/rho)*(A*cexp(I*kz(S,l-1,krho)*z)+B*cexp(-I*kz(S,l-1,krho)*z))/kz(S,l-1,krho);
		case 2:
			if(rho==0.0){
				return 0.5*krho*0.5*(A*cexp(I*kz(S,l-1,krho)*z)+B*cexp(-I*kz(S,l-1,krho)*z))/kz(S,l-1,krho);
			}
			return 0.5*(krho*besselH1(0,krho*rho)*sin(phi)*sin(phi)+besselH1(1,krho*rho)*cos(2*phi)/rho)*(A*cexp(I*kz(S,l-1,krho)*z)+B*cexp(-I*kz(S,l-1,krho)*z))/kz(S,l-1,krho);	
		case 3:
			if(rho==0.0){
				return 0.5*krho*0.5*(A*cexp(I*kz(S,l-1,krho)*z)+B*cexp(-I*kz(S,l-1,krho)*z))/kz(S,l-1,krho);
			}
			return 0.5*(krho*besselH2(0,krho*rho)*sin(phi)*sin(phi)+besselH2(1,krho*rho)*cos(2*phi)/rho)*(A*cexp(I*kz(S,l-1,krho)*z)+B*cexp(-I*kz(S,l-1,krho)*z))/kz(S,l-1,krho);
	}
}

double complex fpxx(double complex krho, void * p){
	ItgParm * P;
	P = (ItgParm *)p;
	double rho=P->rho;
	double phi=P->phi;
	double z=P->z;
	double zs=P->zs;
	Sys *S; S=P->S;
	int num=P->func_num;
	double complex A,B;
	int l=index_finder(S,z);
	int ls=index_finder(S,zs);

	qspace(S,&A,&B,krho,1,2,l,ls,zs);
	if(l==ls && z<zs){
		A=-(A-cexp(-I*kz(S,ls-1,krho)*zs));
		B=-B+cexp(I*kz(S,ls-1,krho)*zs);
	}
	switch(num){
		case 1:
			if(rho==0.0){
				return (z<zs? -1:1)*kz(S,l-1,krho)*krho*0.5*(A*cexp(I*kz(S,l-1,krho)*z)-B*cexp(-I*kz(S,l-1,krho)*z))/(S->k_0)/(S->k_0)/(S->Eps[l-1][1])/(S->Eps[l-1][0]);
			}
			return (z<zs? -1:1)*kz(S,l-1,krho)*(krho*besselJ(0,krho*rho)*cos(phi)*cos(phi)-besselJ(1,krho*rho)*cos(2*phi)/rho)*(A*cexp(I*kz(S,l-1,krho)*z)-B*cexp(-I*kz(S,l-1,krho)*z))/(S->k_0)/(S->k_0)/(S->Eps[l-1][1])/(S->Eps[l-1][0]);
		case 2:
			if(rho==0.0){
				return 0.5*(z<zs? -1:1)*kz(S,l-1,krho)*krho*0.5*(A*cexp(I*kz(S,l-1,krho)*z)-B*cexp(-I*kz(S,l-1,krho)*z))/(S->k_0)/(S->k_0)/(S->Eps[l-1][1])/(S->Eps[l-1][0]);
			}

			return 0.5*(z<zs? -1:1)*kz(S,l-1,krho)*(krho*besselH1(0,krho*rho)*cos(phi)*cos(phi)-besselH1(1,krho*rho)*cos(2*phi)/rho)*(A*cexp(I*kz(S,l-1,krho)*z)-B*cexp(-I*kz(S,l-1,krho)*z))/(S->k_0)/(S->k_0)/(S->Eps[l-1][1])/(S->Eps[l-1][0]);
		case 3:
			if(rho==0.0){
				return 0.5*(z<zs? -1:1)*kz(S,l-1,krho)*krho*0.5*(A*cexp(I*kz(S,l-1,krho)*z)-B*cexp(-I*kz(S,l-1,krho)*z))/(S->k_0)/(S->k_0)/(S->Eps[l-1][1])/(S->Eps[l-1][0]);
			}
			return 0.5*(z<zs? -1:1)*kz(S,l-1,krho)*(krho*besselH2(0,krho*rho)*cos(phi)*cos(phi)-besselH2(1,krho*rho)*cos(2*phi)/rho)*(A*cexp(I*kz(S,l-1,krho)*z)-B*cexp(-I*kz(S,l-1,krho)*z))/(S->k_0)/(S->k_0)/(S->Eps[l-1][1])/(S->Eps[l-1][0]);
	}
}


double complex fsxy(double complex krho, void * data){
	ItgParm * P;
	P = (ItgParm *)data;
	double rho=P->rho;
	double phi=P->phi;
	double z=P->z;
	double zs=P->zs;
	Sys *S; S=P->S;
	int num=P->func_num;
	if(rho==0.0){return 0.0;}
	double complex A, B;
	int l=index_finder(S,z);
	int ls=index_finder(S,zs);
	qspace(S,&A,&B,krho,0,1,l,ls,zs);
	if(l==ls && z<zs){
		A-=cexp(-I*kz(S,ls-1,krho)*zs);
		B+=cexp(I*kz(S,ls-1,krho)*zs);
	}
	switch(num){
		case 1:
			return (-krho*besselJ(0,krho*rho)+besselJ(1,krho*rho)*2/rho)*sin(2*phi)*0.5*(A*cexp(I*kz(S,l-1,krho)*z)+B*cexp(-I*kz(S,l-1,krho)*z))/kz(S,l-1,krho);
		case 2:
			return 0.5*(-krho*besselH1(0,krho*rho)+besselH1(1,krho*rho)*2/rho)*sin(2*phi)*0.5*(A*cexp(I*kz(S,l-1,krho)*z)+B*cexp(-I*kz(S,l-1,krho)*z))/kz(S,l-1,krho);
		case 3:
			return 0.5*(-krho*besselH2(0,krho*rho)+besselH2(1,krho*rho)*2/rho)*sin(2*phi)*0.5*(A*cexp(I*kz(S,l-1,krho)*z)+B*cexp(-I*kz(S,l-1,krho)*z))/kz(S,l-1,krho);
	}
}


double complex fpxy(double complex krho, void * data){
	ItgParm * P;
	P = (ItgParm *)data;
	double rho=P->rho;
	double phi=P->phi;
	double z=P->z;
	double zs=P->zs;
	Sys *S; S=P->S;
	int num=P->func_num;
	double complex A, B;
	int l=index_finder(S,z);
	int ls=index_finder(S,zs);
	if(rho==0.0){return 0.0;}
	qspace(S,&A,&B,krho,1,2,l,ls,zs);
	if(l==ls && z<zs){
		A=-(A-cexp(-I*kz(S,ls-1,krho)*zs));
		B=-B+cexp(I*kz(S,ls-1,krho)*zs);
	}
	switch(num){
		case 1:
			return (z<zs? -1:1)*kz(S,l-1,krho)*(krho*besselJ(0,krho*rho)-besselJ(1,krho*rho)*2/rho)*sin(2*phi)*0.5*(A*cexp(I*kz(S,l-1,krho)*z)-B*cexp(-I*kz(S,l-1,krho)*z))/(S->k_0)/(S->k_0)/(S->Eps[l-1][1])/(S->Eps[l-1][0]);
		case 2:
			return 0.5*(z<zs? -1:1)*kz(S,l-1,krho)*(krho*besselH1(0,krho*rho)-besselH1(1,krho*rho)*2/rho)*sin(2*phi)*0.5*(A*cexp(I*kz(S,l-1,krho)*z)-B*cexp(-I*kz(S,l-1,krho)*z))/(S->k_0)/(S->k_0)/(S->Eps[l-1][1])/(S->Eps[l-1][0]);
		case 3:
			return 0.5*(z<zs? -1:1)*kz(S,l-1,krho)*(krho*besselH2(0,krho*rho)-besselH2(1,krho*rho)*2/rho)*sin(2*phi)*0.5*((A)*cexp(I*kz(S,l-1,krho)*z)-(B)*cexp(-I*kz(S,l-1,krho)*z))/(S->k_0)/(S->k_0)/(S->Eps[l-1][1])/(S->Eps[l-1][0]);
	}
}


inline double complex fsxz(double complex krho, void * p){
	return 0;
}


double complex fpxz(double complex krho,void * data){
	ItgParm * P;
	P = (ItgParm *)data;
	double rho=P->rho;
	double phi=P->phi;
	double z=P->z;
	double zs=P->zs;
	Sys *S; S=P->S;
	int num=P->func_num;
	double complex A, B;
	int l=index_finder(S,z);
	int ls=index_finder(S,zs);
	if(rho==0.0){return 0.0;}
	qspace(S,&A,&B,krho,2,1,l,ls,zs);
	if(l==ls && z<zs){
		A-=cexp(-I*kz(S,ls-1,krho)*zs);
		B+=cexp(I*kz(S,ls-1,krho)*zs);
	}
	switch(num){
		case 1:
			return -I*krho*krho/(S->k_0)/(S->k_0)/(S->Eps[l-1][1])/(S->Eps[l-1][0])*besselJ(1,krho*rho)*cos(phi)*(A*cexp(I*kz(S,l-1,krho)*z)-B*cexp(-I*kz(S,l-1,krho)*z));
		case 2:
			return -0.5*I*krho*krho/(S->k_0)/(S->k_0)/(S->Eps[l-1][1])/(S->Eps[l-1][0])*besselH1(1,krho*rho)*cos(phi)*(A*cexp(I*kz(S,l-1,krho)*z)-B*cexp(-I*kz(S,l-1,krho)*z));
		case 3:
			return -0.5*I*krho*krho/(S->k_0)/(S->k_0)/(S->Eps[l-1][1])/(S->Eps[l-1][0])*besselH2(1,krho*rho)*cos(phi)*(A*cexp(I*kz(S,l-1,krho)*z)-B*cexp(-I*kz(S,l-1,krho)*z));
	}
}


inline double complex fsyx(double complex krho, void * data){
	return fsxy(krho, data);
}

inline double complex fpyx(double complex krho, void * data){
	return fpxy(krho,data);
}

double complex fsyy(double complex krho, void * data){
	ItgParm * P;
	P = (ItgParm *)data;
	double rho=P->rho;
	double phi=P->phi;
	double z=P->z;
	double zs=P->zs;
	Sys *S; S=P->S;
	int num=P->func_num;
	double complex A, B;
	int l=index_finder(S,z);
	int ls=index_finder(S,zs);
	qspace(S,&A,&B,krho,0,1,l,ls,zs);
	if(l==ls && z<zs){
		A-=cexp(-I*kz(S,ls-1,krho)*zs);
		B+=cexp(I*kz(S,ls-1,krho)*zs);
	}
//	if(rho==0.0){
//		return krho*0.5*(A*cexp(I*kz(S,l-1,krho)*z)+B*cexp(-I*kz(S,l-1,krho)*z))/kz(S,l-1,krho);
//	}
	switch(num){
		case 1:
			if(rho==0.0){
				return krho*0.5*(A*cexp(I*kz(S,l-1,krho)*z)+B*cexp(-I*kz(S,l-1,krho)*z))/kz(S,l-1,krho);
			}
			return (krho*besselJ(0,krho*rho)*cos(phi)*cos(phi)-besselJ(1,krho*rho)*cos(2*phi)/rho)*(A*cexp(I*kz(S,l-1,krho)*z)+B*cexp(-I*kz(S,l-1,krho)*z))/kz(S,l-1,krho);
		case 2:
			if(rho==0.0){
				return 0.5*krho*0.5*(A*cexp(I*kz(S,l-1,krho)*z)+B*cexp(-I*kz(S,l-1,krho)*z))/kz(S,l-1,krho);
			}
			return 0.5*(krho*besselH1(0,krho*rho)*cos(phi)*cos(phi)-besselH1(1,krho*rho)*cos(2*phi)/rho)*(A*cexp(I*kz(S,l-1,krho)*z)+B*cexp(-I*kz(S,l-1,krho)*z))/kz(S,l-1,krho);
		case 3:
			if(rho==0.0){
				return 0.5*krho*0.5*(A*cexp(I*kz(S,l-1,krho)*z)+B*cexp(-I*kz(S,l-1,krho)*z))/kz(S,l-1,krho);
			}
			return 0.5*(krho*besselH2(0,krho*rho)*cos(phi)*cos(phi)-besselH2(1,krho*rho)*cos(2*phi)/rho)*(A*cexp(I*kz(S,l-1,krho)*z)+B*cexp(-I*kz(S,l-1,krho)*z))/kz(S,l-1,krho);
	}
}

double complex fpyy(double complex krho, void * data){
	ItgParm * P;
	P = (ItgParm *)data;
	double rho=P->rho;
	double phi=P->phi;
	double z=P->z;
	double zs=P->zs;
	Sys *S; S=P->S;
	int num=P->func_num;
	double complex A, B;
	int l=index_finder(S,z);
	int ls=index_finder(S,zs);
	qspace(S,&A,&B,krho,1,2,l,ls,zs);
	if(l==ls && z<zs){
		A=-(A-cexp(-I*kz(S,ls-1,krho)*zs));
		B=-B+cexp(I*kz(S,ls-1,krho)*zs);
	}
	switch(num){
		case 1:
			if(rho==0.0){
				return (z<zs? -1:1)*kz(S,l-1,krho)*krho*0.5*(A*cexp(I*kz(S,l-1,krho)*z)-B*cexp(-I*kz(S,l-1,krho)*z))/(S->k_0)/(S->k_0)/(S->Eps[l-1][1])/(S->Eps[l-1][0]);
			}
			return (z<zs? -1:1)*kz(S,l-1,krho)*(krho*besselJ(0,krho*rho)*sin(phi)*sin(phi)+besselJ(1,krho*rho)*cos(2*phi)/rho)*(A*cexp(I*kz(S,l-1,krho)*z)-B*cexp(-I*kz(S,l-1,krho)*z))/(S->k_0)/(S->k_0)/(S->Eps[l-1][1])/(S->Eps[l-1][0]);
		case 2:
			if(rho==0.0){
				return 0.5*(z<zs? -1:1)*kz(S,l-1,krho)*krho*0.5*(A*cexp(I*kz(S,l-1,krho)*z)-B*cexp(-I*kz(S,l-1,krho)*z))/(S->k_0)/(S->k_0)/(S->Eps[l-1][1])/(S->Eps[l-1][0]);
			}
			return 0.5*(z<zs? -1:1)*kz(S,l-1,krho)*(krho*besselH1(0,krho*rho)*sin(phi)*sin(phi)+besselH1(1,krho*rho)*cos(2*phi)/rho)*(A*cexp(I*kz(S,l-1,krho)*z)-B*cexp(-I*kz(S,l-1,krho)*z))/(S->k_0)/(S->k_0)/(S->Eps[l-1][1])/(S->Eps[l-1][0]);
		case 3:
			if(rho==0.0){
				return (z<zs? -1:1)*kz(S,l-1,krho)*krho*0.5*(A*cexp(I*kz(S,l-1,krho)*z)-B*cexp(-I*kz(S,l-1,krho)*z))/(S->k_0)/(S->k_0)/(S->Eps[l-1][1])/(S->Eps[l-1][0]);
			}
			return 0.5*(z<zs? -1:1)*kz(S,l-1,krho)*(krho*besselH2(0,krho*rho)*sin(phi)*sin(phi)+besselH2(1,krho*rho)*cos(2*phi)/rho)*(A*cexp(I*kz(S,l-1,krho)*z)-B*cexp(-I*kz(S,l-1,krho)*z))/(S->k_0)/(S->k_0)/(S->Eps[l-1][1])/(S->Eps[l-1][0]);
	}
}


double complex fsyz(double complex krho, void * data){
	return 0;
}


double complex fpyz(double complex krho, void * p){
	ItgParm * P;
	P = (ItgParm *)p;
	double rho=P->rho;
	double phi=P->phi;
	double z=P->z;
	double zs=P->zs;
	Sys *S; S=P->S;
	int num=P->func_num;
	if(rho==0.0){return 0;}
	double complex A, B;
	int l=index_finder(S,z);
	int ls=index_finder(S,zs);
	qspace(S,&A,&B,krho,2,1,l,ls,zs);
	if(l==ls && z<zs){
		A-=cexp(-I*kz(S,ls-1,krho)*zs);
		B+=cexp(I*kz(S,ls-1,krho)*zs);
	}
	if(rho==0.0){return 0;}
	switch(num){
		case 1:
			return -I*krho*krho/(S->k_0)/(S->k_0)/(S->Eps[l-1][1])/(S->Eps[l-1][0])*besselJ(1,krho*rho)*sin(phi)*(A*cexp(I*kz(S,l-1,krho)*z)-B*cexp(-I*kz(S,l-1,krho)*z));
		case 2:
			return -0.5*I*krho*krho/(S->k_0)/(S->k_0)/(S->Eps[l-1][1])/(S->Eps[l-1][0])*besselH1(1,krho*rho)*sin(phi)*(A*cexp(I*kz(S,l-1,krho)*z)-B*cexp(-I*kz(S,l-1,krho)*z));
		case 3:
			return -0.5*I*krho*krho/(S->k_0)/(S->k_0)/(S->Eps[l-1][1])/(S->Eps[l-1][0])*besselH2(1,krho*rho)*sin(phi)*(A*cexp(I*kz(S,l-1,krho)*z)-B*cexp(-I*kz(S,l-1,krho)*z));
	}
}


double complex fszx(double complex krho, void * data){
	return 0;
}


double complex fpzx(double complex krho, void * p){
	ItgParm * P;
	P = (ItgParm *)p;
	double rho=P->rho;
	double phi=P->phi;
	double z=P->z;
	double zs=P->zs;
	Sys *S; S=P->S;
	int num=P->func_num;
	if(rho==0.0){return 0.0;}
	double complex A, B;
	int l=index_finder(S,z);
	int ls=index_finder(S,zs);
	qspace(S,&A,&B,krho,1,2,l,ls,zs);
	if(l==ls && z<zs){
		A=-(A-cexp(-I*kz(S,ls-1,krho)*zs));
		B=-B+cexp(I*kz(S,ls-1,krho)*zs);
	}
	switch(num){
		case 1 :
			return (z<zs? 1:-1)*I*krho*krho/(S->k_0)/(S->k_0)/(S->Eps[l-1][1])/(S->Eps[l-1][0])*besselJ(1,krho*rho)*cos(phi)*(A*cexp(I*kz(S,l-1,krho)*z)+B*cexp(-I*kz(S,l-1,krho)*z));
		case 2 :
			return 0.5*(z<zs? 1:-1)*I*krho*krho/(S->k_0)/(S->k_0)/(S->Eps[l-1][1])/(S->Eps[l-1][0])*besselH1(1,krho*rho)*cos(phi)*(A*cexp(I*kz(S,l-1,krho)*z)+B*cexp(-I*kz(S,l-1,krho)*z));
		case 3 :
			return 0.5*(z<zs? 1:-1)*I*krho*krho/(S->k_0)/(S->k_0)/(S->Eps[l-1][1])/(S->Eps[l-1][0])*besselH2(1,krho*rho)*cos(phi)*(A*cexp(I*kz(S,l-1,krho)*z)+B*cexp(-I*kz(S,l-1,krho)*z));
	}
}


double complex fszy(double complex krho, void * p){
	return 0;
}


double complex fpzy(double complex krho,void * p){
	ItgParm * P;
	P = (ItgParm *)p;
	double rho=P->rho;
	double phi=P->phi;
	double z=P->z;
	double zs=P->zs;
	Sys *S; S=P->S;
	int num=P->func_num;
	double complex A, B;
	int l=index_finder(S,z);
	int ls=index_finder(S,zs);
	if(rho==0.0){return 0.0;}
	qspace(S,&A,&B,krho,1,2,l,ls,zs);
	if(l==ls && z<zs){
		A=-(A-cexp(-I*kz(S,ls-1,krho)*zs));
		B=-B+cexp(I*kz(S,ls-1,krho)*zs);
	}
	switch(num){
		case 1 :
			return (z>=zs? -I: I) *krho*krho/(S->k_0)/(S->k_0)/(S->Eps[l-1][1])/(S->Eps[l-1][0])*besselJ(1,krho*rho)*sin(phi)*(A*cexp(I*kz(S,l-1,krho)*z)+B*cexp(-I*kz(S,l-1,krho)*z));
		case 2 :
			return 0.5*(z>=zs? -I: I) *krho*krho/(S->k_0)/(S->k_0)/(S->Eps[l-1][1])/(S->Eps[l-1][0])*besselH1(1,krho*rho)*sin(phi)*(A*cexp(I*kz(S,l-1,krho)*z)+B*cexp(-I*kz(S,l-1,krho)*z));
		case 3 :
			return 0.5*(z>=zs? -I: I) *krho*krho/(S->k_0)/(S->k_0)/(S->Eps[l-1][1])/(S->Eps[l-1][0])*besselH2(1,krho*rho)*sin(phi)*(A*cexp(I*kz(S,l-1,krho)*z)+B*cexp(-I*kz(S,l-1,krho)*z));
	}
}

double complex fszz(double complex krho, void * p){
	return 0;
}


double complex fpzz(double complex krho, void * p){
	ItgParm * P;
	P = (ItgParm *)p;
	double rho=P->rho;
	double phi=P->phi;
	double z=P->z;
	double zs=P->zs;
	Sys *S; S=P->S;
	int num= P->func_num;
	double complex A, B;
	int l=index_finder(S,z);
	int ls=index_finder(S,zs);
	qspace(S,&A,&B,krho,2,1,l,ls,zs);
	if(l==ls && z<zs){
		A-=cexp(-I*kz(S,ls-1,krho)*zs);
		B+=cexp(I*kz(S,ls-1,krho)*zs);
	}
	switch(num){
		case 1 :
			return (krho)*(krho)*(krho)*besselJ(0,krho*rho)*(A*cexp(I*kz(S,l-1,krho)*z)+B*cexp(-I*kz(S,l-1,krho)*z))/(S->k_0)/(S->k_0)/((S->Eps[l-1][1])*kz(S,l-1,krho))/(S->Eps[l-1][0]);
		case 2 :
			if(rho==0.0){
				return 0.5*(krho)*(krho)*(krho)*(A*cexp(I*kz(S,l-1,krho)*z)+B*cexp(-I*kz(S,l-1,krho)*z))/(S->k_0)/(S->k_0)/((S->Eps[l-1][1])*kz(S,l-1,krho))/(S->Eps[l-1][0]);
			}
			return 0.5*krho*krho*krho*besselH1(0,krho*rho)*(A*cexp(I*kz(S,l-1,krho)*z)+B*cexp(-I*kz(S,l-1,krho)*z))/(S->k_0)/(S->k_0)/((S->Eps[l-1][1])*kz(S,l-1,krho))/(S->Eps[l-1][0]);
		case 3 :
			if(rho==0.0){
				return 0.5*(krho)*(krho)*(krho)*(A*cexp(I*kz(S,l-1,krho)*z)+B*cexp(-I*kz(S,l-1,krho)*z))/(S->k_0)/(S->k_0)/((S->Eps[l-1][1])*kz(S,l-1,krho))/(S->Eps[l-1][0]);
			}
			return 0.5*(krho)*(krho)*(krho)*besselH2(0,krho*rho)*(A*cexp(I*kz(S,l-1,krho)*z)+B*cexp(-I*kz(S,l-1,krho)*z))/(S->k_0)/(S->k_0)/((S->Eps[l-1][1])*kz(S,l-1,krho))/(S->Eps[l-1][0]);
	}
}