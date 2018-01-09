#include "mkj.h"
#define PI 3.1415926535897932384626433832795028841971693993751


double complex angle1(double theta,void * data){
	ScatParm *P; P=(ScatParm *)data;
	double complex krho=P->h;
	double complex delta=P->delta;
	return 2.0*I*ccos(krho*(P->rho)*cos(theta-(P->phi)))/cos(theta)*sin(krho*delta*cos(theta))*sin(krho*delta*sin(theta));
}

double complex angle2(double theta,void * data){
	ScatParm *P; P=(ScatParm *)data;
	double complex krho=P->h;
	double complex delta=P->delta;
	return 2.0*I*csin(krho*(P->rho)*cos(theta-(P->phi)))/sin(theta)*sin(krho*delta*cos(theta))*sin(krho*delta*sin(theta));
}

double complex angle3(double theta,void * data){
	ScatParm *P; P=(ScatParm *)data;
	double complex krho=P->h;
	double complex delta=P->delta;
	return 4.0*(ccos(krho*(P->rho)*cos(theta-(P->phi)))+ccos(-krho*(P->rho)*sin(theta-(P->phi))))/sin(2.0*theta)*sin(krho*delta*cos(theta))*sin(krho*delta*sin(theta));
}

double complex angle4(double theta,void * data){
	ScatParm *P; P=(ScatParm *)data;
	double complex krho=P->h;
	double complex delta=P->delta;
	//printf("%f\n",cabs(krho));
	return 2.0*ccos(krho*(P->rho)*cos(theta-(P->phi)))*tan(theta)*sin(krho*delta*cos(theta))*sin(krho*delta*sin(theta));
}

double complex angle5(double theta,void * data){
	ScatParm *P; P=(ScatParm *)data;
	double complex krho=P->h;
	double complex delta=P->delta;
	return 2.0*ccos(krho*(P->rho)*cos(theta-(P->phi)))/cot(theta)*sin(krho*delta*cos(theta))*sin(krho*delta*sin(theta));
}

double complex ellip_f(double t, void * data){
	ScatParm * P; P = (ScatParm *) data;/////doubt
	Sys *S; S=P->S;
	
	//int num=P->func_num;
	double kmaj, kmax=S->k_0;
	int i=0;
	while(i<S->N){
		kmax=(kmax<(S->k_0)*creal(csqrt(S->Eps[i][1])))? (S->k_0)*creal(csqrt(S->Eps[i][1])):kmax;
		i++;
	}
	kmaj=(S->k_0+kmax)*0.5;
	return MPI()*kmaj*0.5*(cos(MPI()*t*0.5)+I*sin(MPI()*t*0.5)*0.1)*P->f(kmaj+kmaj*sin(MPI()*t*0.5)-I*cos(MPI()*t*0.5)*kmaj*0.1, data);

}

double complex infinite_f(double t, void * data){
	ScatParm * P; P = (ScatParm *) data;/////doubt
	Sys * S; S=P->S;
	double kmaj, kmax=S->k_0;
	int i=0;
	while(i<S->N){
		kmax=(kmax<(S->k_0)*creal(csqrt(S->Eps[i][1])))? (S->k_0)*creal(csqrt(S->Eps[i][1])):kmax;
		i++;
	}
	kmaj=(S->k_0+kmax)*0.5;
	return 2.0/(1.0+t)/(1.0+t)*P->f(2.0*kmaj+(1.0-t)/(1.0+t),data);

}


double complex integral_real(Sys *S, double rho, double phi, double z,double zs, double delta,double h, double complex (*f)(double complex ,void * data)){
	double complex a,b;
	double kmaj, kmax=S->k_0;
	ScatParm P;
	int i=0;
	while(i<S->N){
		kmax=(kmax<(S->k_0)*creal(csqrt(S->Eps[i][1])))? (S->k_0)*creal(csqrt(S->Eps[i][1])):kmax;
		i++;
	}
	kmaj=(S->k_0+kmax)*0.5;
	P=(ScatParm){.rho=rho, .phi=phi, .z=z, .zs=zs, .delta=delta, .h=h, .S=S, .f=f};

	//P.func_num=1;
	a=gauss_legendre_complex(1024,ellip_f, &P, -1.0, 1.0);//printf("%.15f\n",cabs(a));
	b=gauss_legendre_complex(512,infinite_f,&P,-1.0,1.0);
	return I*0.125/MPI()/MPI()*(a+b);

}

double complex Matsxx(double complex krho, void *data){
	ScatParm * P; P = (ScatParm *) data;
	double complex value,angle_part;
	double complex h=P->h;
	double delta=P->delta;
	double z=P->z;
	double zs=P->zs;
	Sys *S; S=P->S;
	P->h=krho;
	int l=index_finder(S,z);
	angle_part=gauss_legendre_complex(512,angle4, P, -PI*0.5, PI*0.5);
	switch(l){
		case 1:
			value=-I*cexp(-I*S->d[0]*(kz(S,0,krho)-kz(S,1,krho)))/kz(S,1,krho)/Gamma(S,0,1,0,krho)/(1-Fres(S,0,1,0,krho)*Fres(S,0,1,2,krho)*cexp(2.0*I*kz(S,1,krho)*(S->d[0]-S->d[1])))*(Fres(S,0,1,2,krho)*(cexp(I*kz(S,1,krho)*(z+zs+h-2.0*S->d[1]))-cexp(I*kz(S,1,krho)*(z+zs-h-2.0*S->d[1])))+cexp(I*kz(S,1,krho)*(z+h-zs))-cexp(I*kz(S,1,krho)*(z-h-zs)));
			break;
		case 2:
			value=-I/kz(S,1,krho)/(1-Fres(S,0,1,0,krho)*Fres(S,0,1,2,krho)*cexp(2.0*I*kz(S,1,krho)*(S->d[0]-S->d[1])))*(Fres(S,0,1,0,krho)*Fres(S,0,1,2,krho)*(cexp(I*kz(S,1,krho)*(2.0*(S->d[0]-S->d[1])+z+h-zs))+cexp(I*kz(S,1,krho)*(2.0*(S->d[0]-S->d[1])-z+h+zs))-cexp(I*kz(S,1,krho)*(2.0*(S->d[0]-S->d[1])+z-h-zs))-cexp(I*kz(S,1,krho)*(2.0*(S->d[0]-S->d[1])-z+zs-h)))+Fres(S,0,1,2,krho)*(cexp(I*kz(S,1,krho)*(zs+z+h-2.0*S->d[1]))-cexp(I*kz(S,1,krho)*(zs+z-h-2.0*S->d[1])))+Fres(S,0,1,0,krho)*(cexp(I*kz(S,1,krho)*(-zs-z+h+2.0*S->d[0]))-cexp(I*kz(S,1,krho)*(-zs-z-h+2.0*S->d[0]))));
			value=value-I/kz(S,1,krho)*(z==zs? 2.0*(cexp(I*kz(S,1,krho)*h)-1):(cexp(I*kz(S,1,krho)*(abs(z-zs)+h))-cexp(I*kz(S,1,krho)*(abs(z-zs)-h))));
			//printf("%f\n",cabs(value));
			break;
		case 3:
			value=-I*cexp(-I*S->d[1]*(kz(S,2,krho)-kz(S,1,krho)))/kz(S,1,krho)/Gamma(S,0,1,2,krho)/(1-Fres(S,0,1,0,krho)*Fres(S,0,1,2,krho)*cexp(2.0*I*kz(S,1,krho)*(S->d[0]-S->d[1])))*(Fres(S,0,1,0,krho)*(cexp(I*kz(S,1,krho)*(-z-zs+h+2.0*S->d[0]))-cexp(I*kz(S,1,krho)*(-z-zs-h+2.0*S->d[0])))+cexp(I*kz(S,1,krho)*(-z+h+zs))-cexp(I*kz(S,1,krho)*(-z-h+zs)));
			break;
	}
	return -1.0*value/kz(S,l-1,krho)/krho*angle_part;
}

double complex Matpxx(double complex krho, void *data){
	ScatParm * P; P = (ScatParm *) data;
	double complex value,angle_part;
	double h=P->h;
	double delta=P->delta;
	double z=P->z;
	double zs=P->zs;
	Sys *S; S=P->S;
	P->h=krho;
	int l=index_finder(S,z);
	angle_part=gauss_legendre_complex(512,angle5, P, 0.0, PI);
	switch(l){
		case 1:
			value=+I*cexp(-I*S->d[0]*(kz(S,0,krho)-kz(S,1,krho)))/kz(S,1,krho)/Gamma(S,1,1,0,krho)/(1-Fres(S,1,1,0,krho)*Fres(S,1,1,2,krho)*cexp(2.0*I*kz(S,1,krho)*(S->d[0]-S->d[1])))*(Fres(S,1,1,2,krho)*(cexp(I*kz(S,1,krho)*(z+zs+h-2.0*S->d[1]))-cexp(I*kz(S,1,krho)*(z+zs-h-2.0*S->d[1])))-cexp(I*kz(S,1,krho)*(z+h-zs))+cexp(I*kz(S,1,krho)*(z-h-zs)));
			break;
		case 2:
			value=-I/kz(S,1,krho)/(1-Fres(S,1,1,0,krho)*Fres(S,1,1,2,krho)*cexp(2.0*I*kz(S,1,krho)*(S->d[0]-S->d[1])))*(Fres(S,1,1,0,krho)*Fres(S,1,1,2,krho)*(cexp(I*kz(S,1,krho)*(2.0*(S->d[0]-S->d[1])+z+h-zs))+cexp(I*kz(S,1,krho)*(2.0*(S->d[0]-S->d[1])-z+h+zs))-cexp(I*kz(S,1,krho)*(2.0*(S->d[0]-S->d[1])+z-h-zs))-cexp(I*kz(S,1,krho)*(2.0*(S->d[0]-S->d[1])-z+zs-h)))-Fres(S,1,1,2,krho)*(cexp(I*kz(S,1,krho)*(zs+z+h-2.0*S->d[1]))-cexp(I*kz(S,1,krho)*(zs+z-h-2.0*S->d[1])))-Fres(S,1,1,0,krho)*(cexp(I*kz(S,1,krho)*(-zs-z+h+2.0*S->d[0]))-cexp(I*kz(S,1,krho)*(-zs-z-h+2.0*S->d[0]))));
			value=value-I/kz(S,1,krho)*(z==zs? 2.0*(cexp(I*kz(S,1,krho)*h)-1):(cexp(I*kz(S,1,krho)*(abs(z-zs)+h))-cexp(I*kz(S,1,krho)*(abs(z-zs)-h))));
			break;
		case 3:
			value=+I*cexp(-I*S->d[1]*(kz(S,2,krho)-kz(S,1,krho)))/kz(S,1,krho)/Gamma(S,1,1,2,krho)/(1-Fres(S,1,1,0,krho)*Fres(S,1,1,2,krho)*cexp(2.0*I*kz(S,1,krho)*(S->d[0]-S->d[1])))*(Fres(S,0,1,0,krho)*(cexp(I*kz(S,1,krho)*(-z-zs+h+2.0*S->d[0]))-cexp(I*kz(S,1,krho)*(-z-zs-h+2.0*S->d[0])))-cexp(I*kz(S,1,krho)*(-z+h+zs))+cexp(I*kz(S,1,krho)*(-z-h+zs)));
			break;
	}
	return -1.0*value*angle_part*kz(S,l-1,krho)/(S->k_0)/(S->k_0)/(S->Eps[l-1][1])/(S->Eps[l-1][0])/krho;
}

double complex Matsxy(double complex krho, void *data){
	ScatParm * P; P = (ScatParm *) data;
	double complex value;
	double h=P->h;
	double delta=P->delta;
	double x=(P->rho)*cos(P->phi);
	double y=(P->rho)*sin(P->phi);
	double z=P->z;
	double zs=P->zs;
	Sys *S; S=P->S;
	int l=index_finder(S,z);
	double rho1=sqrt((x-delta)*(x-delta)+(y-delta)*(y-delta));
	double theta1=abs(y-delta)/(y-delta)*acos((x-delta)/rho1);
	double rho2=sqrt((x+delta)*(x+delta)+(y+delta)*(y+delta));
	double theta2=abs(y+delta)/(y+delta)*acos((x+delta)/rho2);
	double rho3=sqrt((x-delta)*(x-delta)+(y+delta)*(y+delta));
	double theta3=abs(y+delta)/(y+delta)*acos((x-delta)/rho3);
	double rho4=sqrt((x+delta)*(x+delta)+(y-delta)*(y-delta));
	double theta4=abs(y-delta)/(y-delta)*acos((x+delta)/rho4);
	switch(l){
		case 1:
			value=-I*cexp(-I*S->d[0]*(kz(S,0,krho)-kz(S,1,krho)))/kz(S,1,krho)/Gamma(S,0,1,0,krho)/(1-Fres(S,0,1,0,krho)*Fres(S,0,1,2,krho)*cexp(2.0*I*kz(S,1,krho)*(S->d[0]-S->d[1])))*(Fres(S,0,1,2,krho)*(cexp(I*kz(S,1,krho)*(z+zs+h-2.0*S->d[1]))-cexp(I*kz(S,1,krho)*(z+zs-h-2.0*S->d[1])))+cexp(I*kz(S,1,krho)*(z+h-zs))-cexp(I*kz(S,1,krho)*(z-h-zs)));
			break;
		case 2:
			value=-I/kz(S,1,krho)/(1-Fres(S,0,1,0,krho)*Fres(S,0,1,2,krho)*cexp(2.0*I*kz(S,1,krho)*(S->d[0]-S->d[1])))*(Fres(S,0,1,0,krho)*Fres(S,0,1,2,krho)*(cexp(I*kz(S,1,krho)*(2.0*(S->d[0]-S->d[1])+z+h-zs))+cexp(I*kz(S,1,krho)*(2.0*(S->d[0]-S->d[1])-z+h+zs))-cexp(I*kz(S,1,krho)*(2.0*(S->d[0]-S->d[1])+z-h-zs))-cexp(I*kz(S,1,krho)*(2.0*(S->d[0]-S->d[1])-z+zs-h)))+Fres(S,0,1,2,krho)*(cexp(I*kz(S,1,krho)*(zs+z+h-2.0*S->d[1]))-cexp(I*kz(S,1,krho)*(zs+z-h-2.0*S->d[1])))+Fres(S,0,1,0,krho)*(cexp(I*kz(S,1,krho)*(-zs-z+h+2.0*S->d[0]))-cexp(I*kz(S,1,krho)*(-zs-z-h+2.0*S->d[0]))));
			value=value-I/kz(S,1,krho)*(z==zs? 2.0*(cexp(I*kz(S,1,krho)*h)-1):(cexp(I*kz(S,1,krho)*(abs(z-zs)+h))-cexp(I*kz(S,1,krho)*(abs(z-zs)-h))));
			break;
		case 3:
			value=-I*cexp(-I*S->d[1]*(kz(S,2,krho)-kz(S,1,krho)))/kz(S,1,krho)/Gamma(S,0,1,2,krho)/(1-Fres(S,0,1,0,krho)*Fres(S,0,1,2,krho)*cexp(2.0*I*kz(S,1,krho)*(S->d[0]-S->d[1])))*(Fres(S,0,1,0,krho)*(cexp(I*kz(S,1,krho)*(-z-zs+h+2.0*S->d[0]))-cexp(I*kz(S,1,krho)*(-z-zs-h+2.0*S->d[0])))+cexp(I*kz(S,1,krho)*(-z+h+zs))-cexp(I*kz(S,1,krho)*(-z-h+zs)));
			break;
	}
	return 2.0*PI*value*(besselJ(0,krho*rho1)+besselJ(0,krho*rho2)-besselJ(0,krho*rho3)-besselJ(0,krho*rho4))/krho/kz(S,l-1,krho);
}

double complex Matpxy(double complex krho, void *data){
	ScatParm * P; P = (ScatParm *) data;
	double complex value;
	double h=P->h;
	double delta=P->delta;
	double x=(P->rho)*cos(P->phi);
	double y=(P->rho)*sin(P->phi);
	double z=P->z;
	double zs=P->zs;
	Sys *S; S=P->S;
	int l=index_finder(S,z);
	double rho1=sqrt((x-delta)*(x-delta)+(y-delta)*(y-delta));
	double theta1=abs(y-delta)/(y-delta)*acos((x-delta)/rho1);
	double rho2=sqrt((x+delta)*(x+delta)+(y+delta)*(y+delta));
	double theta2=abs(y+delta)/(y+delta)*acos((x+delta)/rho2);
	double rho3=sqrt((x-delta)*(x-delta)+(y+delta)*(y+delta));
	double theta3=abs(y+delta)/(y+delta)*acos((x-delta)/rho3);
	double rho4=sqrt((x+delta)*(x+delta)+(y-delta)*(y-delta));
	double theta4=abs(y-delta)/(y-delta)*acos((x+delta)/rho4);
	switch(l){
		case 0:
			value=+I*cexp(-I*S->d[0]*(kz(S,0,krho)-kz(S,1,krho)))/kz(S,1,krho)/Gamma(S,1,1,0,krho)/(1-Fres(S,1,1,0,krho)*Fres(S,1,1,2,krho)*cexp(2.0*I*kz(S,1,krho)*(S->d[0]-S->d[1])))*(Fres(S,1,1,2,krho)*(cexp(I*kz(S,1,krho)*(z+zs+h-2.0*S->d[1]))-cexp(I*kz(S,1,krho)*(z+zs-h-2.0*S->d[1])))-cexp(I*kz(S,1,krho)*(z+h-zs))+cexp(I*kz(S,1,krho)*(z-h-zs)));
			break;
		case 1:
			value=-I/kz(S,1,krho)/(1-Fres(S,1,1,0,krho)*Fres(S,1,1,2,krho)*cexp(2.0*I*kz(S,1,krho)*(S->d[0]-S->d[1])))*(Fres(S,1,1,0,krho)*Fres(S,1,1,2,krho)*(cexp(I*kz(S,1,krho)*(2.0*(S->d[0]-S->d[1])+z+h-zs))+cexp(I*kz(S,1,krho)*(2.0*(S->d[0]-S->d[1])-z+h+zs))-cexp(I*kz(S,1,krho)*(2.0*(S->d[0]-S->d[1])+z-h-zs))-cexp(I*kz(S,1,krho)*(2.0*(S->d[0]-S->d[1])-z+zs-h)))-Fres(S,1,1,2,krho)*(cexp(I*kz(S,1,krho)*(zs+z+h-2.0*S->d[1]))-cexp(I*kz(S,1,krho)*(zs+z-h-2.0*S->d[1])))-Fres(S,1,1,0,krho)*(cexp(I*kz(S,1,krho)*(-zs-z+h+2.0*S->d[0]))-cexp(I*kz(S,1,krho)*(-zs-z-h+2.0*S->d[0]))));
			value=value-I/kz(S,1,krho)*(z==zs? 2.0*(cexp(I*kz(S,1,krho)*h)-1):(cexp(I*kz(S,1,krho)*(abs(z-zs)+h))-cexp(I*kz(S,1,krho)*(abs(z-zs)-h))));
			break;
		case 2:
			value=+I*cexp(-I*S->d[1]*(kz(S,2,krho)-kz(S,1,krho)))/kz(S,1,krho)/Gamma(S,1,1,2,krho)/(1-Fres(S,1,1,0,krho)*Fres(S,1,1,2,krho)*cexp(2.0*I*kz(S,1,krho)*(S->d[0]-S->d[1])))*(Fres(S,0,1,0,krho)*(cexp(I*kz(S,1,krho)*(-z-zs+h+2.0*S->d[0]))-cexp(I*kz(S,1,krho)*(-z-zs-h+2.0*S->d[0])))-cexp(I*kz(S,1,krho)*(-z+h+zs))+cexp(I*kz(S,1,krho)*(-z-h+zs)));
			break;
	}
	return -2*PI*value*(besselJ(0,krho*rho1)+besselJ(0,krho*rho2)-besselJ(0,krho*rho3)-besselJ(0,krho*rho4))*kz(S,l-1,krho)/(S->k_0)/(S->k_0)/(S->Eps[l-1][1])/(S->Eps[l-1][0])/krho;
}

inline double complex Msxz(double complex krho, void * p){
	return 0;
}

double complex Matpxz(double complex krho, void *data){
	ScatParm * P; P = (ScatParm *) data;
	double complex value,angle_part;
	double h=P->h;
	double delta=P->delta;
	double z=P->z;
	double zs=P->zs;
	Sys *S; S=P->S;
	P->h=krho;
	int l=index_finder(S,z);
	angle_part=gauss_legendre_complex(100,angle2, P, 0.0, PI);
	switch(l){
		case 0:
			value=+I*cexp(-I*S->d[0]*(kz(S,0,krho)-kz(S,1,krho)))/kz(S,1,krho)/Gamma(S,2,1,0,krho)/(1-Fres(S,2,1,0,krho)*Fres(S,2,1,2,krho)*cexp(2.0*I*kz(S,1,krho)*(S->d[0]-S->d[1])))*(Fres(S,2,1,2,krho)*(cexp(I*kz(S,1,krho)*(z+zs+h-2.0*S->d[1]))-cexp(I*kz(S,1,krho)*(z+zs-h-2.0*S->d[1])))+cexp(I*kz(S,1,krho)*(z+h-zs))-cexp(I*kz(S,1,krho)*(z-h-zs)));
			break;
		case 1:
			value=-I/kz(S,1,krho)/(1-Fres(S,2,1,0,krho)*Fres(S,2,1,2,krho)*cexp(2.0*I*kz(S,1,krho)*(S->d[0]-S->d[1])))*(Fres(S,2,1,0,krho)*Fres(S,2,1,2,krho)*(-cexp(I*kz(S,1,krho)*(2.0*(S->d[0]-S->d[1])+z+h-zs))+cexp(I*kz(S,1,krho)*(2.0*(S->d[0]-S->d[1])-z+h+zs))+cexp(I*kz(S,1,krho)*(2.0*(S->d[0]-S->d[1])+z-h-zs))-cexp(I*kz(S,1,krho)*(2.0*(S->d[0]-S->d[1])-z+zs-h)))-Fres(S,2,1,2,krho)*(cexp(I*kz(S,1,krho)*(zs+z+h-2.0*S->d[1]))-cexp(I*kz(S,1,krho)*(zs+z-h-2.0*S->d[1])))+Fres(S,2,1,0,krho)*(cexp(I*kz(S,1,krho)*(-zs-z+h+2.0*S->d[0]))-cexp(I*kz(S,1,krho)*(-zs-z-h+2.0*S->d[0]))));
			value=value-I/kz(S,1,krho)*(z==zs? 0.0:(cexp(I*kz(S,1,krho)*(abs(z-zs)+h))-cexp(I*kz(S,1,krho)*(abs(z-zs)-h))))*(z>zs? -1.0:1.0);
			break;
		case 2:
			value=-I*cexp(-I*S->d[1]*(kz(S,2,krho)-kz(S,1,krho)))/kz(S,1,krho)/Gamma(S,2,1,2,krho)/(1-Fres(S,2,1,0,krho)*Fres(S,2,1,2,krho)*cexp(2.0*I*kz(S,1,krho)*(S->d[0]-S->d[1])))*(Fres(S,2,1,0,krho)*(cexp(I*kz(S,1,krho)*(-z-zs+h+2.0*S->d[0]))-cexp(I*kz(S,1,krho)*(-z-zs-h+2.0*S->d[0])))+cexp(I*kz(S,1,krho)*(-z+h+zs))-cexp(I*kz(S,1,krho)*(-z-h+zs)));
			break;
	}
	return -1.0*value*angle_part/(S->k_0)/(S->k_0)/(S->Eps[l-1][1])/(S->Eps[l-1][0]);
}

inline double complex Matsyx(double complex krho, void * data){
	return fsxy(krho, data);
}

inline double complex Matpyx(double complex krho, void * data){
	return fpxy(krho,data);
}

double complex Matsyy(double complex krho, void *data){
	ScatParm * P; P = (ScatParm *) data;
	double complex value,angle_part;
	double h=P->h;
	double delta=P->delta;
	double z=P->z;
	double zs=P->zs;
	Sys *S; S=P->S;
	P->h=krho;
	int l=index_finder(S,z);
	angle_part=gauss_legendre_complex(100,angle5, P, 0.0, PI);
	switch(l){
		case 0:
			value=-I*cexp(-I*S->d[0]*(kz(S,0,krho)-kz(S,1,krho)))/kz(S,1,krho)/Gamma(S,0,1,0,krho)/(1-Fres(S,0,1,0,krho)*Fres(S,0,1,2,krho)*cexp(2.0*I*kz(S,1,krho)*(S->d[0]-S->d[1])))*(Fres(S,0,1,2,krho)*(cexp(I*kz(S,1,krho)*(z+zs+h-2.0*S->d[1]))-cexp(I*kz(S,1,krho)*(z+zs-h-2.0*S->d[1])))+cexp(I*kz(S,1,krho)*(z+h-zs))-cexp(I*kz(S,1,krho)*(z-h-zs)));
			break;
		case 1:
			value=-I/kz(S,1,krho)/(1-Fres(S,0,1,0,krho)*Fres(S,0,1,2,krho)*cexp(2.0*I*kz(S,1,krho)*(S->d[0]-S->d[1])))*(Fres(S,0,1,0,krho)*Fres(S,0,1,2,krho)*(cexp(I*kz(S,1,krho)*(2.0*(S->d[0]-S->d[1])+z+h-zs))+cexp(I*kz(S,1,krho)*(2.0*(S->d[0]-S->d[1])-z+h+zs))-cexp(I*kz(S,1,krho)*(2.0*(S->d[0]-S->d[1])+z-h-zs))-cexp(I*kz(S,1,krho)*(2.0*(S->d[0]-S->d[1])-z+zs-h)))+Fres(S,0,1,2,krho)*(cexp(I*kz(S,1,krho)*(zs+z+h-2.0*S->d[1]))-cexp(I*kz(S,1,krho)*(zs+z-h-2.0*S->d[1])))+Fres(S,0,1,0,krho)*(cexp(I*kz(S,1,krho)*(-zs-z+h+2.0*S->d[0]))-cexp(I*kz(S,1,krho)*(-zs-z-h+2.0*S->d[0]))));
			value=value-I/kz(S,1,krho)*(z==zs? 2.0*(cexp(I*kz(S,1,krho)*h)-1):(cexp(I*kz(S,1,krho)*(abs(z-zs)+h))-cexp(I*kz(S,1,krho)*(abs(z-zs)-h))));
			break;
		case 2:
			value=-I*cexp(-I*S->d[1]*(kz(S,2,krho)-kz(S,1,krho)))/kz(S,1,krho)/Gamma(S,0,1,2,krho)/(1-Fres(S,0,1,0,krho)*Fres(S,0,1,2,krho)*cexp(2.0*I*kz(S,1,krho)*(S->d[0]-S->d[1])))*(Fres(S,0,1,0,krho)*(cexp(I*kz(S,1,krho)*(-z-zs+h+2.0*S->d[0]))-cexp(I*kz(S,1,krho)*(-z-zs-h+2.0*S->d[0])))+cexp(I*kz(S,1,krho)*(-z+h+zs))-cexp(I*kz(S,1,krho)*(-z-h+zs)));
			break;
	}
	return -1.0*value*angle_part/krho/kz(S,l-1,krho);
}

double complex Matpyy(double complex krho, void *data){
	ScatParm * P; P = (ScatParm *) data;
	double complex value,angle_part;
	double h=P->h;
	double delta=P->delta;
	double z=P->z;
	double zs=P->zs;
	Sys *S; S=P->S;
	P->h=krho;
	int l=index_finder(S,z);
	angle_part=gauss_legendre_complex(100,angle4, P, -PI*0.5, PI*0.5);
	switch(l){
		case 0:
			value=+I*cexp(-I*S->d[0]*(kz(S,0,krho)-kz(S,1,krho)))/kz(S,1,krho)/Gamma(S,1,1,0,krho)/(1-Fres(S,1,1,0,krho)*Fres(S,1,1,2,krho)*cexp(2.0*I*kz(S,1,krho)*(S->d[0]-S->d[1])))*(Fres(S,1,1,2,krho)*(cexp(I*kz(S,1,krho)*(z+zs+h-2.0*S->d[1]))-cexp(I*kz(S,1,krho)*(z+zs-h-2.0*S->d[1])))-cexp(I*kz(S,1,krho)*(z+h-zs))+cexp(I*kz(S,1,krho)*(z-h-zs)));
			break;
		case 1:
			value=-I/kz(S,1,krho)/(1-Fres(S,1,1,0,krho)*Fres(S,1,1,2,krho)*cexp(2.0*I*kz(S,1,krho)*(S->d[0]-S->d[1])))*(Fres(S,1,1,0,krho)*Fres(S,1,1,2,krho)*(cexp(I*kz(S,1,krho)*(2.0*(S->d[0]-S->d[1])+z+h-zs))+cexp(I*kz(S,1,krho)*(2.0*(S->d[0]-S->d[1])-z+h+zs))-cexp(I*kz(S,1,krho)*(2.0*(S->d[0]-S->d[1])+z-h-zs))-cexp(I*kz(S,1,krho)*(2.0*(S->d[0]-S->d[1])-z+zs-h)))-Fres(S,1,1,2,krho)*(cexp(I*kz(S,1,krho)*(zs+z+h-2.0*S->d[1]))-cexp(I*kz(S,1,krho)*(zs+z-h-2.0*S->d[1])))-Fres(S,1,1,0,krho)*(cexp(I*kz(S,1,krho)*(-zs-z+h+2.0*S->d[0]))-cexp(I*kz(S,1,krho)*(-zs-z-h+2.0*S->d[0]))));
			value=value-I/kz(S,1,krho)*(z==zs? 2.0*(cexp(I*kz(S,1,krho)*h)-1):(cexp(I*kz(S,1,krho)*(abs(z-zs)+h))-cexp(I*kz(S,1,krho)*(abs(z-zs)-h))));
			break;
		case 2:
			value=+I*cexp(-I*S->d[1]*(kz(S,2,krho)-kz(S,1,krho)))/kz(S,1,krho)/Gamma(S,1,1,2,krho)/(1-Fres(S,1,1,0,krho)*Fres(S,1,1,2,krho)*cexp(2.0*I*kz(S,1,krho)*(S->d[0]-S->d[1])))*(Fres(S,0,1,0,krho)*(cexp(I*kz(S,1,krho)*(-z-zs+h+2.0*S->d[0]))-cexp(I*kz(S,1,krho)*(-z-zs-h+2.0*S->d[0])))-cexp(I*kz(S,1,krho)*(-z+h+zs))+cexp(I*kz(S,1,krho)*(-z-h+zs)));
			break;
	}
	return -1.0*value*angle_part*kz(S,l-1,krho)/(S->k_0)/(S->k_0)/(S->Eps[l-1][1])/(S->Eps[l-1][0])/krho;
}

double complex Matsyz(double complex krho, void * data){
	return 0;
}

double complex Matszx(double complex krho, void * data){
	return 0;
}

double complex Matpzx(double complex krho, void *data){
	ScatParm * P; P = (ScatParm *) data;
	double complex value,angle_part;
	double h=P->h;
	double delta=P->delta;
	double z=P->z;
	double zs=P->zs;
	Sys *S; S=P->S;
	P->h=krho;
	int l=index_finder(S,z);
	angle_part=gauss_legendre_complex(100,angle2, P, 0.0, PI);
	switch(l){
		case 0:
			value=-I*cexp(-I*S->d[0]*(kz(S,0,krho)-kz(S,1,krho)))/kz(S,1,krho)/Gamma(S,1,1,0,krho)/(1-Fres(S,1,1,0,krho)*Fres(S,1,1,2,krho)*cexp(2.0*I*kz(S,1,krho)*(S->d[0]-S->d[1])))*(Fres(S,1,1,2,krho)*(cexp(I*kz(S,1,krho)*(z+zs+h-2.0*S->d[1]))-cexp(I*kz(S,1,krho)*(z+zs-h-2.0*S->d[1])))-cexp(I*kz(S,1,krho)*(z+h-zs))+cexp(I*kz(S,1,krho)*(z-h-zs)));
			break;
		case 1:
			value=-I/kz(S,1,krho)/(1-Fres(S,1,1,0,krho)*Fres(S,1,1,2,krho)*cexp(2.0*I*kz(S,1,krho)*(S->d[0]-S->d[1])))*(Fres(S,1,1,0,krho)*Fres(S,1,1,2,krho)*(-cexp(I*kz(S,1,krho)*(2.0*(S->d[0]-S->d[1])+z+h-zs))+cexp(I*kz(S,1,krho)*(2.0*(S->d[0]-S->d[1])-z+h+zs))+cexp(I*kz(S,1,krho)*(2.0*(S->d[0]-S->d[1])+z-h-zs))-cexp(I*kz(S,1,krho)*(2.0*(S->d[0]-S->d[1])-z+zs-h)))+Fres(S,1,1,2,krho)*(cexp(I*kz(S,1,krho)*(zs+z+h-2.0*S->d[1]))-cexp(I*kz(S,1,krho)*(zs+z-h-2.0*S->d[1])))-Fres(S,1,1,0,krho)*(cexp(I*kz(S,1,krho)*(-zs-z+h+2.0*S->d[0]))-cexp(I*kz(S,1,krho)*(-zs-z-h+2.0*S->d[0]))));
			value=value-I/kz(S,1,krho)*(z==zs? 0.0:(cexp(I*kz(S,1,krho)*(abs(z-zs)+h))-cexp(I*kz(S,1,krho)*(abs(z-zs)-h))))*(z>zs? -1.0:1.0);
			break;
		case 2:
			value=+I*cexp(-I*S->d[1]*(kz(S,2,krho)-kz(S,1,krho)))/kz(S,1,krho)/Gamma(S,1,1,2,krho)/(1-Fres(S,1,1,0,krho)*Fres(S,1,1,2,krho)*cexp(2.0*I*kz(S,1,krho)*(S->d[0]-S->d[1])))*(Fres(S,0,1,0,krho)*(cexp(I*kz(S,1,krho)*(-z-zs+h+2.0*S->d[0]))-cexp(I*kz(S,1,krho)*(-z-zs-h+2.0*S->d[0])))-cexp(I*kz(S,1,krho)*(-z+h+zs))+cexp(I*kz(S,1,krho)*(-z-h+zs)));
			break;
	}
	return -1.0*value*angle_part*krho/(S->k_0)/(S->k_0)/(S->Eps[l-1][1])/(S->Eps[l-1][0]);
}


double complex Matszy(double complex krho, void * data){
	return 0;
}


double complex Matpzy(double complex krho, void *data){
	ScatParm * P; P = (ScatParm *) data;
	double complex value,angle_part;
	double h=P->h;
	double delta=P->delta;
	double z=P->z;
	double zs=P->zs;
	Sys *S; S=P->S;
	P->h=krho;
	int l=index_finder(S,z);
	angle_part=gauss_legendre_complex(100,angle1, P, -0.5*PI, 0.5*PI);
	switch(l){
		case 0:
			value=-I*cexp(-I*S->d[0]*(kz(S,0,krho)-kz(S,1,krho)))/kz(S,1,krho)/Gamma(S,1,1,0,krho)/(1-Fres(S,1,1,0,krho)*Fres(S,1,1,2,krho)*cexp(2.0*I*kz(S,1,krho)*(S->d[0]-S->d[1])))*(Fres(S,1,1,2,krho)*(cexp(I*kz(S,1,krho)*(z+zs+h-2.0*S->d[1]))-cexp(I*kz(S,1,krho)*(z+zs-h-2.0*S->d[1])))-cexp(I*kz(S,1,krho)*(z+h-zs))+cexp(I*kz(S,1,krho)*(z-h-zs)));
			break;
		case 1:
			value=-I/kz(S,1,krho)/(1-Fres(S,1,1,0,krho)*Fres(S,1,1,2,krho)*cexp(2.0*I*kz(S,1,krho)*(S->d[0]-S->d[1])))*(Fres(S,1,1,0,krho)*Fres(S,1,1,2,krho)*(-cexp(I*kz(S,1,krho)*(2.0*(S->d[0]-S->d[1])+z+h-zs))+cexp(I*kz(S,1,krho)*(2.0*(S->d[0]-S->d[1])-z+h+zs))+cexp(I*kz(S,1,krho)*(2.0*(S->d[0]-S->d[1])+z-h-zs))-cexp(I*kz(S,1,krho)*(2.0*(S->d[0]-S->d[1])-z+zs-h)))+Fres(S,1,1,2,krho)*(cexp(I*kz(S,1,krho)*(zs+z+h-2.0*S->d[1]))-cexp(I*kz(S,1,krho)*(zs+z-h-2.0*S->d[1])))-Fres(S,1,1,0,krho)*(cexp(I*kz(S,1,krho)*(-zs-z+h+2.0*S->d[0]))-cexp(I*kz(S,1,krho)*(-zs-z-h+2.0*S->d[0]))));
			value=value-I/kz(S,1,krho)*(z==zs? 0.0:(cexp(I*kz(S,1,krho)*(abs(z-zs)+h))-cexp(I*kz(S,1,krho)*(abs(z-zs)-h))))*(z>zs? -1.0:1.0);
			break;
		case 2:
			value=+I*cexp(-I*S->d[1]*(kz(S,2,krho)-kz(S,1,krho)))/kz(S,1,krho)/Gamma(S,1,1,2,krho)/(1-Fres(S,1,1,0,krho)*Fres(S,1,1,2,krho)*cexp(2.0*I*kz(S,1,krho)*(S->d[0]-S->d[1])))*(Fres(S,0,1,0,krho)*(cexp(I*kz(S,1,krho)*(-z-zs+h+2.0*S->d[0]))-cexp(I*kz(S,1,krho)*(-z-zs-h+2.0*S->d[0])))-cexp(I*kz(S,1,krho)*(-z+h+zs))+cexp(I*kz(S,1,krho)*(-z-h+zs)));
			break;
	}
	return -1.0*value*angle_part*krho/(S->k_0)/(S->k_0)/(S->Eps[l-1][1])/(S->Eps[l-1][0]);
}

double complex Matszz(double complex krho, void * data){
	return 0;
}

double complex Matpzz(double complex krho, void *data){
	ScatParm * P; P = (ScatParm *) data;
	double complex value,angle_part;
	double h=P->h;
	double delta=P->delta;
	double z=P->z;
	double zs=P->zs;
	Sys *S; S=P->S;
	P->h=krho;
	int l=index_finder(S,z);
	angle_part=gauss_legendre_complex(100,angle3, P, 0.0, 0.5*PI);
	switch(l){
		case 0:
			value=-I*cexp(-I*S->d[0]*(kz(S,0,krho)-kz(S,1,krho)))/kz(S,1,krho)/Gamma(S,2,1,0,krho)/(1-Fres(S,2,1,0,krho)*Fres(S,2,1,2,krho)*cexp(2.0*I*kz(S,1,krho)*(S->d[0]-S->d[1])))*(Fres(S,2,1,2,krho)*(cexp(I*kz(S,1,krho)*(z+zs+h-2.0*S->d[1]))-cexp(I*kz(S,1,krho)*(z+zs-h-2.0*S->d[1])))+cexp(I*kz(S,1,krho)*(z+h-zs))-cexp(I*kz(S,1,krho)*(z-h-zs)));
			break;
		case 1:
			value=-I/kz(S,1,krho)/(1-Fres(S,2,1,0,krho)*Fres(S,2,1,2,krho)*cexp(2.0*I*kz(S,1,krho)*(S->d[0]-S->d[1])))*(Fres(S,2,1,0,krho)*Fres(S,2,1,2,krho)*(cexp(I*kz(S,1,krho)*(2.0*(S->d[0]-S->d[1])+z+h-zs))+cexp(I*kz(S,1,krho)*(2.0*(S->d[0]-S->d[1])-z+h+zs))-cexp(I*kz(S,1,krho)*(2.0*(S->d[0]-S->d[1])+z-h-zs))-cexp(I*kz(S,1,krho)*(2.0*(S->d[0]-S->d[1])-z+zs-h)))+Fres(S,2,1,2,krho)*(cexp(I*kz(S,1,krho)*(zs+z+h-2.0*S->d[1]))-cexp(I*kz(S,1,krho)*(zs+z-h-2.0*S->d[1])))+Fres(S,2,1,0,krho)*(cexp(I*kz(S,1,krho)*(-zs-z+h+2.0*S->d[0]))-cexp(I*kz(S,1,krho)*(-zs-z-h+2.0*S->d[0]))));
			value=value-I/kz(S,1,krho)*(z==zs? 2.0*(cexp(I*kz(S,1,krho)*h)-1):(cexp(I*kz(S,1,krho)*(abs(z-zs)+h))-cexp(I*kz(S,1,krho)*(abs(z-zs)-h))));
			break;
		case 2:
			value=-I*cexp(-I*S->d[1]*(kz(S,2,krho)-kz(S,1,krho)))/kz(S,1,krho)/Gamma(S,2,1,2,krho)/(1-Fres(S,2,1,0,krho)*Fres(S,2,1,2,krho)*cexp(2.0*I*kz(S,1,krho)*(S->d[0]-S->d[1])))*(Fres(S,2,1,0,krho)*(cexp(I*kz(S,1,krho)*(-z-zs+h+2.0*S->d[0]))-cexp(I*kz(S,1,krho)*(-z-zs-h+2.0*S->d[0])))+cexp(I*kz(S,1,krho)*(-z+h+zs))-cexp(I*kz(S,1,krho)*(-z-h+zs)));
			break;
	}
	return -1.0*value*angle_part/(S->k_0)/(S->k_0)/((S->Eps[l-1][1])*kz(S,l-1,krho))/(S->Eps[l-1][0])*krho;
}
