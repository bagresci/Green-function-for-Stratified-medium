#include <stdio.h>
#include <stdlib.h>
#include "math.h"
#include "mkj.h"
#ifndef PI
	#define PI 3.1415926535897932384626433832795028841971693993751
#endif
/*
void main(){
	FILE *fxx,*fzx,*fzz;
	fxx=fopen("metalfxx.txt","w");fzz=fopen("metalfzz.txt","w");fzx=fopen("metalfzx.txt","w");
	double depth[2]={50.0,0.0};
	Sys *S; S=createSys((parm){.N=3,.d=depth, .k0=2*PI/500.0});
	S->Eps[0][0]=1.0;S->Eps[0][1]=1.96;S->Eps[0][2]=1.96;
	S->Eps[1][0]=1.0;S->Eps[1][1]=-15.99+0.8*I;S->Eps[1][2]=-15.99+0.8*I;
	S->Eps[2][0]=1.0;S->Eps[2][1]=1.0;S->Eps[2][2]=1.0;
	double complex Gxx,Gzz,Gs;
	for(int i=0;i<10000;i++){
		S->k_0=2*PI/(500.0+10.0*(i/200));
		Gxx=0.25*I/PI*integral_imag(S,0.0,0.0,-75.0+1.0*(i%200),-75.0+1.0*(i%200),fpxx);
		Gzz=0.25*I/PI*integral_imag(S,0.0,0.0,-75.0+1.0*(i%200),-75.0+1.0*(i%200),fpzz);
		Gs=0.25*I/PI*integral_imag(S,0.0,0.0,-75.0+1.0*(i%200),-75.0+1.0*(i%200),fsxx);
		fprintf(fxx,"%d,%d,%.15f\n",-75+(i%200),-75+(i%200),cimag(Gxx));
		fprintf(fzz,"%d,%d,%.15f\n",-75+(i%200),-75+(i%200),cimag(Gzz));
		fprintf(fzx,"%d,%d,%.15f\n",-75+(i%200),-75+(i%200),cimag(Gs));
	}
	fclose(fxx);
	fclose(fzz);
	fclose(fzx);
}*/
/*
void main(){
	double depth[2]={50.0,0.0};
	Sys *S; S=createSys((parm){.N=3,.d=depth, .k0=2*PI/500.0});
	S->Eps[0][0]=1.0;S->Eps[0][1]=1.4;S->Eps[0][2]=1.4;
	S->Eps[1][0]=1.0;S->Eps[1][1]=0.1+4.0*I;S->Eps[1][2]=0.1+4.0*I;
	S->Eps[2][0]=1.0;S->Eps[2][1]=1.0;S->Eps[2][2]=1.0;
	double complex Gxx,Gzz,Gs;
	for(int i=0;i<12;i++){
		Gxx=0.25*I/PI*integral_real(S,0.0,0.0,-10.0+i*5.0,-10.0+i*5.0,fpxx);
		Gzz=0.25*I/PI*integral_real(S,0.0,0.0,-10.0+i*5.0,-10.0+i*5.0,fpzz);
		Gs=0.25*I/PI*integral_real(S,0.0,0.0,-10.0+i*5.0,-10.0+i*5.0,fsxx);
		printf("%d,%.15f\n",-10+i*5,cimag(Gxx));
		printf("%d,%.15f\n",-10+i*5,cimag(Gzz));
		printf("%d,%.15f\n",-10+i*5,cimag(Gs));
	}
}
*/
/*
void main(){
	double depth[2]={50.0,0.0};
	Sys *S; S=createSys((parm){.N=3,.d=depth, .k0=2*PI/500.0});
	S->Eps[0][0]=1.0;S->Eps[0][1]=1.96;S->Eps[0][2]=1.96;
	S->Eps[1][0]=1.0;S->Eps[1][1]=0.1+4.0*I;S->Eps[1][2]=0.1+4.0*I;
	S->Eps[2][0]=1.0;S->Eps[2][1]=1.0;S->Eps[2][2]=1.0;
	ItgParm P;
	P.func_num=1;P.rho=0;P.phi=0;P.z=55.0;P.zs=55.0;P.S=S;P.f=fpxx;
	double complex value1,value2;
	for(int i=0;i<100;i++){
		value1=fsxx(2.4*S->k_0*(i*0.1+100.0),&P);
		value2=fpzz(2.4*S->k_0*(i*0.1+100.0),&P);
		printf("%f,%.15f\n",2.4*S->k_0*(i*0.1+1000),creal(value1));
		printf("%f,%.15f\n",2.4*S->k_0*(i*0.1+1000),creal(value2));
	}
	

}*/
/*
void main(){
	FILE *f;
	f=fopen("finding_plasmony.txt","w");
	double depth[2]={50.0,0.0};
	Sys *S; S=createSys((parm){.N=3,.d=depth, .k0=2*PI/370.0});
	S->Eps[0][0]=1.0;S->Eps[0][1]=1.96;S->Eps[0][2]=1.96;
	S->Eps[1][0]=1.0;S->Eps[1][1]=-2.475+0.5054*I;S->Eps[1][2]=-2.475+0.5054*I;
	S->Eps[2][0]=1.0;S->Eps[2][1]=1.0;S->Eps[2][2]=1.0;
	double zs=55.0;
	double rho;double z;
	double complex v1,v2,v3,v4,v5,v6,w1,w2,w3,w4,w5,w6;
	double Esquare1,Esquare2,Esquare3,Esquare4;
	for(int i=0;i<120000;i++){
		rho=10.0+5.0*(i/600);
		z=-300.0+1.0*(i%600);
		v1=0.25*I/PI*(integral_imag(S,rho,0.0,z,zs,fsxx)+integral_imag(S,rho,0.0,z,zs,fpxx));
		v2=0.25*I/PI*integral_imag(S,rho,0.0,z,zs,fpzx);
		v3=0.25*I/PI*(integral_imag(S,rho,0.0,z,zs,fsyx)+integral_imag(S,rho,0.0,z,zs,fpyx));
		v4=0.25*I/PI*(integral_imag(S,rho,PI,z,zs,fsxx)+integral_imag(S,rho,PI,z,zs,fpxx));
		v5=0.25*I/PI*integral_imag(S,rho,PI,z,zs,fpzx);
		v6=0.25*I/PI*(integral_imag(S,rho,PI,z,zs,fsyx)+integral_imag(S,rho,PI,z,zs,fpyx));
		w1=0.25*I/PI*(integral_imag(S,rho+1000.0,0.0,z,zs,fsxx)+integral_imag(S,rho+1000.0,0.0,z,zs,fpxx));
		w2=0.25*I/PI*integral_imag(S,rho+1000.0,0.0,z,zs,fpzx);
		w3=0.25*I/PI*(integral_imag(S,rho+1000.0,0.0,z,zs,fsyx)+integral_imag(S,rho+1000.0,0.0,z,zs,fpyx));	
		w4=0.25*I/PI*(integral_imag(S,rho+1000.0,PI,z,zs,fsxx)+integral_imag(S,rho+1000.0,PI,z,zs,fpxx));
		w5=0.25*I/PI*integral_imag(S,rho+1000.0,PI,z,zs,fpzx);
		w6=0.25*I/PI*(integral_imag(S,rho+1000.0,PI,z,zs,fsyx)+integral_imag(S,rho+1000.0,PI,z,zs,fpyx));	

		Esquare1=creal(v1)*creal(v1)+creal(v2)*creal(v2)+creal(v3)*creal(v3);
		Esquare2=creal(w1)*creal(w1)+creal(w2)*creal(w2)+creal(w3)*creal(w3);

		Esquare3=creal(v4)*creal(v4)+creal(v5)*creal(v5)+creal(v6)*creal(v6);
		Esquare4=creal(w4)*creal(w4)+creal(w5)*creal(w5)+creal(w6)*creal(w6);
		fprintf(f,"%f,%f,%.15f\n%f,%f,%.15f\n",rho,z,Esquare1,rho+1000.0,z,Esquare2);
		fprintf(f,"%f,%f,%.15f\n%f,%f,%.15f\n",-1.0*rho,z,Esquare3,-1.0*rho-1000.0,z,Esquare4);
	
	}
	fclose(f);

}*/
/*
void main(){
	FILE *f;
	f=fopen("finding_plasmon2d.txt","w");
	double depth[2]={50.0,0.0};
	Sys *S; S=createSys((parm){.N=3,.d=depth, .k0=2*PI/370.0});
	S->Eps[0][0]=1.0;S->Eps[0][1]=1.96;S->Eps[0][2]=1.96;
	S->Eps[1][0]=1.0;S->Eps[1][1]=-2.475+0.5054*I;S->Eps[1][2]=-2.475+0.5054*I;
	S->Eps[2][0]=1.0;S->Eps[2][1]=1.0;S->Eps[2][2]=1.0;
	double zs=55.0;
	double rho;double z=50.0;
	double complex v1,v2,v3,v4;
	for(int i=0;i<6500;i++){
		rho=0.1*i;
		//rho=2.0*(i/150);
		//z=-50.0+(i%150);
		v1=0.25*I/PI*(integral_imag(S,rho,0.0,z,zs,fsxx)+integral_imag(S,rho,0.0,z,zs,fpxx));
		//v2=0.25*I/PI*integral_imag(S,rho,0.0,z,zs,fpzx);
		v3=0.25*I/PI*(integral_imag(S,rho,PI,z,zs,fsxx)+integral_imag(S,rho,PI,z,zs,fpxx));
		//v4=0.25*I/PI*integral_imag(S,rho,PI,z,zs,fpzx);
		fprintf(f,"%f,%.15f,%.15f\n",rho,creal(v1),cimag(v1));
		fprintf(f,"%f,%.15f,%.15f\n",-1.0*rho,creal(v3),cimag(v3));
	}
	fclose(f);

}*/

/*
void main(){
	double complex mScatter[475][475];
	FILE *f;
	f=fopen("testing.txt","w");
	double complex v1,v2,v3;
	double depth[2]={10.0,0.0};
	Sys *S; S=createSys((parm){.N=3,.d=depth, .k0=2*PI/370.0});
	S->Eps[0][0]=1.0;S->Eps[0][1]=1.0;S->Eps[0][2]=1.0;
	S->Eps[1][0]=1.0;S->Eps[1][1]=-2.475+0.5054*I;S->Eps[1][2]=-2.475+0.5054*I;
	S->Eps[2][0]=1.0;S->Eps[2][1]=1.0;S->Eps[2][2]=1.0;
	
	for(int i=0;i<125;i++){
		for(int j=0;j<125;j++){
			v1=(i==j? 1.0:0.0)-(S->k_0)*(S->k_0)*(integral_real(S, 2.0*(i/25)-1.0,2.0*(j/25)-1.0,2.0*(i%5)-1.0,2.0*(j%5)-1.0,2.0*((i/5)%5)-1.0,2.0*((j/5)%5)-1.0,1.0,1.0,Msxx)+integral_real(S,  2.0*(i/25)-1.0,2.0*(j/25)-1.0,2.0*(i%5)-1.0,2.0*(j%5)-1.0,2.0*((i/5)%5)-1.0,2.0*((j/5)%5)-1.0,1.0,1.0,Mpxx));
			v2=-(S->k_0)*(S->k_0)*(integral_real(S,  2.0*(i/25)-1.0,2.0*(j/25)-1.0,2.0*(i%5)-1.0,2.0*(j%5)-1.0,2.0*((i/5)%5)-1.0,2.0*((j/5)%5)-1.0,1.0,1.0,Msxy)+integral_real(S, 2.0*(i/25)-1.0,2.0*(j/25)-1.0,2.0*(i%5)-1.0,2.0*(j%5)-1.0,2.0*((i/5)%5)-1.0,2.0*((j/5)%5)-1.0,1.0,1.0,Mpxy));
			v3=-(S->k_0)*(S->k_0)*integral_real(S, 2.0*(i/25)-1.0,2.0*(j/25)-1.0,2.0*(i%5)-1.0,2.0*(j%5)-1.0,2.0*((i/5)%5)-1.0,2.0*((j/5)%5)-1.0,9.0,9.0,1.0,1.0,Mpxz);
			mScatter[3*i][3*j]=v1;mScatter[3*i][3*j+1]=v2;mScatter[3*i][3*j+2]=v3;
			printf("%f+I%f\n",creal(v1),cimag(v1));
			printf("%f+I%f\n",creal(v2),cimag(v2));
			printf("%f+I%f\n",creal(v3),cimag(v3));
		}
		for(int j=0;j<125;j++){
			v1=-(S->k_0)*(S->k_0)*(integral_real(S, 2.0*(i/25)-1.0,2.0*(j/25)-1.0,2.0*(i%5)-1.0,2.0*(j%5)-1.0,2.0*((i/5)%5)-1.0,2.0*((j/5)%5)-1.0,1.0,1.0,Msyx)+integral_real(S,  2.0*(i/25)-1.0,2.0*(j/25)-1.0,2.0*(i%5)-1.0,2.0*(j%5)-1.0,2.0*((i/5)%5)-1.0,2.0*((j/5)%5)-1.0,1.0,1.0,Mpyx));
			v2=(i==j? 1.0:0.0)-(S->k_0)*(S->k_0)*(integral_real(S,  2.0*(i/25)-1.0,2.0*(j/25)-1.0,2.0*(i%5)-1.0,2.0*(j%5)-1.0,2.0*((i/5)%5)-1.0,2.0*((j/5)%5)-1.0,1.0,1.0,Msyy)+integral_real(S, 2.0*(i/25)-1.0,2.0*(j/25)-1.0,2.0*(i%5)-1.0,2.0*(j%5)-1.0,2.0*((i/5)%5)-1.0,2.0*((j/5)%5)-1.0,1.0,1.0,Mpyy));
			v3=-(S->k_0)*(S->k_0)*integral_real(S, 2.0*(i/25)-1.0,2.0*(j/25)-1.0,2.0*(i%5)-1.0,2.0*(j%5)-1.0,2.0*((i/5)%5)-1.0,2.0*((j/5)%5)-1.0,9.0,9.0,1.0,1.0,Mpyz);
			mScatter[3*i+1][3*j]=v1;mScatter[3*i+1][3*j+1]=v2;mScatter[3*i+1][3*j+2]=v3;
			printf("%f+I%f\n",creal(v1),cimag(v1));
			printf("%f+I%f\n",creal(v2),cimag(v2));
			printf("%f+I%f\n",creal(v3),cimag(v3));
		}
		for(int j=0;j<125;j++){
			v1=-(S->k_0)*(S->k_0)*integral_real(S,  2.0*(i/25)-1.0,2.0*(j/25)-1.0,2.0*(i%5)-1.0,2.0*(j%5)-1.0,2.0*((i/5)%5)-1.0,2.0*((j/5)%5)-1.0,1.0,1.0,Mpzx);
			v2=-(S->k_0)*(S->k_0)*integral_real(S, 2.0*(i/25)-1.0,2.0*(j/25)-1.0,2.0*(i%5)-1.0,2.0*(j%5)-1.0,2.0*((i/5)%5)-1.0,2.0*((j/5)%5)-1.0,1.0,1.0,Mpzy);
			v3=(i==j? 1.0:0.0)*(1.0+8.0/(S->k_0)/(S->k_0)/;S->Eps[1][2])-(S->k_0)*(S->k_0)*integral_real(S, 2.0*(i/25)-1.0,2.0*(j/25)-1.0,2.0*(i%5)-1.0,2.0*(j%5)-1.0,2.0*((i/5)%5)-1.0,2.0*((j/5)%5)-1.0,9.0,9.0,1.0,1.0,Mpzz);
			mScatter[3*i+2][3*j]=v1;mScatter[3*i+2][3*j+1]=v2;mScatter[3*i+2][3*j+2]=v3;
			printf("%f+I%f\n",creal(v1),cimag(v1));
			printf("%f+I%f\n",creal(v2),cimag(v2));
			printf("%f+I%f\n",creal(v3),cimag(v3));
		}

	}
	fclose(f);
}
*/

void main(){
	FILE *f,*g;
	f=fopen("metalfxx.txt","w");
	g=fopen("metalfzz.txt","w");
	double complex v1,v2;
	double depth[2]={10.0,0.0};
	Sys *S; S=createSys((parm){.N=3,.d=depth, .k0=2*PI/633.0});
	S->Eps[0][0]=1.0;S->Eps[0][1]=1.0;S->Eps[0][2]=1.0;
	//S->Eps[1][0]=1.0;S->Eps[1][1]=-2.475+0.5054*I;S->Eps[1][2]=-2.475+0.5054*I;
	S->Eps[1][0]=1.0;S->Eps[1][1]=1.0;S->Eps[1][2]=1.0;
	S->Eps[2][0]=1.0;S->Eps[2][1]=1.0;S->Eps[2][2]=1.0;
	//for(int i=0;i<5;i++){
	//	v1=integral_real(S,1.0,1.0,1.0,1.0,3.0,3.0,1.0-0.2*i,1.0-0.2*i,Matsxx)+integral_real(S,1.0,1.0,1.0,1.0,3.0,3.0,1.0-0.2*i,1.0-0.2*i,Matpxx);
	//	printf("%.15f,%.15f\n",creal(v1),cimag(v1));
	//}
	//v1=integral_real(S,1.0,1.0,1.0,1.0,4.0,3.0,1.0,0.5,Matsxx)+integral_real(S,1.0,1.0,1.0,1.0,2.0,3.0,1.0,0.5,Matsxx)+integral_real(S,1.0,1.0,1.0,1.0,4.0,3.0,1.0,0.5,Matpxx)+integral_real(S,1.0,1.0,1.0,1.0,2.0,3.0,1.0,0.5,Matpxx);
	v1=(integral_real(S,0.0,0.0,1.0,1.0,0.0000000000001,0.0000000000001,Matsxx)+integral_real(S,0.0,0.0,1.0,1.0,0.0000000000001,0.0000000000001,Matpxx))/S->k_0/S->k_0;
	printf("%f,%f\n",creal(v1),cimag(v1));
/*
	for(int i=0;i<10;i++){
		v1=integral_real(S,0.0,0.0,1.0,1.0,1.0/pow(10.0,i),1.0/pow(10.0,i),Matsxx)+integral_real(S,0.0,0.0,1.0,1.0,1.0/pow(10.0,i),1.0/pow(10.0,i),Matpxx);
		v2=integral_real(S,0.0,0.0,1.0,1.0,10./pow(10.0,i),1.0/pow(10.0,i),Matsyy)+integral_real(S,0.0,0.0,1.0,1.0,1.0/pow(10.0,i),1.0/pow(10.0,i),Matpyy);
		fprintf(f,"%d,%.15f,%.15f\n",i,creal(v1),cimag(v1));
		fprintf(g,"%d,%.15f,%.15f\n",i,creal(v2),cimag(v2));	
	}
	*/
	fclose(f);
	fclose(g);
}