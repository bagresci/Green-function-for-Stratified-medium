#include "mkj.h"

inline double complex besselJ(int order, double complex z)
{
    // Input values for Fortran subroutines.
    double zr = creal(z);
    double zi = cimag(z);
    double nu = fabs((double) order);
    int kode = 1;
    int N = 1;

    // Output values.
    double cyr,cyi;
    int nz,ierr;

    // External function call
    zbesj_wrap(zr,zi,nu,kode,N,&cyr,&cyi,&nz,&ierr);    // Call Fortran subroutine.
    double complex answer = cyr + _Complex_I * cyi;               // Placeholder for output.

    // If order is negative, then we must apply the reflection formula.
    if (order < 0)
    {
        answer *= (order & 1 ? -1.0 : 1.0);
    }

    // If the return code is not normal, we print the error code.
    if (ierr!=0){
    	printf("besselJ: Error code %d\n",ierr) ;
    } 

    return answer;
}

inline double complex besselY(int order, double complex z)
{
    // Input values for Fortran subroutines
    double zr = creal(z);
    double zi = cimag(z);
    double nu = fabs((double) order);
    int kode = 1;
    int N = 1;

    // Output and temporary varibles
    double cyr, cyi,cwrkr,cwrki;
    int nz, ierr;

    // External function call
    zbesy_wrap(zr,zi,nu,kode,N,&cyr,&cyi,&nz,&cwrkr,&cwrki,&ierr); // Call Fortran subroutine.
    
    // In passing from C++ to FORTRAN, the exact zero becomes the numerical zero (10^(-14)). 
    // The limiting form of Y_nu(z) for high order, -Gamma(nu)/pi*Re(z)^(-nu)*(1-i*nu*Im(z)/Re(z)),
    // leads to product of the form zero*infinity, which destroys numerical precision. We hence
    // manually set the imaginary part of the answer to zero is the imaginary part of the input
    // is zero. 
    if (zi == 0.0) cyi=0.0;
    double complex answer=cyr+_Complex_I *cyi;                           // Placeholder for output

    // If order is negative, we must apply the reflection formula.
    if (order < 0)
    {
        answer *= (order & 1 ? -1.0 : 1.0);
    }

    // If the return code is not normal, we print the error code.
    if (ierr!=0){
        printf("besselY: Error code %d\n",ierr);
    } 
    return answer;
}

inline double complex besselH1(int order, double complex z)
{
    // Input values.
    double zr = creal(z);
    double zi = cimag(z);
    double nu = fabs((double) order);
    int kode = 1;
    int N = 1;
    int kind = 1;

//    // Output values
    double cyr,cyi;
    int nz,ierr;

//    // External function call.
    zbesh_wrap(zr,zi,nu,kode,kind,N,&cyr,&cyi,&nz,&ierr);
    return cyr + _Complex_I*cyi;

//    // Reflection formula if order is negative.
//    if (order < 0.0)
//    {
//        // Compute complex exponential.
//        std::complex<double> i(0.0,1.0);
//        i = std::exp(arma::datum::pi*nu*i);
//        answer *= i;
//    }

 ///   double complex i(0.0,1.0);
//    return besselJ(order,z)+_Complex_I*besselY(order,z);
}

inline double complex besselH2(int order, double complex z){
//    // Input values.
    double zr = creal(z);
    double zi = cimag(z);
    double nu = fabs((double) order);
    int kode = 1;
    int N = 1;
    int kind = 2;

//    // Output values
    double cyr,cyi;
    int nz,ierr;

//    // External function call.
    zbesh_wrap(zr,zi,nu,kode,kind,N,&cyr,&cyi,&nz,&ierr);
    return cyr+_Complex_I*cyi;

//    // Reflection formula if order is negative.
//    if (order < 0.0 )
//    {
//        std::complex<double> i(0.0,1.0);
//        i = std::exp(-arma::datum::pi*nu*i);
//        answer *= i;
//    }

///  std::complex<double> i(0.0,1.0);
//    return besselJ(order,z)-_Complex_I*besselY(order,z);
}