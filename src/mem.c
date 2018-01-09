#include "mkj.h"

double complex ** makefield_doublecomplex(int row, int col){
	double complex ** state; void * pvoid;
    if(NULL==(pvoid = calloc(row,sizeof(double complex *)))){
        printf("Fail : creating makefield_double \n");
        exit(1); //return EXIT_FAILURE;
    }else state = (double complex **)pvoid;
    if(NULL==(pvoid =calloc(col*row,sizeof(double complex)))){
        printf("Fail : (bulk) creating states \n");
        exit(1); //return EXIT_FAILURE;
    } else state[0] = (double complex *)pvoid;
    for( int i=1;i<row;i++) state[i] = state[0] + col*i;
    return state;
}