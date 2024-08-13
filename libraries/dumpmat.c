/******************  dumpmat.c  ******************************
*									*
*  void dumpmat(matrix *mat )          					*
*  print out a 4x4 complex matrix					*
*/
#include "../include/config.h"
#include <stdio.h>
#include "../include/complex.h"
#include "../include/sp.h"

void dumpmat(matrix *m ){
int i,j;
    for(i=0;i<DIMF;i++){
	for(j=0;j<DIMF;j++)printf("(%.6g,%.6g)\t",
	    m->e[i][j].real,m->e[i][j].imag);
    printf("\n");	
    }
    printf("\n");
}
