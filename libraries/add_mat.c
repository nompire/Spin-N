/********************  add_mat.c            *****************************
*									*
*  Add two matrices a += b 						*
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/sp.h"

void add_mat( matrix *a, matrix *b) {
register int i,j;
    for(i=0;i<DIMF;i++)for(j=0;j<DIMF;j++){
	CSUM( a->e[i][j], b->e[i][j]);
    }
}
