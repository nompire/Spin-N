/********************  sub_mat.c            *****************************
*									*
*  Subtract two  matrices 						*
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/sp.h"

void sub_mat(matrix *a, matrix *b) {
register int i,j;
    for(i=0;i<DIMF;i++)for(j=0;j<DIMF;j++){
	CDIF( a->e[i][j], b->e[i][j]);
    }
}
