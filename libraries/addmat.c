/********************  addmat.c (in su3.a)  *****************************
*									*
*  Add two  matrices 						*
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/sp.h"

void add_matrix( matrix *a, matrix *b, matrix *c ) {
register int i,j;
    for(i=0;i<DIMF;i++)for(j=0;j<DIMF;j++){
	CADD( a->e[i][j], b->e[i][j], c->e[i][j] );
    }
}
