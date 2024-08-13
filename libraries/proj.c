/*****************   proj.c ******************************
*									*
* void su2_projector( su2_vector *a, su2_vector *b, su2_matrix *c )	*
* C  <- outer product of A and B					*
*  C_ij = A_i * B_adjoint_j						*
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/sp.h"


void proj(vector *a,vector *b, matrix *c ){
register int i,j;

    for(i=0;i<DIMF;i++)for(j=0;j<DIMF;j++){
	CMUL_SUM(a->c[i], b->c[j], c->e[i][j]);
    }
}


