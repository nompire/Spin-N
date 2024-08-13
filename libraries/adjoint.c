/******************  su2_adjoint.c  (in su2.a) **************************
*									*
* B  <- A_adjoint,  adjoint of an su2 matrix 				*
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/sp.h"

/* adjoint of an  matrix */
void adjoint(matrix *a,matrix *b ){
register int i,j;
    for(i=0;i<DIMF;i++)for(j=0;j<DIMF;j++){
	CONJG( a->e[j][i], b->e[i][j] );
    }
}
