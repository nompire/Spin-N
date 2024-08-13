/***************** conjg.c                     **************************
*									*
* B  <- A*,  adjoint of an su2 matrix 				*
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/sp.h"

/* conjugate of an  matrix */
void conjug(matrix *a,matrix *b ){
register int i,j;
    for(i=0;i<DIMF;i++){
    for(j=0;j<DIMF;j++){
	CONJG( a->e[i][j], b->e[i][j] );
    }
  }
}
