/****************  m_mat_na.c             *******************************
*									*
* matrix multiply, second matrix is adjoint 				*
* C  <-  A*B_adjoint							*
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/sp.h"

void mult_mat_na(matrix *a,matrix *b,matrix *c ){
register int i,j,k;
register complex x,y;
    for(i=0;i<DIMF;i++)for(j=0;j<DIMF;j++){
	x.real=x.imag=0.0;
	for(k=0;k<DIMF;k++){
	    CMUL_J( a->e[i][k] , b->e[j][k] , y );
	    CSUM( x , y );
	}
	c->e[i][j] = x;
    }
}


