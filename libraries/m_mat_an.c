/******************  m_mat_an.c *****************************
*									*
* matrix multiply, first matrix is adjoint 				*
* C <-  A_adjoint*B							*
*/


#include "../include/config.h"
#include "../include/complex.h"
#include "../include/sp.h"

void mult_mat_an( matrix *a, matrix *b, matrix *c ){
register int i,j,k;
register complex x,y;
    for(i=0;i<DIMF;i++)for(j=0;j<DIMF;j++){
	x.real=x.imag=0.0;
	for(k=0;k<DIMF;k++){
	    CMULJ_( a->e[k][i] , b->e[k][j], y );
	    CSUM( x , y );
	}
	c->e[i][j] = x;
    }
}


