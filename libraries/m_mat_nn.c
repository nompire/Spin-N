/*******************  m_mat_nn.c             ****************************
*									*
*                                        				*
* matrix multiply, no adjoints 						*
* C  <-  A*B								*
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/sp.h"

//#ifndef FAST
void mult_mat_nn(matrix *a, matrix *b, matrix *c ){
register int i,j,k;
register complex x,y;
    for(i=0;i<DIMF;i++)for(j=0;j<DIMF;j++){
	x.real=x.imag=0.0;
	for(k=0;k<DIMF;k++){
	    CMUL( a->e[i][k] , b->e[k][j] , y );
	    CSUM( x , y );
	}
	c->e[i][j] = x;
    }
}


