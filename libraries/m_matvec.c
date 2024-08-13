/****************  m_matvec.c  (in su2.a) *******************************
*									*
* matrix times vector multiply, no adjoints 				*
*  C  <-  A*B								*
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/sp.h"


void mult_mat_vec(matrix *a,vector *b,vector *c  ){
register int i,j;
register complex x,y;
    for(i=0;i<DIMF;i++){
	x.real=x.imag=0.0;
	for(j=0;j<DIMF;j++){
	    CMUL( a->e[i][j] , b->c[j] , y )
	    CSUM( x , y );
	}
	c->c[i] = x;
    }
}


