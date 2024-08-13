/****************  s_m_a_mat.c             ******************************
*									*
* 					*
* C <- A + s*B								*
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/sp.h"

/* c <- a + s*b, matrices */
void scalar_mult_add_matrix(matrix *a,matrix *b,double s,
	matrix *c){


register int i,j;
    for(i=0;i<DIMF;i++){
    for(j=0;j<DIMF;j++){
	c->e[i][j].real = a->e[i][j].real + s*b->e[i][j].real;
	c->e[i][j].imag = a->e[i][j].imag + s*b->e[i][j].imag;
       }
    }

}
