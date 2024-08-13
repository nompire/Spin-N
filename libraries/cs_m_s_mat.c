/****************  cs_m_s_mat.c             *****************************
*									*
*                                                              		*
*					                                *
* C <- A - s*B,   A,B and C matrices 					*
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/sp.h"

/* c <- a - s*b, matrices */
void c_scalar_mult_sub_mat( matrix *a, matrix *b, complex *s,matrix *c){


register int i,j;
register double sr,si;
sr = (*s).real; si = (*s).imag;
 
    for(i=0;i<DIMF;i++){
    for(j=0;j<DIMF;j++){
	
         c->e[i][j].real = a->e[i][j].real -  (sr*b->e[i][j].real -  si*b->e[i][j].imag);
	 c->e[i][j].imag = a->e[i][j].imag -  (sr*b->e[i][j].imag +  si*b->e[i][j].real); 

    }

  }

}
