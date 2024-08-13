/****************  m_matvec_s.c *****************************
*                                                           *
* C  <-  C + A*B                                            *
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/sp.h"

//#ifndef FAST
/* matrix times vector multiply and add to another vector */
/* c  <-  A*b+c */
void mult_mat_vec_sum(matrix *a,vector *b, vector *c ){
register int i,j;
register complex x,y;
    for(i=0;i<DIMF;i++){
  x.real=x.imag=0.0;
  for(j=0;j<DIMF;j++){
      CMUL( a->e[i][j] , b->c[j] , y )
      CSUM( x , y );
  }
  c->c[i].real += x.real;
  c->c[i].imag += x.imag;
    }
}


