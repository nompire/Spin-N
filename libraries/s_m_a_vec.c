/****************  s_m_a_vec.c             ******************************
*									*
* C <- A + s*B,   A,B and C vectors 					*
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/sp.h"

/* c <- a + s*b, vectors */

void scalar_mult_add_vec(vector *a, vector *b, double s,vector *c){

register int i;
  for(i=0;i<DIMF;i++){
    c->c[i].real = a->c[i].real + s*b->c[i].real;
    c->c[i].imag = a->c[i].imag + s*b->c[i].imag;
  }
  
}
