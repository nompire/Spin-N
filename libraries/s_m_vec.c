/******************  s_m_vec.c  (in su2.a) ******************************
*									*
*  C <- s*A,  A and C vectors 						*
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/sp.h"

/* c <- s*a, vectors */
void scalar_mult_vec( vector *a, double s, vector *c){


register int i;
    for(i=0;i<DIMF;i++){
	c->c[i].real = s*a->c[i].real;
	c->c[i].imag = s*a->c[i].imag;
    }


}
