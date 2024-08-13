/******************  cs_m_a_mat.c             ***************************
*									*
*                                         				*
*  multiply an  matrix by a complex scalar and add it to another	*
*  matrix:   m3 <- m1 + number*m2 					*
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/sp.h"

void c_scalar_mult_add_mat( matrix *m1, matrix *m2,
	complex *phase, matrix *m3){


register int i,j;
register double sr,si;
sr = (*phase).real; si = (*phase).imag;

    for(i=0;i<DIMF;i++){
    for(j=0;j<DIMF;j++){
	 m3->e[i][j].real = m1->e[i][j].real  +   (sr*m2->e[i][j].real - si*m2->e[i][j].imag);
	 m3->e[i][j].imag = m1->e[i][j].imag +   (sr*m2->e[i][j].imag +  si*m2->e[i][j].real);
        }


   }
}
