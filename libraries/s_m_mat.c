/******************  s_m_mat.c             ******************************
*									*
* B <- s*A								*
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/sp.h"

/* b <- s*a, matrices */
void scalar_mult_matrix( matrix *a, double s,matrix *b ){

register int i,j;
    for(i=0;i<DIMF;i++)for(j=0;j<DIMF;j++){
	b->e[i][j].real = s*a->e[i][j].real;
	b->e[i][j].imag = s*a->e[i][j].imag;
    }
}


