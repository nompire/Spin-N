/****************  s_m_a_mat.c             ******************************
*									*
*            								*
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/sp.h"

/* c <- c+ s*a, matrices */
void scalar_mult_sum_matrix(matrix *a,double s,
	matrix *c){


register int i,j;
    for(i=0;i<DIMF;i++){
    for(j=0;j<DIMF;j++){
	c->e[i][j].real += s*a->e[i][j].real ;
	c->e[i][j].imag += s*a->e[i][j].imag ; 
       }
    }

}
