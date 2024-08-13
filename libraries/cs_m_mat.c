/****************  cs_m_mat.c             *******************************
*									*
*                                                               	*
* C <- s*B,   B and C matrices 						*
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/sp.h"

/* c <- s*b, matrices */
void c_scalar_mult_mat(matrix *b, complex *s,matrix *c ){

register double sr,si;
sr = (*s).real; si = (*s).imag;
register int i,j;
    for(i=0;i<DIMF;i++)
    {for(j=0;j<DIMF;j++){
	//CMUL(&b->e[i][j], s,c->e[i][j]);
	 c->e[i][j].real = sr*b->e[i][j].real -  si*b->e[i][j].imag;
	 c->e[i][j].imag = sr*b->e[i][j].imag +  si*b->e[i][j].real; 
    }


}

}
