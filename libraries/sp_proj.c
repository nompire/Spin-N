/*****************  su2_proj.c  (in su2.a) ******************************
*									*
* void su2_projector( su2_vector *a, su2_vector *b, su2_matrix *c )	*
* C  <- outer product of A and B					*
*  C_ij = A_i * B_adjoint_j						*
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/sp.h"


void projector(vector *u,vector *v, matrix *w, int flag){
register int a,b;
    for(a=0;a<DIMF;a++)for(b=0;b<DIMF;b++){
	if(flag == 1){    
	  w->e[b][a].real = u->c[b].real * v->c[a].real + u->c[b].imag * v->c[a].imag;
	  w->e[b][a].imag = u->c[b].imag * v->c[a].real - u->c[b].real * v->c[a].imag; 
	}
	else{
	  w->e[a][b].real = u->c[b].real * v->c[a].real + u->c[b].imag * v->c[a].imag;
          w->e[a][b].imag = u->c[b].imag * v->c[a].real - u->c[b].real * v->c[a].imag;
  
	}
    }
}


