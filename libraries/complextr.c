/******************  complextr.c****************************
*									*
* complex complextrace_(matrix *a,*b)				*
* return Tr( A_adjoint*B )   						*
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/sp.h"
#include <stdio.h>

complex complextrace( matrix *a,matrix *b) {
register int i,j;
register Real sumr, sumi;
complex tr;
    sumr=0.0,sumi=0.0;
    for(i=0;i<DIMF;i++)for(j=0;j<DIMF;j++){
     sumr+= a->e[i][j].real*b->e[i][j].real + a->e[i][j].imag*b->e[i][j].imag;
     sumi+= a->e[i][j].real*b->e[i][j].imag - a->e[i][j].imag*b->e[i][j].real;
    }
    
   tr.real= sumr; tr.imag=sumi; 
   return(tr); 
}
