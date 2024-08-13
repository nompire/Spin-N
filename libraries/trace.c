/*******************  trace.c                 ***************************
*									*
* 				                                   	*
* return complex trace of a  matrix 				        *
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/sp.h"

/* Complex trace of a matrix */
complex trace(matrix *a) {
    complex tr = cmplx(0.0,0.0);

    register int i;
    
    for(i=0;i<DIMF;i++){
      
	    CSUM(tr,a->e[i][i]);
    } 
    
    return(tr);
   
    
}
