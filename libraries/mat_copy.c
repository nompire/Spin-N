/*****************  su2mat_copy.c  (in su2.a) ***************************
*									*
* Copy an su2 matrix:  B <- A   						*
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/sp.h"

/* Copy a  matrix:  b <- a   */
void mat_copy(matrix *a,matrix *b ){
register int i,j;
    for(i=0;i<DIMF;i++){
     for(j=0;j<DIMF;j++){
	b->e[i][j].real = a->e[i][j].real;
	b->e[i][j].imag = a->e[i][j].imag;
      }
   }
}
