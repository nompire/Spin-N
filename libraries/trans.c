/******************  trans.c   **************************
*									*
* B  <- A_trans,  transpose of a matrix 				*
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/sp.h"

/* transpose of an  matrix */
void trans(matrix *a, matrix *b ){
register int i,j;
    for(i=0;i<DIMF;i++){
    for(j=0;j<DIMF;j++){
        b->e[i][j].real = a->e[j][i].real;
        b->e[i][j].imag = a->e[j][i].imag;
    }
  }
}
