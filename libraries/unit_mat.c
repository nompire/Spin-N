/********************  unit_mat.c          ********************
*

* dest  <-  unit_matrix
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/sp.h"

void unit_mat(matrix *dest ){
register int i,j;
    for(i=0;i<DIMF;i++){
    for(j=0;j<DIMF;j++){
    if( i != j){dest->e[i][j].real = 0.0;  dest->e[i][j].imag = 0.0;}
    else{dest->e[i][j].real = 1.0; dest->e[i][j].imag = 0.0;}
    }
 }

}
