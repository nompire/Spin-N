/********************  clear_mat.c  ********************
*/

#include "../include/config.h"
#include "../include/complex.h"
#include "../include/sp.h"

void clear_mat(matrix *dest ){
register int i,j;
    for(i=0;i<DIMF;i++)for(j=0;j<DIMF;j++){
	dest->e[i][j].real = dest->e[i][j].imag = 0.0;
    }
}
