/******************  det.c                 ******************************
*									*
*                                         				*
* Complex determinant of matrix 					*
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/sp.h"
#include <stdio.h>

/* FIX THIS - more efficient to take cross product of first two
   rows, dot with third. */
void det(matrix *a) {
register complex cc,dd,s;
 double d;
    CMUL(a->e[0][0],a->e[1][1],cc);
    CMUL(a->e[0][1],a->e[1][0],dd);
    CSUB(cc,dd,s);
    d=s.real;
    //d = (a->e[0][0].real * a->e[1][1].real) - ( a->e[0][0].imag * a->e[1][1].imag) - ( a->e[0][1].real * a->e[1][0].real) + (a->e[0][1].imag * a->e[1][0].imag); 
     
    printf("Det = %.8g\n",d);
    

    
}
