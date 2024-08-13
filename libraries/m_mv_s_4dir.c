/****************  m_mv_s_4dir.c            *****************************
*									*
* void mult_mat_vec_sum_4dir( matrix *a, vector *b[0123],*c )*
* Multiply the elements of an array of four su2_matrices by the		*
* four su2_vectors, and add the results to				*
* produce a single su2_vector.						*
* C  <-  A[0]*B[0]+A[1]*B[1]+A[2]*B[2]+A[3]*B[3]			*
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/sp.h"

void mult_mat_vec_sum_4dir(matrix *a, vector *b0,
       vector *b1, vector *b2, vector *b3, vector *c  ){
    mult_mat_vec( a+0,b0,c );
    mult_mat_vec_sum( a+1,b1,c );
    mult_mat_vec_sum( a+2,b2,c );
    mult_mat_vec_sum( a+3,b3,c );
}


