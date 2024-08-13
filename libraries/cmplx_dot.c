/******************   ******************************
*									*
* complex cmplx_dot(vector *a, vector *b )			*
* return dot product of two vectors: a^dagger b			*
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/sp.h"

complex cmplx_dot(vector *a,vector *b ){


complex temp1,temp2,temp3,temp4,t,u,s;
    CMULJ_(a->c[0],b->c[0],temp1); 
    CMULJ_(a->c[1],b->c[1],temp2);
    CMULJ_(a->c[2],b->c[2],temp3);
    CMULJ_(a->c[3],b->c[3],temp4);
    CADD(temp1,temp2,t)
    CADD(temp3,temp4,u);
    CADD(t,u,s);
    return(s);


}
