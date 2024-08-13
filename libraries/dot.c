// -----------------------------------------------------------------
// Return the dot product of two vectors: adag b
#include "../include/config.h"
#include "../include/sp.h"
#include "../include/complex.h"
#include <stdio.h>
#include <stdlib.h>

double dot(vector *a, vector *b) {
register complex x,y;
x=cmplx(0.0,0.0);
register int i;
  
  for(i = 0; i < DIMF; i++){
           
       CMULJ_(a->c[i],b->c[i],y);
       CSUM(x,y);
          
  }

  return(x.real);  
}
// -----------------------------------------------------------------
