// -----------------------------------------------------------------
// Print the given vector
#include "../include/config.h"
#include <stdio.h>
#include "../include/complex.h"
#include "../include/sp.h"

void dumpvec(vector *v) {

  register int i;
  
  for (i = 0; i < DIMF; i++){
    printf("  %.4g", v->c[i].real);
    printf("  %.4g", v->c[i].imag);
    printf("\n");

}}
// -----------------------------------------------------------------
