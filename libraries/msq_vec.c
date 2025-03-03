// -----------------------------------------------------------------
// Squared magnitude of vector
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/sp.h"

double magsq_vec(vector *a) {
    register double sum ;

    sum=0.0;
    register int i;
    for(i=0;i<DIMF;i++){sum += a->c[i].real*a->c[i].real  + a->c[i].imag*a->c[i].imag;}
    

   return(sum);
}
// -----------------------------------------------------------------
