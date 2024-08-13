/*************************** ranmom.c *******************************/
/* Produce Gaussian random momenta for the gauge fields. */

#include "generic_includes.h"
//#include <defines.h>                 /* For SITERAND */


void ranmom() {
  register int i, j, mu;
  register site *s;
  double grn;

  FORALLSITES(i, s) {
    for(mu = XUP ; mu <=TUP ; mu++) {
      clear_mat(&(s->mom[mu]));
      for (j = 0; j < NUMGEN; j++) {
#ifdef SITERAND
        grn = gaussian_rand_no(&(s->site_prn));
       
#else
        
        grn = gaussian_rand_no(&node_prn);
        
#endif
        scalar_mult_add_matrix(&(s->mom[mu]) ,&(Lambda[j]), grn, &(s->mom[mu]));
      }
    }
  }




  FORALLSITES(i,s){
	          clear_mat(&(s->p_sigma));
		  for(j=0 ; j < NUMYUK ; j++) {
		  #ifdef SITERAND
		      grn =  gaussian_rand_no(&(s->site_prn));
	          #else
		      grn =  gaussian_rand_no(&node_prn);
		  #endif
	              scalar_mult_add_matrix(&(s->p_sigma) ,&(Lambda2[j]), grn, &(s->p_sigma));
                  }
  }

}


