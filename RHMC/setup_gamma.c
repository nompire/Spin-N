#include "../include/complex.h"
#include "lattice.h"
#include "math.h"
#include "sp_includes.h"


#define IMAG_TOL 1.0e-8

#define DEBUC_CHECK


void setup_gamma(){

//setting up Euclidean Dirac gamma matrices//
matrix tmp,tmat,tmat1;
clear_mat(&tmp);
clear_mat(&tmat);
clear_mat(&tmat1);
register int i,j,l;
complex tr,tc;
complex half = cmplx(0.0,0.5);
//initialize to zero

for(i=0;i<DIMF+1;i++){

     clear_mat(&(gamma_dirac[i]));
}
for(i=0;i<NUMGEN;i++){
     clear_mat(&(Lambda[i]));
}    	
for(j=0;j<NUMYUK+1;j++){
     clear_mat(&(Lambda2[j]));
}	
node0_printf("Computing generators for SP(N)\n");
// initialize non-zero elements of gamma matrices


  gamma_dirac[0].e[0][2] = cmplx(1.0,0.0);
  gamma_dirac[0].e[2][0] = cmplx(1.0,0.0);   
  gamma_dirac[0].e[1][3] = cmplx(1.0,0.0);
  gamma_dirac[0].e[3][1] = cmplx(1.0,0.0);

   
  gamma_dirac[1].e[0][3] = cmplx(0.0,1.0);
  gamma_dirac[1].e[1][2] = cmplx(0.0,1.0); 
  gamma_dirac[1].e[2][1] = cmplx(0.0,-1.0);
  gamma_dirac[1].e[3][0] = cmplx(0.0,-1.0);


  gamma_dirac[2].e[0][3] = cmplx(-1.0,0.0);
  gamma_dirac[2].e[1][2] = cmplx(1.0,0.0); 
  gamma_dirac[2].e[2][1] = cmplx(1.0,0.0);
  gamma_dirac[2].e[3][0] = cmplx(-1.0,0.0);

  gamma_dirac[3].e[0][2] = cmplx(0.0,1.0);
  gamma_dirac[3].e[2][0] = cmplx(0.0,-1.0); 
  gamma_dirac[3].e[1][3] = cmplx(0.0,-1.0);
  gamma_dirac[3].e[3][1] = cmplx(0.0,1.0); 
 
 
//  gamma_5 
  mult_mat_nn(&(gamma_dirac[0]),&(gamma_dirac[1]),&tmp);
  mult_mat_nn(&tmp,&(gamma_dirac[2]),&tmat);
  mult_mat_nn(&tmat,&(gamma_dirac[3]),&(gamma_dirac[4]));
//print gamma matrices
for(i=0;i<DIMF+1;i++){
        node0_printf("Gamma[%d]\n",i);
        dumpmat(&(gamma_dirac[i]));
}

// check trace of gamma matrices
for(i=0;i<DIMF+1;i++){
    mult_mat_nn(&(gamma_dirac[i]),&(gamma_dirac[i]),&tmp);
    tr = trace(&tmp);
    node0_printf(" Tr(gamma_mu * gamma_mu): %.8g %.8g \n",tr.real,tr.imag);	    
    } 

// generators for Spin(4)  

   l=0;
   
   int ii,jj;
    for(ii=0;ii<DIMF;ii++) {
      for(jj=ii+1;jj<DIMF;jj++) {
       
          mult_mat_nn(&(gamma_dirac[ii]),&(gamma_dirac[jj]),&tmp);
          scalar_mult_matrix(&tmp,0.5,&(Lambda[l]));
          l++;  	  
	
      }
       
    }

  



// generators for Spin(5)
if(NUMGEN == 10){
  
   j=6;
   for (i=0;i<DIMF;i++){
        
  	 c_scalar_mult_mat(&(gamma_dirac[i]), &half, &(Lambda[i+j]));
   }
}



// generators for Spin(6)
if ((NUMGEN) == 15){
    
  j = 10;
  for (i = 0; i < DIMF; i++){

     mult_mat_nn(&(gamma_dirac[4]), &(gamma_dirac[i]), &tmp);
     scalar_mult_matrix(&tmp,0.5,&(Lambda[i+j]));
}
  c_scalar_mult_mat(&(gamma_dirac[4]),&half, &(Lambda[NUMGEN-1]));
}

node0_printf("Generators for SP(N)\n");
for(i=0;i<NUMGEN;i++){
        node0_printf("Lambda[%d]\n",i);
	dumpmat(&(Lambda[i]));
}

// generators for Yukawa terms


mult_mat_nn(&(gamma_dirac[1]),&(gamma_dirac[3]),&tmp);
scalar_mult_matrix(&tmp,0.5,&(Lambda2[0]));

c_scalar_mult_mat(&(gamma_dirac[1]),&half,&(Lambda2[1]));

c_scalar_mult_mat(&(gamma_dirac[3]),&half,&(Lambda2[2]));


mult_mat_nn(&(gamma_dirac[4]),&(gamma_dirac[0]),&tmp);
c_scalar_mult_mat(&tmp,&half,&(Lambda2[3]));

mult_mat_nn(&(gamma_dirac[4]),&(gamma_dirac[2]),&tmp);
c_scalar_mult_mat(&tmp,&half,&(Lambda2[4]));

mult_mat_nn(&(gamma_dirac[0]),&(gamma_dirac[2]),&tmp);
c_scalar_mult_mat(&tmp,&half,&(Lambda2[5]));

if ((NUMYUK) == 1){
mult_mat_nn(&(gamma_dirac[0]),&(gamma_dirac[2]),&tmp);
c_scalar_mult_mat(&tmp,&half,&(Lambda2[5]));}



//Just using the sp4 invariant term
mat_copy(&(Lambda2[5]),&(Lambda2[0]));

node0_printf("Generators for SP(N) invariant Yukawas\n");
for(i=0;i<NUMYUK;i++){
	node0_printf("Lambda2[%d]\n",5);
	dumpmat(&(Lambda2[5]));
}

//pseudo-real property of SP(N) generators and Yukawa generators


//for(i=0;i<NUMYUK;i++){
/*
	conjug(&(Lambda2[5]),&tmp);
        mult_mat_nn(&(Lambda2[5]),&(Lambda2[5]),&tmat);
      mult_mat_nn(&(Lambda2[5]),&tmat,&tmat1);
	scalar_mult_matrix(&tmat1,4.0,&tmat);
       	sub_matrix(&tmat,&tmp,&tmat1);
	dumpmat(&tmat1);
*/
//}
// Test orthogonality of products of Lambdas

node0_printf("Trace of Lambdas\n");
for (i = 0; i < NUMGEN; i++) {       
	 for (j = 0; j < NUMGEN; j++) {
                    mult_mat_nn(&(Lambda[i]), &(Lambda[j]), &tmat);
                    tc = trace(&tmat);
                        if (tc.real * tc.real > IMAG_TOL)
                                 node0_printf("Tr[T_%d T_%d] = (%.4g, %.4g)\n",
                                                      i, j, tc.real, tc.imag);                                                          
	 }
                                                            
 }
node0_printf("Trace of Yukawas\n");
for (i = 0; i < NUMYUK; i++) {       
	 for (j = 0; j < NUMYUK; j++) {
                    mult_mat_nn(&(Lambda2[i]), &(Lambda2[j]), &tmat);
                    tc = trace(&tmat);
                        if (tc.real * tc.real > IMAG_TOL)
                                 node0_printf("Tr[T_%d T_%d] = (%.4g, %.4g)\n",
                                                      i, j, tc.real, tc.imag);                                                          
	 }
                                                            
 }

node0_printf("Gamma matrices setup completed\n");
}


