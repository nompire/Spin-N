#include "../include/config.h"
#include "../include/complex.h"
#include "../include/sp.h"
#include "../include/macros.h"
#include "../RHMC/sp_includes.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
void exp_matrix(matrix *u, matrix *v){

        matrix *del = malloc(sizeof(*del));
        matrix *prod = malloc(sizeof(*prod));
	    matrix *tmat = malloc(sizeof(*tmat));
	    double fac=1.0;
	    
        int i=1;
        int j,k,l;
        for( j=0; j < DIMF ; j++){
           for( k=0; k < DIMF ; k++){
             prod->e[j][k] = cmplx(0.0,0.0);  
             v->e[j][k] =cmplx(0.0,0.0); 
           }
        }
 
        for(l=0 ; l < DIMF ; l++){
             prod->e[l][l] = cmplx(1.0,0.0);
             v->e[l][l] = cmplx(1.0,0.0);
        }
        
        static int sum=0,counter=0;
	
        do{
        fac=fac*(double)i;
        mult_mat_nn(prod,u,tmat);
        //prod=prod*u;
        scalar_mult_matrix(tmat, (1.0/fac),del);
        //del=prod*(1.0/fac);
        add_mat(v,del);
        mat_copy(tmat,prod);
        //v=v+del;
        i++;
        }
        	
        while(sqrt(realtrace(del,del))>GAUGETOL);
        //node0_printf("No of terms computed for matrix exponential : %d\n", i);
        sum+=i;
        counter++;
        if(counter==100000){ 
           node0_printf("mean no. of terms in exp()=%.6g\n" ,(double)(sum/counter) ); counter=0;sum=0;}
        free(del);
        free(prod);
        free(tmat);
    

}
