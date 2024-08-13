// -----------------------------------------------------------------
// Update the momentum matrices
#include "sp_includes.h"



// -----------------------------------------------------------------



// -----------------------------------------------------------------

// -----------------------------------------------------------------
double gauge_force(Real eps) {
  register int i,k, dir1, dir2;
  register site *st;
  
  msg_tag *tag0, *tag1, *tag2;
  int start;
  matrix ftemp,tmat1, tmat2,tmat3,tmat4,tmat5;
  complex ctmp ;
  double norm = 0.0;
  


   



  // Loop over directions, update mom[dir1]
  for (dir1 = XUP; dir1 <= TUP; dir1++) {
    start = 1; // Indicates staple sum not initialized
    FORALLSITES(i, st)
      clear_mat(&(st->staple));    // Initialize staple

    // Loop over other directions
    // Compute force from plaquettes in the dir1, dir2 plane
    for (dir2 = XUP; dir2 <= TUP; dir2++) {
      if (dir2 != dir1) {
        // Get link[dir2] from direction dir1
        tag0 = start_gather_site(F_OFFSET(link[dir2]),
                                 sizeof(matrix),
                                 dir1, EVENANDODD, gen_pt[0]);

        // Start gather for the "upper staple"
       tag2 = start_gather_site(F_OFFSET(link[dir1]),
                                sizeof(matrix),
                                dir2, EVENANDODD, gen_pt[2]);

        // Begin the computation "at the dir2DOWN point"
        // We will later gather the intermediate result "to the home point"
        wait_gather(tag0);
        FORALLSITES(i, st) {
          mult_mat_an(&(st->link[dir2]), &(st->link[dir1]), &tmat1);
          mult_mat_nn(&tmat1, (matrix *)gen_pt[0][i],
                      (matrix *)&(st->tempmat1));
        }

        // Gather lower staple "up to home site"
        tag1 = start_gather_site(F_OFFSET(tempmat1), sizeof(matrix),
                                 OPP_DIR(dir2), EVENANDODD, gen_pt[1]);

        // The "upper" staple
        // One of the links has already been gathered,
        // since it was used in computing
        // the "lower" staple of the site above (in dir2)
        wait_gather(tag2);
        if (start) {  // This is the first contribution to staple
          FORALLSITES(i, st) {
            mult_mat_nn(&(st->link[dir2]), (matrix *)gen_pt[2][i], &tmat1);
            mult_mat_na(&tmat1, (matrix *)gen_pt[0][i], &(st->staple));
            scalar_mult_matrix(&(st->staple), 1.0 ,&(st->staple));
          }
          start = 0;
        }
        else {
          FORALLSITES(i, st) {
            mult_mat_nn(&(st->link[dir2]), (matrix *)gen_pt[2][i], &tmat1);
            mult_mat_na(&tmat1, (matrix *)gen_pt[0][i], &tmat2);
            scalar_mult_add_matrix(&(st->staple),&tmat2, 1.0 ,&(st->staple));
          }
        }

        wait_gather(tag1);
        FORALLSITES(i, st) {
          
          scalar_mult_add_matrix(&(st->staple),(matrix *)gen_pt[1][i], 1.0 ,&(st->staple));
        }
        cleanup_gather(tag0);
        cleanup_gather(tag1);
        cleanup_gather(tag2);
      }
    } // End of loop over other directions

    // Now multiply the staple sum by the link, then update momentum
    FORALLSITES(i, st) {
      
      clear_mat(&ftemp);
      for(k=0;k<NUMGEN;k++){
      
	      mult_mat_an(&(st->staple),&(Lambda[k]),&tmat3);
	      mult_mat_nn(&tmat3,&(st->link[dir1]),&tmat5);
              ctmp = trace(&tmat5);
              scalar_mult_add_matrix(&ftemp,&(Lambda[k]),2.0*ctmp.real,&ftemp);	      

      }
           
      /*mult_mat_na(&(st->staple), &(st->link[dir1]) ,&tmat3);
      adjoint(&(st->staple) ,&tmat5);
      mult_mat_nn(&(st->link[dir1]), &tmat5, &tmat1);
      
      scalar_mult_add_matrix(&tmat1 ,&tmat3 , -1.0 , &(st->staple));
      unit_mat(&tmat4);
      scalar_mult_matrix(&tmat4 ,-0.25 ,&tmat4);
      ctmp =trace(&(st->staple)) ;
      c_scalar_mult_add_mat(&(st->staple),&tmat4,&ctmp,&(st->staple));

      scalar_mult_matrix(&(st->staple) , -1.0*(BETA/8.0) ,&tmat1);
      */
      
      //Finally updating the momenta
      scalar_mult_matrix(&ftemp,(BETA/8.0),&tmat1); 
      
      scalar_mult_add_matrix(&(st->mom[dir1]), &tmat1, eps, &(st->mom[dir1]));
      norm +=  (double)realtrace(&tmat1, &tmat1);
    }
  } // End of loop over dir1

  g_doublesum(&norm);
  
  return (eps * sqrt(norm) / volume);
}
// -----------------------------------------------------------------



double scalar_force(Real eps) {
      register int i;
      register site *s;
      double dt = -1.0 * (double) eps;   // To subtract
      double returnit = 0.0;
      //matrix tmp,tmat;
      //complex tr=cmplx(0.0,0.0);
      //int c;
      FORALLSITES(i, s){
           /*clear_mat(&tmp);
           for(c=0 ; c < NUMYUK; c++){

	        mult_mat_nn(&(s->sigma),&Lambda2[c],&tmat);
		tr = trace(&tmat);
		scalar_mult_add_matrix(&tmp,&Lambda2[c],tr.real,&tmp);

	   }*/

	   mat_copy(&(s->sigma),&(s->f_sigma));

       }
       // Just subtract sigma from the momenta and compute average force
       FORALLSITES(i, s) {
	     scalar_mult_add_matrix(&(s->p_sigma), &(s->f_sigma), dt, &(s->p_sigma));
	     returnit += realtrace(&(s->f_sigma),&(s->f_sigma)); 
       }
       g_doublesum(&returnit);
       return (eps * sqrt(returnit) / volume);
}
			  

//Assume CG has been run and the solution is in sol[n]
double fermion_force(Real eps, vector *src, vector **sol){

Real ferm_epsilon =  eps;
register site *st;
register int i,dir;
int n;
double norm   =0.0;
double norm_s =0.0;
//Zero the fermion force collectors

for (dir = XUP; dir <= TUP; dir++) {
     FORALLSITES(i, st) {
      
      clear_mat(&(sigma[dir][i]));
      
        }

}
FORALLSITES(i,st){

clear_mat(&(sigma_phi[i]));

}


//Hit the solution vectors with the fermion operator


for (n = 0; n < Norder; n++) {
//Hit the solution vector with M
fermion_op(sol[n], tempvec , PLUS); 

//Compute the fermion force

compute_fermion_force(tempvec,sol[n]);

for(dir = XUP ; dir <=TUP ; dir++){
FORALLSITES(i,st) { 

scalar_mult_add_matrix(&(sigma[dir][i]),&(st->f_U[dir]),-1.0 * amp4[n], &(sigma[dir][i]));
 //clear_mat(&(sigma[dir][i])); //Uncommenting this line sets the fermion force to zero

         }

    }
FORALLSITES(i,st) {
  scalar_mult_add_matrix(&(sigma_phi[i]),&(st->f_sigma),-1.0 * amp4[n], &(sigma_phi[i])); 
  //clear_mat(&(sigma_phi[i]));
  }
		  

}

for(dir = XUP ; dir <=TUP ; dir++){

FORALLSITES(i,st){     
//Update the momenta with the fermion force 
scalar_mult_add_matrix(&(st->mom[dir]), &sigma[dir][i],ferm_epsilon, &(st->mom[dir]));
      
norm += (double)realtrace(&sigma[dir][i], &sigma[dir][i]);
     
    }

}
 FORALLSITES(i,st){

   scalar_mult_add_matrix(&(st->p_sigma), &sigma_phi[i],ferm_epsilon, &(st->p_sigma));
   norm_s +=(double)realtrace(&sigma_phi[i],&sigma_phi[i]);
 }
 



  g_doublesum(&norm_s);
  g_doublesum(&norm);
  
  return (ferm_epsilon * sqrt(norm) / volume );
  
}
// -----------------------------------------------------------------

void compute_fermion_force(vector *psol,vector *sol){


 register int i, dir;
 register site *s;
 int dumcoords[NDIMS];
 double tr,real;
 complex dum,Tr;
 int c,d;
 vector tvec_psol, tvec_sol;
 matrix tmat, tmat1,tmat2,temp,temp1,tmp,tmp2;
 msg_tag *tag0[NDIMS] ,*tag1[NDIMS];
 //Real link_mass =  site_mass;



for (dir = XUP; dir <= TUP; dir++) {
FORALLSITES(i,s) { 
         
       clear_mat(&(s->f_U[dir])); 
      
      } 

}

FORALLSITES(i,s){

clear_mat(&(s->f_sigma));
}

FORALLSITES(i,s){

vec_copy(&(psol[i]),&(src[i]));
vec_copy(&(sol[i]),&(dest[i]));

}






for (dir = XUP; dir <= TUP; dir++) {

  //Start gathers for psol[n] = M sol[n]  , and sol[n]
   

   
   tag0[dir]= start_gather_field(dest, sizeof(vector), dir,
                                  EVENANDODD, gen_pt[dir]); 

   tag1[dir] = start_gather_field(src, sizeof(vector), dir,
                                  EVENANDODD, gen_pt[4+dir]); 

}


 

for (dir = XUP; dir <= TUP; dir++) {
 wait_gather(tag0[dir]);
 wait_gather(tag1[dir]);

 FORALLSITES(i,s){

  clear_mat(&tmp);
  dumcoords[0]=s->x;
  dumcoords[1]=s->y;
  dumcoords[2]=s->z;
  dumcoords[3]=s->t;

  vec_copy((vector *)gen_pt[dir][i], &tvec_sol);
      
  if (dir == TUP && PBC < 0 && s->t == nt - 1){  scalar_mult_vec(&tvec_sol, -1.0, &tvec_sol);}
  projector(&tvec_sol, &(psol[i]),&temp,PLUS);
 
 
   
   
 //Generator version for forces
 //-----------------------------------------------------------------
 for(c = 0; c < NUMGEN; c++){
  
 if(s->parity == EVEN){
  mult_mat_nn(&(s->link[dir]),&temp,&tmat);
  mult_mat_nn(&(Lambda[c]),&tmat,&tmat1);
  dum = trace(&tmat1);
  }
  else{ 
  mult_mat_cn(&(s->link[dir]),&temp,&tmat);
  mult_mat_cn(&(Lambda[c]), &tmat,&tmat1);
  dum =trace(&tmat1);}
  
  if(dumcoords[dir]%2==0){
  tr = -0.5 * dum.real * (s->phase[dir] + link_mass * s->phase[dir]); 
  }
  else{
  tr = -0.5 * dum.real * (s->phase[dir] - link_mass * s->phase[dir]);
  }
  
  scalar_mult_add_matrix(&tmp,&(Lambda[c]),tr,&tmp);
 }
 //-----------------------------------------------------------------
  vec_copy((vector *)gen_pt[4+dir][i], &tvec_psol);

  if (dir == TUP && PBC < 0 && s->t == nt - 1){scalar_mult_vec(&tvec_psol, -1.0, &tvec_psol); }
  projector(&(sol[i]),&tvec_psol,&temp,MINUS);

  

  
  for(d = 0; d < NUMGEN; d++){
  if(s->parity == EVEN){
  mult_mat_nn(&(s->link[dir]),&temp,&tmat);
  mult_mat_nn(&(Lambda[d]),&tmat,&tmat1);
  dum = trace(&tmat1);}

  else{ 
  mult_mat_cn(&(s->link[dir]),&temp,&tmat);
  mult_mat_cn(&(Lambda[d]),&tmat,&tmat1);
  dum = trace(&tmat1);}
  

  if(dumcoords[dir]%2==0){
  tr = 0.5 * dum.real * (s->phase[dir] + link_mass * s->phase[dir]);
  }
 
  else{
  tr = 0.5 * dum.real * (s->phase[dir] - link_mass * s->phase[dir]);
  }
  



   scalar_mult_add_matrix(&tmp,&(Lambda[d]),tr,&tmp);
   }
   
  scalar_mult_matrix(&tmp,2.0,&(s->f_U[dir]));
  
 }
  cleanup_gather(tag0[dir]);
  cleanup_gather(tag1[dir]);
} 

  //-------------------------------f_sigma = dS/dsigma-----------------------------//
     FORALLSITES(i,s){
       	dumcoords[0] = s->x ;
        dumcoords[1] = s->y ;
        dumcoords[2] = s->z ;
        dumcoords[3] = s->t ;
        clear_mat(&tmp);
               // Fermion force contribution for scalar 

 	      projector(&(sol[i]),&(psol[i]),&tmat,PLUS);
	      for(c=0;c<NUMYUK;c++){
	           mult_mat_nn(&(Lambda2[c]),&tmat,&tmat2);
             Tr = trace(&tmat2);
             real = 0.5 * G * Tr.real;
             scalar_mult_add_matrix(&tmp,&(Lambda2[c]),real*parity(dumcoords),&tmp);
	      }
  
        scalar_mult_matrix(&tmp,-2.0,&tmat);
        mat_copy(&tmat,&(s->f_sigma));
  
  
     }
  


}
  

