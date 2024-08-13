
// -----------------------------------------------------------------
// Measurements for SP(N) gauge model with red. stagg'd fermions correlator, bilinear condensate and one-link condensate
#include "sp_includes.h"
#include "../include/sp.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Set dest to unit source at given point in given SP(N) index
// Then set src = Ddag dest
// so that (Ddag.D)^(-1).src will give D^(-1).pnt_src
// Return the number of iterations from the inversion
void pnt_src(int *pnt, int index) {
  register int i;
  register site *s;

  // Set dest to unit source at given point in given SU(2) index
  FORALLSITES(i, s){
    clearvec(&(dest[i]));
}
  
  if (node_number(pnt[0], pnt[1], pnt[2], pnt[3]) == mynode()) {
    i = node_index(pnt[0], pnt[1], pnt[2], pnt[3]);
    dest[i].c[index] = cmplx(1.0,0.0) ; 
  }

  adj_fermion_op(dest, src);
  
  
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// For some reason this doesn't work in the libraries
Real Z2_rand_no(double_prn *prn_pt) {
  if (myrand(prn_pt) > 0.5)
    return 1.0;
  else
    return -1.0;
}

// Set dest to random Z2 source at all sites
// Then set src = Ddag dest
// so that (Ddag.D)^(-1).src will give D^(-1).vol_src
void vol_src() {
  register int i,j;
  double inv_sqrt = 1.0/sqrt(2.0);
  register site *s;

  FORALLSITES(i,s){
#ifdef SITERAND
	  for(j = 0; j < DIMF; j++){
           
		  dest[i].c[j].real =  Z2_rand_no(&(s->site_prn));
		  dest[i].c[j].imag = 0.0;//inv_sqrt * Z2_rand_no(&(s->site_prn));
	  
	  }

#else
	  for(j = 0; j < DIMF; j++){

                  dest[i].c[j].real =  Z2_rand_no(&node_prn);
                  dest[i].c[j].imag = 0.0;//inv_sqrt * Z2_rand_no(&node_prn);

          }

#endif
  }
  adj_fermion_op(dest, src);
}
// -----------------------------------------------------------------


// -----------------------------------------------------------------
// Measure bilinear condensate
// Use Z2 stochastic sources
// Return total number of iterations
int condensates() {
  register int i,dir;
  register site *s;
  int a, b, c, d,  j, iters, tot_iters = 0, sav = Norder;
  
  Real size_r;
  double norm = 1.0 / (double)Nstoch;
  int dumcoords[NDIMS];
  double xi ;
  complex id = cmplx(1.0,0.0);
  complex four = cmplx(0.0,0.0);
  complex bilin = cmplx(0.0,0.0);
  complex one_link = cmplx(0.0,0.0);
  double dtime;
  vector **psim;
  vector tvec_dir , tvec_opp;
  msg_tag *mtag[NDIMS] , *tag[2 * NDIMS];
  matrix *tm, tmat;

  // Initialize stochastic propagators
  FORALLSITES(i, s) {
    for (a = 0; a < DIMF; a++) {
      for (b = 0; b < DIMF; b++) {
        prop[i].e[a][b] = cmplx(0.0,0.0);  //prop is a matrix 
        prop2[i].e[a][b] = cmplx(0.0,0.0);
        prop3[i].e[a][b] = cmplx(0.0,0.0);
      }
    }
  }

   for (dir =XUP ;dir <=TUP ; dir++){
   FORALLSITES(i,s){
    clear_mat(&(s->temp_link[dir]));
      
    if(lattice[i].parity == EVEN){mat_copy(&(s->link[dir]),&(s->temp_link[dir])) ;}
     
    else{conjug(&(s->link[dir]), &(s->temp_link[dir]));} 
    
   
   }

}
  // Hack a basic CG out of the multi-mass CG
  Norder = 1;
  psim = malloc(sizeof(**psim));
  psim[0] = malloc(sites_on_node * sizeof(vector));
  shift[0] = 0;


  // Construct stochastic propagator
  //   D_{ab}^{-1}(x) = (1/N) sum_N psi_a(x) dest_b(x)
  // where dest are random Z2 sources
  // Hit each dest with Mdag to get src_j, invert to get D_{kj}^{-1} dest_j
  for (j = 0; j < Nstoch; j++) {
    

    dtime = -dclock();
    vol_src();
    iters = congrad_multi(src, psim, niter, rsqmin, &size_r);
    dtime += dclock();
    tot_iters += iters;
    node0_printf("Inversion %d-1 of %d took %d iters and %.4g seconds\n",
                 j + 1, Nstoch, iters, dtime);


    // Copy psim into f[k][j]
    FORALLSITES(i, s) {
      
       proj(&(dest[i]), &(psim[0][i]), &(prop[i]));  
     }
    


    dtime = -dclock();
    vol_src();
    iters = congrad_multi(src, psim, niter, rsqmin, &size_r);
    dtime += dclock();
    tot_iters += iters;
    node0_printf("Inversion %d-2 of %d took %d iters and %.4g seconds\n",
                 j + 1, Nstoch, iters, dtime);

    // Copy psim into f[k][j]
    FORALLSITES(i, s) {
      
       proj(&(dest[i]), &(psim[0][i]) , &(prop2[i]));   
   }
   }
  // Normalize stochastic propagator by norm = 1 / Nstoch
  FORALLSITES(i, s) {
    
     scalar_mult_matrix(&(prop[i]), norm, &(prop[i]));
     scalar_mult_matrix(&(prop2[i]), norm, &(prop2[i]));  
  }
  // Four fermion calculation
  for(a = 0; a < DIMF; a++){
       for(b = 0; b < DIMF; b++){
             //if(b == a)
	       // continue;
		for(c = 0; c < DIMF; c++){
		   //if(c == b || c == a)
		     // continue;
		      for(d = 0; d < DIMF; d++){
		          //if(d == c || d == b || d == a)
			  //  continue;

			    FORALLSITES(i,s){

			       CMULREAL_SUM(prop[i].e[a][b],prop2[i].e[c][d], perm[a][b][c][d], four)
			    }
			    
		       }
		       
		 }
       
         }
    }


  g_complexsum(&four);
  four.real /= (double) volume;
  four.imag /= (double) volume;

  c_scalar_mult_mat(&(Lambda2[0]), &id, &tmat);
  // Add the contributions to bilin
  for (a = 0; a < DIMF; a++) {
    for (b = 0; b < DIMF; b++) {
     
      FORALLSITES(i, s) {

          dumcoords[0] = s->x ;
          dumcoords[1] = s->y ;
          dumcoords[2] = s->z ;
          dumcoords[3] = s->t ;

          bilin.real += parity(dumcoords)*(prop[i].e[a][b].real * tmat.e[a][b].real - prop[i].e[a][b].imag * tmat.e[a][b].imag);
          bilin.imag += parity(dumcoords)*(prop[i].e[a][b].imag * tmat.e[a][b].real + prop[i].e[a][b].real * tmat.e[a][b].imag);
	  bilin.real += parity(dumcoords)*(prop2[i].e[a][b].real * tmat.e[a][b].real - prop2[i].e[a][b].imag * tmat.e[a][b].imag);
          bilin.imag += parity(dumcoords)*(prop2[i].e[a][b].imag * tmat.e[a][b].real + prop2[i].e[a][b].real * tmat.e[a][b].imag);

         
        
      }
    }
  }
//------------------------------------------------------------------------------------------------------

    
  
      g_complexsum(&bilin);
      bilin.real /= 2.0 * (double)volume;
      bilin.imag /= 2.0 * (double)volume;

    
  

  // prop2/prop3 will store the propagator M^{-1}_(x,x +/- \mu) ...This requires a gather call 
//---------------------------------------------------------------------------------------------------------------------------------------------------------
   
 
  
 for (j = 0; j < Nstoch; j++) {
    

    dtime = -dclock();
    vol_src();
    iters = congrad_multi(src, psim, niter, rsqmin, &size_r);
    dtime += dclock();
    tot_iters += iters;
    node0_printf("Inversion %d-1 of %d took %d iters and %.4g seconds\n",
                 j + 1, Nstoch, iters, dtime);


    for (dir = XUP; dir <= TUP; dir++) {
      
      tag[dir] = start_gather_field(dest, sizeof(vector), dir,
                                    EVENANDODD, gen_pt[dir]);
      tag[OPP_DIR(dir)] = start_gather_field(dest, sizeof(vector), OPP_DIR(dir),EVENANDODD, gen_pt[OPP_DIR(dir)]);
    

      mtag[dir] = start_gather_site(F_OFFSET(temp_link[dir]), sizeof(matrix),OPP_DIR(dir), EVENANDODD, gen_pt[11-dir]);

  }
    
    // Copy psim into f[k][j]

    for (dir = XUP; dir <= TUP; dir++) {
   
    
    wait_gather(tag[dir]);
    wait_gather(tag[OPP_DIR(dir)]);
    wait_gather(mtag[dir]);

     

    FORALLSITES(i, s) {
      
        dumcoords[0]=s->x;
        dumcoords[1]=s->y;
        dumcoords[2]=s->z;
        dumcoords[3]=s->t;
      
       tm = (matrix *)(gen_pt[11-dir][i]);

       vec_copy((vector *)gen_pt[dir][i], &tvec_dir);
       vec_copy((vector *)gen_pt[OPP_DIR(dir)][i], &tvec_opp); 

      if (dir == TUP && PBC < 0 && s->t == nt - 1)       //Boundary conditions
 
     { scalar_mult_vec(&tvec_dir, -1.0, &tvec_dir); } 
      else if (dir == TUP && PBC < 0 && s->t == 0)
     { scalar_mult_vec(&tvec_opp, -1.0, &tvec_opp); }
       
       proj(&(psim[0][i]) , &tvec_dir, &(prop2[i]));   
       proj(&(psim[0][i]) , &tvec_opp ,&(prop3[i]));

       scalar_mult_matrix(&(prop2[i]), norm, &(prop2[i]));
       scalar_mult_matrix(&(prop3[i]), norm, &(prop3[i]));
      
      xi=1.0;
      if(dumcoords[dir] % 2 != 0 ) { xi = -1.0 ;}
                  
                  for(a=0 ; a < DIMF ; a++){
                  for(b=0 ; b < DIMF ; b++){
      
                 CMULREAL_SUM(s->temp_link[dir].e[a][b], prop2[i].e[a][b], xi * (double) s->phase[dir], one_link);
      
                 CMULREAL_SUM(tm->e[b][a], prop3[i].e[a][b], xi * (double) s->phase[dir] , one_link); 
            
               }
      
            }



        }
      
    cleanup_gather(tag[dir]);
    cleanup_gather(tag[OPP_DIR(dir)]);
    cleanup_gather(mtag[dir]);



   }


}

      g_complexsum(&(one_link));
      
      one_link.real /= 2*(double) volume  ;
      one_link.imag /= 2*(double) volume ;
      





//----------------------------------------------------------------------------------------------------------------------------------------------------------------


  // Print condensates 
  node0_printf("STOCH FOUR %.6g %.6g\n", four.real , four.imag);
  node0_printf("STOCH BILIN %.6g %.6g\n",bilin.real,bilin.imag);
  node0_printf("STOCH ONELINK  %.6g %.6g\n",one_link.real,one_link.imag);
  // Reset multi-mass CG and clean up
  Norder = sav;
  free(psim[0]);
  free(psim);
  return tot_iters;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Measure two-point  correlators
// Return total number of iterations
int correlators(int *pnt) {
  register int i;
  register site *s;
  
  int a, b, BC,iters, tot_iters = 0, sav = Norder, dir;
  int L[NDIMS] = {nx, ny, nz, nt};
  Real size_r;
  int t;
  complex corr[nt],tmp,tmp2,tmp3;
  complex id = cmplx(1.0,0.0);
  complex bilin = cmplx(0.0,0.0);//bilin has to defined complex as well
  double dtime;
  complex one_link = cmplx(0.0,0.0);//defined complex
  vector **psim;
  matrix *tm,*tm2,*tm3, tmat;
  //static int first_time=1;
  msg_tag *tag[2 * NDIMS];
  msg_tag *mtag[NDIMS];
  
  int dumcoords[NDIMS];
  int j;
  for(j=0 ; j < NDIMS ; j++){
      
    dumcoords[j] = pnt [j] ;
  }
  
  
  
  for(t=0 ; t< nt; t++){
  
  corr[t]= cmplx(0.0,0.0);
  
  }
  
   
FORALLSITES(i,s){

      clear_mat(&(prop[i]));
    }




for (dir =XUP ;dir <=TUP ; dir++){
   FORALLSITES(i,s){
    clear_mat(&(s->temp_link[dir]));
      
    if(lattice[i].parity == EVEN){mat_copy(&(s->link[dir]),&(s->temp_link[dir])) ;}
     
    else{conjug(&(s->link[dir]), &(s->temp_link[dir]));} 
    
   
   }

}
  // Make sure pnt stays within lattice volume
  for (dir = XUP; dir <= TUP; dir++) {
    i = pnt[dir] % L[dir];
    pnt[dir] = i;
  }

  // Hack a basic CG out of the multi-mass CG
  Norder = 1;
  psim = malloc(sizeof(**psim));
  psim[0] = malloc(sites_on_node * sizeof(vector));
  shift[0] = 0;

  for (a = 0; a < DIMF; a++) {
    dtime = -dclock();
    pnt_src(pnt, a);
    iters = congrad_multi(src, psim, niter, rsqmin, &size_r);
    dtime += dclock();
    tot_iters += iters;
    node0_printf("Inversion %d of %d took %d iters and %.4g seconds\n",
                 a + 1, DIMF, iters, dtime);

    // Copy psim into f[j][k]
    FORALLSITES(i, s) {
      for (b = 0; b < DIMF; b++)
        prop[i].e[a][b] = psim[0][i].c[b];
    }
  }

  
 
  c_scalar_mult_mat(&(Lambda2[5]),&id,&tmat);
  // Compute bilinear condensate

  for (a = 0; a < DIMF; a++) {
    for (b = 0; b < DIMF; b++)
  if (node_number(pnt[0], pnt[1], pnt[2], pnt[3]) == mynode()) {
    i = node_index(pnt[0], pnt[1], pnt[2], pnt[3]);
    bilin.real += parity(dumcoords)*(prop[i].e[a][b].real * tmat.e[a][b].real - prop[i].e[a][b].imag * tmat.e[a][b].imag);
    bilin.imag += parity(dumcoords)*(prop[i].e[a][b].imag * tmat.e[a][b].real + prop[i].e[a][b].real * tmat.e[a][b].imag);
  }
}
 
   
  g_complexsum(&bilin);
  
  
//-------------------------------------------------------------------------------------------------

  for (a = 0; a < DIMF; a++) {
    for (b = 0; b < DIMF; b++) {
      for(c = 0; c < DIMF; c++){
        for(d = 0; d < DIMF; d++){
           FORALLSITES(i, s) {
       // Two-point function
             t = (pnt[3] - (s->t) + nt) % nt;
       
             CMUL(prop[i].e[a][d],prop[i].e[b][c],tmp);
             CMUL(Lambda2[5].e[a][b],Lambda2[5].e[c][d],tmp2)      
             CMULREAL(tmp,parity(pnt)*parity(dumcoords),tmp);
             CMUL(tmp,tmp2,tmp3)
             CDIF(corr[t],tmp3);
             CMUL(prop[i].e[a][c],prop[i].e[b][d],tmp);
             CMUL(Lambda2[5].e[a][b],Lambda2[5].e[c][d],tmp2)
             CMULREAL(tmp,parity(pnt)*parity(dumcoords),tmp);
             CMUL(tmp,tmp2,tmp3)
             CSUM(corr[t],tmp3);
           }
         }
       }
     }
   }

        
  for (t = 0; t < nt; t++)
   g_complexsum(&(corr[t]));
        
        
  

  
  

 // Gather prop from -mu for one-link condensates and gauge-links from ( x -\mu)...temp_link[dir] represents curly link matrices in reduced staggered formalism
  for (dir = XUP; dir <= TUP; dir++) {
      tag[dir] = start_gather_field(prop, sizeof(matrix), dir,
                                    EVENANDODD, gen_pt[dir]);
      tag[OPP_DIR(dir)] = start_gather_field(prop, sizeof(matrix), OPP_DIR(dir),EVENANDODD, gen_pt[OPP_DIR(dir)]);

      mtag[dir] = start_gather_site(F_OFFSET(temp_link[dir]), sizeof(matrix),OPP_DIR(dir), EVENANDODD, gen_pt[11-dir]);
  }
  // Compute one-link condensates

 
  
  for (dir = XUP; dir <= TUP; dir++) {
  
     
    wait_gather(tag[dir]);
    wait_gather(tag[OPP_DIR(dir)]);
    wait_gather(mtag[dir]);
  if (node_number(pnt[0], pnt[1], pnt[2], pnt[3]) == mynode()) {
      i = node_index(pnt[0], pnt[1], pnt[2], pnt[3]);
    
       tm =  (matrix *)(gen_pt[dir][i]);
       
      
       BC=1;
       if (dir == TUP && pnt[TUP] == nt-1 && PBC < 0) {BC=-1;}
      
      for(a=0 ; a < DIMF ; a++){
      for(b=0 ; b < DIMF ; b++){
     
       if (dumcoords[dir]%2 == 0){
       
       CMULREAL_SUM(lattice[i].temp_link[dir].e[a][b], tm->e[a][b],(double) BC * lattice[i].phase[dir], one_link);
  
                 }
       
       else {

       CMULREAL_DIF(lattice[i].temp_link[dir].e[a][b], tm->e[a][b],(double) BC * lattice[i].phase[dir], one_link);
            }
       
       }

   }
       tm2 = (matrix *)(gen_pt[OPP_DIR(dir)][i]);
       tm3 = (matrix *)(gen_pt[11-dir][i]);

       BC=1;
       if(dir == TUP && pnt[TUP] == 0 && PBC < 0) {BC=-1;}

      for(a=0 ; a < DIMF ; a++){
      for(b=0 ; b < DIMF ; b++){

       if(dumcoords[dir]%2 == 0){

       CMULREAL_SUM(tm3->e[b][a], tm2->e[a][b], (double) BC  * lattice[i].phase[dir] , one_link); 
       
                   }
        
       else {

       CMULREAL_DIF(tm3->e[b][a], tm2->e[a][b], (double) BC * lattice[i].phase[dir] , one_link); 

            }    

          
       }
    
    
    }
     

 }

    cleanup_gather(tag[dir]);
    cleanup_gather(tag[OPP_DIR(dir)]);
    cleanup_gather(mtag[dir]);
}
  
  g_complexsum(&(one_link));


    
  
  // Print condensates
 
  
  node0_printf("PNT BILIN %d %d %d %d %.6g %.6g %d\n",
               pnt[0], pnt[1], pnt[2], pnt[3],
               bilin.real,bilin.imag,tot_iters);
  
  node0_printf("PNT ONELINK %d %d %d %d %.6g %.6g %d\n",
               pnt[0], pnt[1], pnt[2], pnt[3],
               one_link.real, one_link.imag,tot_iters);
  
  for (t = 0; t < nt; t++)
    node0_printf("CORR %d %lg %lg\n", t, 4*corr[t].real,4*corr[t].imag);



  // Reset multi-mass CG and clean up
  Norder = sav;
  free(psim[0]);
  free(psim);
  return tot_iters;
}
// -----------------------------------------------------------------
