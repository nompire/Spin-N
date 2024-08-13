// -----------------------------------------------------------------
// Measure total action, as needed by the hybrid Monte Carlo algorithm
// When this routine is called the CG should already have been run,
// so that the vector **sol contains (M_adjoint*M+shift[n])^(-1) * src
#include "sp_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------

// -----------------------------------------------------------------
double hmom_action(){
register int i, mu;
  register site *s;
  matrix tmat;
  complex tr;
  double sum = 0.0;
  double sum_sig = 0.0;
  for(mu = XUP ; mu <=TUP ; mu ++){
  FORALLSITES(i, s) {
      mult_mat_nn(&(s->mom[mu]),&(s->mom[mu]),&tmat);
      tr = trace(&tmat);
      sum += -0.5 * tr.real;
        }
  }
  g_doublesum(&sum);
  
  //--------------------Scalar momenta contribution------------------
  FORALLSITES(i, s) {
          mult_mat_an(&(s->p_sigma),&(s->p_sigma),&tmat);
                tr = trace(&tmat);
                      sum_sig += 0.5 * tr.real;
                        }
  
  
   g_doublesum(&sum_sig);
  // -----------------------------------------------------------------
  
  
  
  
  
  
  
  
  
  
  
  
  return (sum + sum_sig);
}



// -----------------------------------------------------------------




// -----------------------------------------------------------------
// Fermion contribution to the action
// Include the ampdeg term to allow sanity check that the fermion action
// is 4*volume on average
// Since the pseudofermion src is fixed throughout the trajectory,
// ampdeg actually has no effect on Delta S (checked)
// sol, however, depends on the gauge fields through the CG
double fermion_action(vector *src, vector **sol) {
  register int i, j;
  register site *s;
  double sum = 0.0;
  double tr;
#ifdef DEBUG_CHECK
  double im = 0.0;
#endif

  FORALLSITES(i, s) {
    sum += ampdeg4 * magsq_vec(&(src[i]));
    for (j = 0; j < Norder; j++) {
      tr = dot(&(src[i]), &(sol[j][i]));   // src^dag.sol[j]
      sum += (amp4[j] * tr);
    }
  }
  g_doublesum(&sum);
  
  return sum;
}
// -----------------------------------------------------------------


// Adds adjoint plaquette term
// Use tempmat for temporary storage
void plaquette_a(double *ss_plaq, double *st_plaq) {
  register int i, dir, dir2;
  register site *s;
  register matrix *m1, *m4;
  
  double ss_sum = 0.0, st_sum = 0.0, td;
  complex tr;
  msg_tag *mtag, *mtag2;
  matrix tmat;
  
FORALLSITES(i, s) {
  clear_mat(&(tempmat[i]));
}
 
  
  for (dir = YUP; dir <= TUP; dir++) {
    for (dir2 = XUP; dir2 < dir; dir2++) {
      mtag = start_gather_site(F_OFFSET(link[dir2]), sizeof(matrix),
                               dir, EVENANDODD, gen_pt[0]);
      mtag2 = start_gather_site(F_OFFSET(link[dir]), sizeof(matrix),
                                dir2, EVENANDODD, gen_pt[1]);

      FORALLSITES(i, s) {
        m1 = &(s->link[dir]);
        m4 = &(s->link[dir2]);
        
        mult_mat_an(m4, m1, &(tempmat[i]));
      }
      

      wait_gather(mtag);
      wait_gather(mtag2);
      FORALLSITES(i, s) {
        m1 = (matrix *)(gen_pt[0][i]);
        m4 = (matrix *)(gen_pt[1][i]);
        
        mult_mat_nn(&(tempmat[i]), m1, &tmat);
        tr = complextrace(m4, &tmat);
        td = tr.real;  
                                  
                     
        
        
        if (dir == TUP)
          st_sum += td;
          
        else
          ss_sum += td;
      }
      cleanup_gather(mtag);
      cleanup_gather(mtag2);
      
    }
  }
  
  g_doublesum(&ss_sum);
  g_doublesum(&st_sum);
  
  
  *ss_plaq = ss_sum / ((double)volume);
  *st_plaq = st_sum / ((double)volume);

   
}

double scalar_action(){

	register int i;
	register site *s;
	double sum = 0.0;
	matrix tmat;
	complex tr = cmplx(0.0,0.0);
	FORALLSITES(i,s) {

		  mult_mat_an(&(s->sigma),&(s->sigma),&tmat);
		  tr = trace(&tmat);

		  sum += 0.5*tr.real;
	}
	g_doublesum(&sum);

	return(sum);
}



// -----------------------------------------------------------------
// Print out total action and individual contributions
double action(vector **src, vector ***sol) {
  
  int n;

  double h_act, f_act, g_action, s_action;
  double total = 0.0;
  double ssplaq=0.0;
  double stplaq=0.0;
  plaquette_a(&ssplaq,&stplaq); 
  
 // Fermion action
  for (n = 0; n < Nroot; n++) {
    f_act = fermion_action(src[n], sol[n]);
    node0_printf("fermion%d %.8g ", n, f_act);
    //f_act = 0.0;
    total += f_act;
  }

  h_act = hmom_action(); // Calculates the momentum contribution to the Hamiltonian
  g_action = -1.0 * (BETA /4.0) * volume * (ssplaq + stplaq); // Gauge Action
  s_action = scalar_action();  
  node0_printf("momenta %.8g ", h_act);
  node0_printf("plaquette %.8g ", g_action);
  node0_printf("scalar action %.8g",s_action);
  total += h_act;
  total += g_action;
  total += s_action;
  node0_printf("sum %.8g\n", total);

  
  node0_printf("Wilson Plaquette %.8g\n", (ssplaq + stplaq)/(DIMF));
  node0_printf("\n");
 
  
  return total;
}
// -----------------------------------------------------------------
