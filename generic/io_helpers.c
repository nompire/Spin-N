// -----------------------------------------------------------------
// General purpose IO-related routines
#include "generic_includes.h"
#include "../include/io_lat.h"
#include "../RHMC/lattice.h"
#include "../include/sp.h"
#include <stdlib.h>
// -----------------------------------------------------------------



// -----------------------------------------------------------------

// -----------------------------------------------------------------



// -----------------------------------------------------------------
gauge_file *save_lattice(int flag, char *filename) {
  double dtime;
  gauge_file *gf = NULL;

  //d_sigmasum(&sigmasum);

  dtime = -dclock();
  switch(flag) {
    case SAVE_SERIAL:
      gf = save_serial(filename);
      break;
    case FORGET:
      gf = NULL;
      break;
    default:
      node0_printf("\nsave_lattice: ERROR: unknown type for saving lattice\n");
      terminate(1);
  }
  dtime += dclock();
  if (flag != FORGET)
    node0_printf("Time to save = %e\n", dtime);
//#if PRECISION == 1
  //node0_printf("CHECK SIGMA SUM: %e\n", sigmasum);
//#else             // Double precision
  //node0_printf("CHECK SIGMA SUM: %.16e\n", sigmasum);
//#endif
  return gf;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Set sigma to unity
void coldlat() {
  
  register int i;
  register int dir;
  register site *s;
  register int k,l,m;
  FORALLSITES(i, s) {
  for(dir= XUP; dir <= TUP ; dir++){
  for(k=0; k < DIMF ; k++){
  for(l=0; l < DIMF ; l++){
   
   s->link[dir].e[k][l] =cmplx(0.0,0.0);
                    }
         }

   for(m=0 ;m < DIMF ; m++){
   
   s->link[dir].e[m][m] = cmplx(1.0,0.0);

     }
  
  }
 
 
 
 
}
    FORALLSITES(i,s){

	      clear_mat(&(s->sigma));
	      for(m=0 ; m < NUMYUK ; m++){
	      #ifdef SITERAND

	            scalar_mult_add_matrix(&(s->sigma),&(Lambda2[m]),gaussian_rand_no(&(s->site_prn)),&(s->sigma));
	      #else

	            scalar_mult_add_matrix(&(s->sigma),&(Lambda2[m]), gaussian_rand_no(&node_prn),&(s->sigma));
	      #endif

	       }

    }



  node0_printf("unit scalar configuration loaded\n");
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Set sigma to random gaussian numbers
void randomlat() {
  
  register int i;
  register site *s;
  register int j,m ;
  FORALLSITES(i, s) {
#ifdef SITERAND
   
  for (j=0  ;  j< 4  ; j++){
  for (m=0  ;  m < NUMGEN ; m++){
  
     
    clear_mat(&(s->temp_link[j]));
    scalar_mult_add_matrix(&(s->temp_link[j]),&(Lambda[m]),  0.5 * gaussian_rand_no(&(s->site_prn)) , &(s->temp_link[j])); 
   }
   exp_matrix(&(s->temp_link[j]),&(s->link[j]));
   
}
#else
  for(j=0 ; j < 4 ; j++){
  for(m=0 ; m < NUMGEN ; m++){
  
    clear_mat(&(s->temp_link[j]));
    scalar_mult_add_matrix(&(s->temp_link[j]), &(Lambda[m]),  0.5 * gaussian_rand_no(&node_prn) , &(s->temp_link[j]));
  }
   exp_matrix(&(s->temp_link[j]),&(s->link[j]));
   
}
#endif
   
}

 FORALLSITES(i,s){

	   clear_mat(&(s->sigma));
	   for(m=0 ; m < NUMYUK ; m++){
           #ifdef SITERAND

	      scalar_mult_add_matrix(&(s->sigma),&(Lambda2[m]), gaussian_rand_no(&(s->site_prn)),&(s->sigma));
	   #else

	      scalar_mult_add_matrix(&(s->sigma),&(Lambda2[m]),gaussian_rand_no(&node_prn),&(s->sigma));
	   #endif

	   }

 }


node0_printf("random gauge configuration loaded\n");

}
// -----------------------------------------------------------------






// -----------------------------------------------------------------
// Reload a lattice in binary format, set to unity or keep current
gauge_file *reload_lattice(int flag, char *filename) {
  double dtime;
  gauge_file *gf = NULL;

  dtime = -dclock();
  switch(flag) {
    case CONTINUE:        // Do nothing
      gf = NULL;
      break;
    case FRESH:           // Cold lattice
      coldlat();
      gf = NULL;
      break;
    case RANDOM:          // Random (hot) lattice
      randomlat();
      gf = NULL;
      break;
    case RELOAD_SERIAL:   // Read binary lattice serially
      gf = restore_serial(filename);
      break;
    default:
      node0_printf("reload_lattice: Bad startflag %d\n", flag);
      terminate(1);
  }
  dtime += dclock();
  if (flag != FRESH && flag != RANDOM && flag != CONTINUE)
    node0_printf("Time to reload gauge configuration = %e\n", dtime);

  //d_sigmasum(&sigmasum);
//#if PRECISION == 1
//  node0_printf("CHECK SIGMA SUM: %e\n", sigmasum);
//#else             // Double precision
//  node0_printf("CHECK SIGMA SUM: %.16e\n", sigmasum);
//#endif
  fflush(stdout);
  dtime = -dclock();
  return gf;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Find out kind of starting lattice to use, and lattice name if necessary
// This routine is only called by node 0
int ask_starting_lattice(FILE *fp, int prompt, int *flag, char *filename) {
  char savebuf[256];
  int status;

  if (prompt!=0)
    printf("enter 'continue', 'fresh', 'random' or 'reload_serial'\n");
  status = fscanf(fp, "%s", savebuf);
  if (status == EOF) {
    printf("ask_starting_lattice: EOF on STDIN.\n");
    return 1;
  }
  if (status != 1) {
    printf("\nask_starting_lattice: ERROR IN INPUT: ");
    printf("can't read starting lattice option\n");
    return 1;
  }

  printf("%s", savebuf);
  if (strcmp("fresh", savebuf) == 0) {
    *flag = FRESH;
    printf("\n");
  }
  else if (strcmp("random", savebuf) == 0) {
    *flag = RANDOM;
    printf("\n");
  }
  else if (strcmp("continue", savebuf) == 0) {
    *flag = CONTINUE;
    printf("\n");
  }
  else if (strcmp("reload_serial", savebuf) == 0)
    *flag = RELOAD_SERIAL;
  else {
    printf(" is not a valid starting lattice option. INPUT ERROR.\n");
    return 1;
  }

  // Read name of file and load it
  if (*flag != FRESH && *flag != RANDOM && *flag != CONTINUE) {
    if (prompt != 0)
      printf("enter name of file containing lattice\n");
    status = fscanf(fp, " %s", filename);
    if (status != 1) {
      printf("\nask_starting_lattice: ERROR IN INPUT: ");
      printf("error reading file name\n");
      return 1;
    }
    printf(" %s\n", filename);
  }
  return 0;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Find out what do to with lattice at end, and lattice name if necessary
// This routine is only called by node 0
int ask_ending_lattice(FILE *fp, int prompt, int *flag, char *filename) {
  char savebuf[256];
  int status;

  if (prompt!=0)
    printf("'forget' lattice at end or 'save_serial'\n");
  status = fscanf(fp,"%s", savebuf);
  if (status != 1) {
    printf("\nask_ending_lattice: ERROR IN INPUT: error reading ending lattice command\n");
    return 1;
  }
  printf("%s", savebuf);
  if (strcmp("save_serial", savebuf) == 0)
    *flag = SAVE_SERIAL;
  else if (strcmp("forget", savebuf) == 0) {
    *flag = FORGET;
    printf("\n");
  }
  else {
    printf(" is not a save lattice command. INPUT ERROR\n");
    return 1;
  }

  if (*flag != FORGET) {
    if (prompt != 0)
      printf("enter filename\n");
    status = fscanf(fp, "%s", filename);
    if (status != 1) {
      printf("\nask_ending_lattice: ERROR IN INPUT: ");
      printf("error reading filename\n");
      return 1;
    }
    printf(" %s\n", filename);
  }
  return 0;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Read and echo the next tag.  Echo any intervening comments
// Comments begin with # and apply to the rest of the line
// Verify that the input tag agrees with the expected tag
static int get_tag(FILE *fp, char *tag, char *myname) {
  static char checktag[80];
  char line[512];
  int s;

  while (1) {
    s = fscanf(fp, "%s", checktag);
    if (s == EOF) {
      printf("%s(%d): EOF on input\n", myname, this_node);
      return 1;
    }
    if (s == 0) {
      printf("%s(%d) Error reading %s\n", myname, this_node, tag);
      return 1;
    }
    if (strchr(checktag, '#') != NULL) {
      printf("%s", checktag);
      if (fgets(line, 512, fp) == NULL) {
        printf("%s(%d) EOF on input.\n", myname, this_node);
        return 1;
      }
      printf("%s", line);
    }
    else {
      if (strcmp(checktag, tag) != 0) {
        printf("\n%s: ERROR IN INPUT: expected %s but found %s\n",
               myname, tag, checktag);
        return 1;
      }
      printf("%s ", tag);
      return 0;
    }
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Check return value of scanf
static int check_read(int s, char *myname, char *tag) {
  if (s == EOF) {
    printf("\n%s: Expecting value for %s but found EOF\n", myname, tag);
    return 1;
  }
  else if (s == 0) {
    printf("\n%s: Format error reading value for %s\n", myname, tag);
    return 1;
  }
  else
    return 0;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Get a floating point number
// If prompt is non-zero, ask for the input value with tag
// If prompt is zero, require that tag precede the input value
// Return the value and exit on error
int get_f(FILE *fp, int prompt, char *tag, Real *value) {
  int s;
  char checkvalue[80];
  char myname[] = "get_f";

  if (prompt) {
    s = 0;
    while (s != 1) {
      printf("enter %s ", tag);
      (void)fscanf(fp,"%s", checkvalue);
#if PRECISION == 1
      s = sscanf(checkvalue, "%e", value);
#else
      s = sscanf(checkvalue, "%le", value);
#endif
      if (s == EOF)
        return 1;
      if (s == 0)
        printf("Data format error.\n");
      else printf("%s %g\n", tag, *value);
    }
  }
  else {
    if (get_tag(fp, tag, myname) == 1)
      return 1;

#if PRECISION == 1
    s = fscanf(fp, "%e", value);
#else
    s = fscanf(fp, "%le", value);
#endif
    if (check_read(s, myname, tag) == 1)
      return 1;

    printf("%g\n", *value);
  }
  return 0;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Get an integer with same behavior as get_f
int get_i(FILE *fp, int prompt, char *tag, int *value) {
  int s;
  char checkvalue[80];
  char myname[] = "get_i";

  if (prompt) {
    s = 0;
    while (s != 1) {
      printf("enter %s ", tag);
      (void)fscanf(fp, "%s", checkvalue);
      s=sscanf(checkvalue,"%d", value);
      if (s == EOF) return 1;
      if (s == 0) printf("Data format error\n");
      else printf("%s %d\n", tag, *value);
    }
  }
  else {
    if (get_tag(fp, tag, myname) == 1)
      return 1;

    s = fscanf(fp, "%d", value);
    if (check_read(s, myname, tag) == 1)
      return 1;
    printf("%d\n", *value);
  }
  return 0;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Read a single word as a string with same behavior as get_f
int get_s(FILE *fp, int prompt, char *tag, char *value) {
  int s;
  char myname[] = "get_s";

  if (prompt) {
    s = 0;
    while (s != 1) {
      printf("enter %s ", tag);
      s = fscanf(fp, "%s", value);
      if (s == EOF) return 1;
      if (s == 0) printf("Data format error\n");
      else printf("%s %s\n", tag, value);
    }
  }
  else {
    if (get_tag(fp, tag, myname) == 1) return 1;

    s = fscanf(fp, "%s", value);
    if (check_read(s, myname, tag) == 1) return 1;
    printf("%s\n", value);
  }
  return 0;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Read a vector of integers with same behavior as get_f
int get_vi(FILE* fp, int prompt, char *tag, int *value, int nvalues) {
  int s, i;
  char myname[] = "get_vi";

  if (prompt) {
    s = 0;
    printf("enter %s with %d values", tag, nvalues);
    for (i = 0; i < nvalues; i++) {
      while (s != 1) {
        printf("\n[%d] ", i);
        s = fscanf(fp, "%d", value + i);
        if (s == EOF) return 1;
        if (s == 0) printf("Data format error\n");
        printf("%s %d\n", tag, value[i]);
      }
    }
  }
  else {
    if (get_tag(fp, tag, myname) == 1) return 1;

    for (i = 0; i < nvalues - 1; i++) {
      s = fscanf(fp, "%d", value + i);
      if (check_read(s, myname, tag) == 1) return 1;
      printf("%d ", value[i]);
    }
    s = fscanf(fp, "%d", value + nvalues - 1);
    if (check_read(s, myname, tag) == 1) return 1;
    printf("%d\n", value[nvalues - 1]);
  }
  return 0;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Read a vector of reals with same behavior as get_f
int get_vf(FILE* fp, int prompt, char *tag, Real *value, int nvalues) {
  int s, i;
  char myname[] = "get_vf";

  if (prompt) {
    s = 0;
    printf("enter %s with %d values", tag, nvalues);
    for (i = 0; i < nvalues; i++) {
      while (s != 1) {
        printf("\n[%d] ", i);
#if PRECISION == 1
        s = scanf("%e", value + i);
#else
        s = scanf("%le", value + i);
#endif
        if (s == EOF)
          return 1;
        if (s == 0)
          printf("Data format error\n");

        printf("%s %g\n", tag, *(value + i));
      }
    }
  }
  else {
    if (get_tag(fp, tag, myname) == 1)
      return 1;

    for (i = 0; i < nvalues; i++) {
#if PRECISION == 1
      s = fscanf(fp, "%e", value + i);
#else
      s = fscanf(fp, "%le", value + i);
#endif
      if (check_read(s, myname, tag) == 1)
        return 1;
      printf("%g ", value[i]);
    }
    printf("\n");
  }
  return 0;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Get the initial value of prompt:
// 0 for reading from file, 1 prompts for input from terminal
// Should be called only by node 0
// Return 0 if successful, 1 if failure
int get_prompt(FILE *fp, int *prompt) {
  char initial_prompt[512];
  int status;
  char myname[] = "get_prompt";

  *prompt = -1;
  printf("type 0 for no prompts or 1 for prompts\n");
  while (1) {
    status = fscanf(fp, "%s", initial_prompt);
    if (status != 1) {
      printf("\n%s: Can't read input\n", myname);
      terminate(1);
    }
    if (strchr(initial_prompt,'#') == NULL) break;
    // Provide for comment lines with # before "prompt"
    else {
      printf("%s", initial_prompt);
      if (fgets(initial_prompt, 512, fp) == NULL) {
        printf("%s(%d) EOF on input.\n", myname, this_node);
        return 1;
      }
      printf("%s", initial_prompt);
    }
  }
  if (strcmp(initial_prompt, "prompt") == 0)
    (void)fscanf(fp, "%d", prompt);
  else if (strcmp(initial_prompt, "0") == 0)
    *prompt = 0;
  else if (strcmp(initial_prompt, "1") == 0)
    *prompt = 1;

  if (*prompt == 0 || *prompt == 1)
    return 0;
  else {
    printf("\n%s: ERROR IN INPUT: initial prompt\n", myname);
    return 1;
  }
}
// -----------------------------------------------------------------
