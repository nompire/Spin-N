// -----------------------------------------------------------------
// Defines and subroutine declarations
// for the  reduced staggered-fermion system with SP(N) gauge symmetry
#ifndef _SP_H
#define _SP_H

#include "../include/random.h"
#include "../include/complex.h"

// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Fermions are SP(N) spinors
//
#define DIMF 4
#define GAUGETOL 0.00000000000001

typedef struct {fcomplex e[DIMF][DIMF];} fmatrix;
typedef struct {dcomplex e[DIMF][DIMF];} dmatrix;
typedef struct {
  fcomplex m01;
  float m00im, m11im;
  float space;
} fanti_hermitmat;
typedef struct {
  dcomplex m01;
  double m00im, m11im;
  double space;
} danti_hermitmat;

typedef struct { fcomplex c[DIMF]; } fvector;

typedef struct { dcomplex c[DIMF]; } dvector;

#if (PRECISION == 1)
#define matrix     fmatrix
#define anti_hermitmat fanti_hermitmat
#define vector         fvector

#else
#define matrix  dmatrix
#define anti_hermitmat danti_hermitmat
#define vector         dvector
#endif

#define PLUS 1          // Flags for selecting D or D_adjoint
#define MINUS -1
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Subroutine definitions
// Vector operations
// In file clearvec.c
void clearvec(vector *s);

// In file vec_copy.c
void vec_copy(vector *a, vector *b);

// In file dumpvec.c
void dumpvec(vector *v);

// In file subvec.c
void sub_vec(vector *a, vector *b, vector *c);

void add_vec(vector *a, vector *b, vector *c);

// In file msq_vec.c
double magsq_vec(vector *a);

// In file dot.c
 
void proj(vector *a,vector *b, matrix *c );

double dot(vector *a, vector *b);

complex cmplx_dot(vector *a,vector *b );


complex link_op(matrix *a , matrix *b);


// In file s_m_vec.c
void scalar_mult_vec(vector *src, double s, vector *dest);

// In file s_m_a_vec.c
void scalar_mult_add_vec(vector *a, vector *b, double s, vector *c);
// -----------------------------------------------------------------

void scalar_mult_sum_matrix(matrix *a,double s,
		        matrix *c);

// Matrix operations
// In file clear_mat.c
void clear_mat(matrix *m);
//In file dumpmat.c
void dumpmat(matrix *m);

void add_matrix(matrix *a, matrix *b, matrix *c );
void add_mat(matrix *a,matrix *b);
void sub_mat(matrix *a,matrix *b);
void sub_matrix(matrix *a,matrix *b,matrix *c );


//-------------------- SP(4) generators----------------------------

#define NUMGEN 10
#define NUMYUK 1
void setup_gamma();


//------------------------------------------------------------------

// ROUTINES FOR SP(N) MATRIX OPERATIONS

void mult_mat_nn(  matrix *a, matrix *b, matrix *c  );

void mult_mat_na( matrix *a, matrix *b, matrix *c );


void mult_mat_an( matrix *a, matrix *b, matrix *c );

void mult_mat_cn( matrix *a, matrix *b, matrix *c );

double realtrace(  matrix *a, matrix *b );

//void trace_sp( matrix *a,complex s);
//* file "trace_sp.c"
complex  complextrace( matrix *a, matrix *b);


void scalar_mult_matrix( matrix *a, double s, matrix *b );

void scalar_mult_add_matrix( matrix *a, matrix *b,double s, matrix *c);

void c_scalar_mult_sub_mat( matrix *m1, matrix *m2, complex *phase, matrix *m3);
// file "cs_m_s_mat.c"

void c_scalar_mult_add_mat( matrix *m1, matrix *m2, complex *phase, matrix *m3);

void c_scalar_mult_mat( matrix *b, complex *s, matrix *c );

void adjoint( matrix *a, matrix *b );
// file "adjoint.c"

void make_anti_hermitian( matrix *m3,  anti_hermitmat *ah3 );
// file "make_ahmat.c"

void random_anti_hermitian( anti_hermitmat *mat_antihermit, double_prn *prn_pt );
// (prn_pt passed through to myrand())
// file "rand_ahmat.c"

void uncompress_anti_hermitian( anti_hermitmat *mat_anti, matrix *mat );


void mat_copy(matrix *a,matrix *b );


complex trace( matrix *a);


void conjug( matrix *a, matrix *b );

void trans(matrix *a, matrix *b);

void det(matrix *a);

void exp_matrix(matrix *u,matrix *v);

void unit_mat( matrix *dest );


//-------------------------------------------------------------------------------------------//

// ROUTINES FOR vector OPERATIONS ( 4 COMPONENT COMPLEX )

void scalar_mult_vec( vector *a, double s, vector *c);

void projector(vector *u,vector *v, matrix *c,int flag);

void mult_mat_vec( matrix *a, vector *b, vector *c );

void mult_mat_vec_sum( matrix *a, vector *b, vector *c );


void mult_mat_vec_sum_4dir( matrix *a, vector *b0,
 vector *b1, vector *b2, vector *b3, vector *c );


void vec_conjug(vector *v,vector *u);


//---------------------------------------------------------------------------------------------//



// -----------------------------------------------------------------
// Miscellaneous routines
// In file gaussrand.c
Real gaussian_rand_no(double_prn *prn_pt);

#include "../include/int32type.h"
void byterevn(int32type w[], int n);
void byterevn64(int32type w[], int n);

#endif
// -----------------------------------------------------------------
