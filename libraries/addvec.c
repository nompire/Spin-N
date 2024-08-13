// -----------------------------------------------------------------
// Add two vectors
// c <-- a + b
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/sp.h"

void add_vec(vector *a, vector *b, vector *c) {

  CADD(a->c[0],b->c[0],c->c[0]);
  CADD(a->c[1],b->c[1],c->c[1]);
  CADD(a->c[2],b->c[2],c->c[2]);
  CADD(a->c[3],b->c[3],c->c[3]);
  
}
// -----------------------------------------------------------------
