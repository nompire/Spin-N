// -----------------------------------------------------------------
// Set up generator matrices and epsilon^{ijklm}
#include "sp_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Set up four-index totally anti-symmetric tensor
// Initialize swap to avoid optimization dependence!!!
Real order(int i, int j, int k, int l) {
  int seq[4] = {i, j, k, l};
  int swap = 1, tmp, p, permutation = 1;
  while (swap > 0) {
    swap = 0;
    for (p = 0; p < 3; p++) {
      if (seq[p] > seq[p + 1]) {
        tmp = seq[p];
        seq[p] = seq[p + 1];
        seq[p + 1] = tmp;
        swap++;
        permutation *= -1;
      }
    }
  }
  return (Real)permutation;
}


void epsilon() {
  int i, j, k, l;
  
  for (i = 0; i < DIMF; i++) {
    for (j = 0; j < DIMF; j++) {
      for (k = 0; k < DIMF; k++) {
        for (l = 0; l < DIMF; l++)
          perm[i][j][k][l] = 0;
      }
    }
  }

  for (i = 0; i < DIMF; i++) {
    for (j = 0; j < DIMF; j++) {
      if (j == i)
        continue;
      for (k = 0; k < DIMF; k++) {
        if (k == j || k == i)
          continue;
        for (l = 0; l < DIMF; l++) {
          if (l == k || l == j || l == i)
            continue;
          perm[i][j][k][l] = order(i, j, k, l);
#ifdef DEBUG_CHECK
          if (perm[i][j][k][l] * perm[i][j][k][l] > 1.0e-4)
            node0_printf("PERM %d%d%d%d = %.4g\n",
                         i, j, k, l, perm[i][j][k][l]);
#endif
        }
      }
    }
  }
  return;
}
// -----------------------------------------------------------------
