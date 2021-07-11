#define R_NO_REMAP
#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <R.h>
#include <Rinternals.h>

// Calculate global level associated with bounds given in b_vec
// Parameter bound_id indicates if bound is upper or lower
// Uses speedup for equal local levels case 
void jointlevel_twosided_ell_speedup(double *b_vec, int *bound_id, int *num_points, double *out) {

  int n = *num_points;
  double *b_vec_prev;
  double *lgamma_arr;
  // Pre-compute vector of lgamma values to avoid repeat computation
  lgamma_arr = (double *) malloc((n + 2) * sizeof(double));
  lgamma_arr[0] = 0;
  lgamma_arr[1] = 0;
  for(int i = 2; i <= (n + 1); i = i + 1) {

      lgamma_arr[i] = lgamma_arr[i - 1] + log(i - 1);

  }

  b_vec_prev = (double *) malloc((n + 1) * sizeof(double));
  b_vec_prev[0] = pow((1 - b_vec[0]), n);

  double *b_vec_next;
  b_vec_next = (double *) malloc((n + 1) * sizeof(double));

  int j_lower = 0;
  int j_upper = 0;
  int l_lower = 0;
  int l_upper = 0;

  for(int k = 1; k <= n; k = k + 1) {

    l_lower = j_lower;
    l_upper = j_upper;

    if (bound_id[k - 1] == 0) {

      j_upper = j_upper + 1;

    }

    if(bound_id[k] == 1) {

      // checks the current bound id
      j_lower = j_lower + 1;

    }

    for(int j = j_lower; j <= j_upper; j = j + 1) {

      b_vec_next[j] = 0; // initialize the probability value
      int l_loop_upper_bound = MIN(l_upper, j);

      for(int l = l_lower; l <= l_loop_upper_bound; l = l + 1) {

        b_vec_next[j] = b_vec_next[j] + exp(log(b_vec_prev[l]) + (j - l) * log(b_vec[k] -
          b_vec[k - 1]) + (n - j) * log(1 - b_vec[k]) - (n - l) * log(1 - b_vec[k - 1]) +
          lgamma_arr[n - l + 1] - lgamma_arr[j - l + 1] - lgamma_arr[n - j + 1]);

      }

    }

    // Only copy over what is necessary
    if(k < n) memcpy(b_vec_prev + j_lower, b_vec_next + j_lower, (j_upper - j_lower + 1) * sizeof(double));

  }

  *out = 0;

  for(int l = l_lower; l <= l_upper; l = l + 1) {

    *out = *out + exp(log(b_vec_prev[l]) + log(b_vec_next[n - l]) + lgamma_arr[l + 1] + 
      lgamma_arr[n - l + 1] - lgamma_arr[n + 1] - l * log(b_vec[n - 1]) - 
      (n - l) * log(1 - b_vec[n - 1]));

  }

  free(b_vec_prev);
  free(b_vec_next);
  free(lgamma_arr);

}

