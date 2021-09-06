#define R_NO_REMAP

#include <math.h>
#include <stdlib.h>
#include <R.h>
#include <Rinternals.h>

// Calculate global level associated with one-sided bounds given in crit_vals
void jointlevel_onesided(double *crit_vals, int *num_points, int *checkint,
                         int *firstcheck, double *relerr, double *out){

  // Allocate all memory for cc and cn
  int n = *num_points;

  double *cc;
  cc = (double *) malloc((n + 1) * sizeof(double));

  double *cn;
  cn = (double *) malloc((n + 1) * sizeof(double));

  cc[0] = pow((1 - crit_vals[1]), n);
  cc[1] = exp(log(n) + log(crit_vals[1] - crit_vals[0]) + (n - 1) * log(1 - crit_vals[1]));
  
  // Pre-compute vector of lgamma values to avoid repeat computation
  double *lgamma_arr;
  lgamma_arr = (double *) malloc((n + 2) * sizeof(double));
  lgamma_arr[0] = 0;
  lgamma_arr[1] = 0;
  for(int i = 2; i <= (n + 1); i = i + 1) {
    
    lgamma_arr[i] = lgamma_arr[i - 1] + log(i - 1);
    
  }

  memcpy(cn, cc, 2 * sizeof(double));

  int skip = -1;
  int mem_skip = 0; // This value will cary around the number of values that are being skipped to make memcopys more efficient
  double accerr = 0;

  for(int k = 2; k <= (n - 1); k = k + 1) {

    if (skip == -1) {

      cn[0] = pow(1 - crit_vals[k], n);

    } else {

      cn[skip + 1] = exp(log(cc[skip + 1]) + (n - (skip + 1)) *
        log(1 - crit_vals[k]) - (n - (skip + 1)) * log(1 - crit_vals[k - 1]));

    }

    if (k >= (skip + 3)) {

      for(int j = (skip + 2); j <= (k - 1); j = j + 1){

        cn[j] = exp(log(cc[j]) + (n - j) * log(1 - crit_vals[k]) - (n - j) * log(1 - crit_vals[k - 1]));

        for(int l = (skip + 1); l <= (j - 1); l = l + 1) {

          cn[j] = cn[j] + exp(log(cc[l]) + (j - l) * log(crit_vals[k] - crit_vals[k - 1]) +
            (n - j) * log(1 - crit_vals[k]) - (n - l) * log(1 - crit_vals[k - 1]) +
            lgamma_arr[n - l + 1] - lgamma_arr[j - l + 1] - lgamma_arr[n - j + 1]);

        }

      }

    }

    cn[k] = 0;
    for(int l = (skip + 1); l <= (k - 1); l = l + 1) {

      cn[k] = cn[k] + exp(log(cc[l]) + (k - l) * log(crit_vals[k] - crit_vals[k - 1]) +
        (n - k) * log(1 - crit_vals[k]) - (n - l) * log(1 - crit_vals[k - 1]) +
        lgamma_arr[n - l + 1] - lgamma_arr[k - l + 1] - lgamma_arr[n - k + 1]);

    }
    

    memcpy(cc + mem_skip, cn + mem_skip, (k + 1 - mem_skip) * sizeof(double));

    if((k > *firstcheck && (k % *checkint) == 0 && k < n-1) || (k == *firstcheck)) {

      double tmp1 = *relerr - exp(log(1 + *relerr) + log(accerr));
      double possible_error = 0;

      for(int i = (skip + 1); i <= k; i = i + 1) {

        possible_error = possible_error + cc[i];

      }

      tmp1 = tmp1 - exp(log(*relerr) + log(possible_error));

      double tmp2 = cc[skip + 1];
      int tmp3 = skip + 1;

      while(tmp2 <= tmp1){

        tmp3 = tmp3 + 1;
        tmp2 = tmp2 + cc[tmp3];

      }

      accerr = accerr + tmp2 - cc[tmp3];
      skip = tmp3 - 1;
      if(skip > -1) {

        mem_skip = skip;

      }

    }

  }

  cn[n] = 0;
  for(int l = (skip + 1); l <= (n - 1); l = l + 1) {

    cn[n] = cn[n] + cc[l];

  }

  *out = 1 - cn[n];
  free(cc);
  free(cn);
  free(lgamma_arr);

}
