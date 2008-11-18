/* Helper functions used in envfit for assessing goodness of fit in
 * permuation tests. The functions are really minimal, but they are
 * the bottlenecks that need be hidden in C code. 
 */

/*
 * Goodness of fit for a fitted factor: Sum of weighted variances
 * within factor levels. 
 */

#include <R.h>

void goffactor(double *ord, int *f, double *w, int *nrow, int *ndim, int *nlev, 
	       double *sw, double *swx, double *swxx, double *var)
{
     int i, j, k;
     
     for (i = 0; i < (*nlev); i++) {
	  sw[i] = 0.0;
     }
      for (i = 0; i < (*nrow); i++) {
	  sw[f[i]] += w[i];
     }

     var[0] = 0.0;
     for (k = 0; k < (*ndim); k++) {
	  for (i = 0; i < (*nlev) ; i++) {
	       swx[i] = swxx[i] = 0.0;
	  }
	  for (i = k * (*nrow), j = 0; j < (*nrow); j++, i++) {
	       swx[f[j]] += w[j]*ord[i];
	       swxx[f[j]] += w[j]*ord[i]*ord[i];
	  }
	  for (j = 0; j < (*nlev); j++) {
	       if (sw[j] > 0)
		    var[0] += swxx[j] - swx[j]*swx[j]/sw[j];
	  }
     }
     return;
} 

/* Weighted centring of a matrix. Needed in many places in
 * vectorfit. R implementations has one apply and two sweeps, and is
 * probably much too slow. 
 */

#include <math.h> /* sqrt */

void wcentre(double *x, double *w, int *nr, int *nc)
{
     int i, j, ij;
     double sw, swx;

     for (i = 0, sw = 0.0; i < (*nr); i++)
	  sw += w[i];

     for (j = 0; j < (*nc) ; j++) {
	  for (i = 0, swx = 0.0, ij = (*nr)*j; i < (*nr); i++, ij++) {
	       swx += w[i] * x[ij];
	  }
	  swx /= sw;
	  for (i = 0,  ij = (*nr)*j; i < (*nr); i++, ij++) {
	       x[ij] -= swx;
	       x[ij] *= sqrt(w[i]);
	  }
     }
}
