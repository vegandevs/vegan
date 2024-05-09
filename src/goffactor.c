/* Helper functions used in envfit for assessing goodness of fit in
 * permuation tests. The functions are really minimal, but they are
 * the bottlenecks that need be hidden in C code. 
 */

/*
 * Goodness of fit for a fitted factor: Sum of weighted variances
 * within factor levels. 
 */

#include <R.h>

static void goffactor(double *ord, int *f, double *w, int *nrow, int *ndim,
		      int *nlev,  double *sw, double *swx, double *swxx,
		      double *var)
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

void wcentre(double *x, double *retx, double *w, int *nr, int *nc)
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
	       retx[ij] = x[ij] - swx;
	       retx[ij] *= sqrt(w[i]);
	  }
     }
}

/* .Call interfaces */

#include <Rinternals.h>

SEXP do_wcentre(SEXP x, SEXP w)
{
    int nr = nrows(x), nc = ncols(x);
    if (length(w) != nr)
	error("weights 'w' and data do not match");
    if (TYPEOF(x) != REALSXP)
	x  = coerceVector(x, REALSXP);
    SEXP rx = PROTECT(duplicate(x));
    if (TYPEOF(w) != REALSXP)
	w = coerceVector(w, REALSXP);
    w = PROTECT(duplicate(w));
    SEXP retx = PROTECT(allocMatrix(REALSXP, nr, nc));
    wcentre(REAL(rx), REAL(retx), REAL(w), &nr, &nc);
    /* set dimnames */
    SEXP dnames = getAttrib(x, R_DimNamesSymbol);
    if (!isNull(dnames)) {
        SEXP dimnames = PROTECT(allocVector(VECSXP, 2));
        SET_VECTOR_ELT(dimnames, 0, duplicate(VECTOR_ELT(dnames, 0)));
        SET_VECTOR_ELT(dimnames, 1, duplicate(VECTOR_ELT(dnames, 1)));
        setAttrib(retx, R_DimNamesSymbol, dimnames);
        UNPROTECT(1);
    }
    UNPROTECT(3);
    return retx;
}

SEXP do_goffactor(SEXP x, SEXP factor, SEXP nlevels, SEXP w)
{
    int i, nr = nrows(x), nc = ncols(x), nl = asInteger(nlevels);
    /* check dimensions (or segfault) */
    if (length(factor) != nr)
	error("dimensions of data and factor do not match");
    if (length(w) != nr)
	error("dimensions of data and weights (w) do not match");
    /* return variance */
    SEXP var = PROTECT(allocVector(REALSXP, 1));
    if (TYPEOF(factor) != INTSXP)
	factor = coerceVector(factor, INTSXP);
    PROTECT(factor);
    SEXP fac1 = PROTECT(duplicate(factor));
    if (TYPEOF(x) != REALSXP)
	x = coerceVector(x, REALSXP);
    PROTECT(x);
    if (TYPEOF(w) != REALSXP)
	w = coerceVector(w, REALSXP);
    PROTECT(w);
    /* goffactor assume factors start from 0 */
    for (i = 0; i < nr; i++)
	INTEGER(fac1)[i]--;
    double *work1 = (double *) R_alloc(nl, sizeof(double));
    double *work2 = (double *) R_alloc(nl, sizeof(double));
    double *work3 = (double *) R_alloc(nl, sizeof(double));
    goffactor(REAL(x), INTEGER(fac1), REAL(w), &nr, &nc, &nl,
	      work1, work2, work3, REAL(var));
    UNPROTECT(5);
    return var;
}
