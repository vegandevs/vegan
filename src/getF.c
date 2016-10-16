/* Function to evaluate F-value in permutest.cca. For instance, in R
 * CMD check in Macbook Air this function uses 1/4 of computing time,
 * and in applications with constrained ordination this is the major
 * function. Even small speed-up in this function will have a
 * considerable impact in running time.
 */


/* Function called form getF to evaluate the sum of all eigenvalues */

double getEV(double *x, int nr, int nc, int isDB)
{
    int i, ii;
    double sumev;
    switch(isDB) {
    case 1:
	for(i = 0, sumev = 0; i < nr; i++) {
	    ii = i * nr + i;
	    sumev += x[ii] * x[ii];
	}
	break;
    case 0:
	for(i = 0, sumev = 0; i < nr * nc; i++)
	    sumev += x[i] * x[i];
	break;
    }
    return sumev;
}

/* do_getF: get the F value. Input data will be

 */
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Linpack.h> /* QR */
#include <R_ext/Lapack.h>  /* SVD */

SEXP do_getF(SEXP perms, SEXP E)
{
    int i, j, k, ki, nperm = nrows(perms),
	nr = nrows(E), nc = ncols(E);
    /* SEXP ans = PROTECT(allocMatrix(REALSXP, nperm, 3)); */
    SEXP Y = PROTECT(duplicate(E));

    int *iperm = INTEGER(perms);

    /* double *wtake = (double *) R_alloc(nr, sizeof(double)); */

    /* Elements for LINPACK QR decomposition */

    /* permutations to zero base */
    for(i = 0; i < nperm * nr; i++)
	iperm[i]--;

    for (k = 0; k < nperm; k++) {
	for (i = 0; i < nr; i++) {
	    ki = iperm[k + nperm * i];
	    for(j = 0; j < nc; j++) {
		REAL(Y)[i + nr * j] = REAL(E)[ki + nr * j];
	    }
	}
    }
    UNPROTECT(1);
    return Y; /* return to check permutations */
}
