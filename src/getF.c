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

/* do_getF: get the F value. At the first pass *only* for RDA.

 */
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Linpack.h> /* QR */
#include <R_ext/Lapack.h>  /* SVD */

SEXP do_getF(SEXP perms, SEXP E, SEXP QR)
{
    int i, j, k, ki, nperm = nrows(perms),
	nr = nrows(E), nc = ncols(E);
    SEXP ans = PROTECT(allocVector(REALSXP, nperm));
    double *rans = REAL(ans);
    SEXP Y = PROTECT(duplicate(E));
    double *rY = REAL(Y);

    /* pointers and new objects to the QR decomposition */

    double *qr = REAL(VECTOR_ELT(QR, 0));
    int qrank = asInteger(VECTOR_ELT(QR, 1));
    double *qraux = REAL(VECTOR_ELT(QR, 2));
    /* int *pivot = INTEGER(VECTOR_ELT(QR, 3)); */

    double *fitted = (double *) R_alloc(nr * nc, sizeof(double));
    /* int *ny = (int *) R_alloc(nc, sizeof(int)); */
    double dummy;
    int info, qrfit = 1;

    /* double *wtake = (double *) R_alloc(nr, sizeof(double)); */

    /* permutation matrix must be duplicated */
    SEXP dperms = PROTECT(duplicate(perms));
    int *iperm = INTEGER(dperms);

    /* permutations to zero base */
    for(i = 0; i < nperm * nr; i++)
	iperm[i]--;

    /* loop over rows of permutation matrix */
    for (k = 0; k < nperm; k++) {
	/* Y will be permuted data */
	for (i = 0; i < nr; i++) {
	    ki = iperm[k + nperm * i];
	    for(j = 0; j < nc; j++) {
		rY[i + nr * j] = REAL(E)[ki + nr * j];
	    }
	}

	/* Partial models not yet implemented: put them here */

	/* qr.fitted(QR, Y) with LINPACK */
	for (i = 0; i < nc; i++)
	    F77_CALL(dqrsl)(qr, &nr, &nr, &qrank, qraux, rY + i*nr, &dummy,
			    rY + i*nr, &dummy, &dummy, fitted + i*nr, &qrfit,
			    &info);

	/* Eigenvalues: only sum of all, first ev not yet implemented */

	rans[k] = getEV(fitted, nr, nc, 0);

    } /* end permutation loop */

    UNPROTECT(3);
    return ans; /* return to check permutations */
}
