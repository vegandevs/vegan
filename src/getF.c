/* Function to evaluate F-value in permutest.cca. For instance, in R
 * CMD check in Macbook Air this function uses 1/4 of computing time,
 * and in applications with constrained ordination this is the major
 * function. Even small speed-up in this function will have a
 * considerable impact in running time.
 */

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Linpack.h> /* QR */
#include <R_ext/Lapack.h>  /* SVD, eigen */

#include <math.h> /* sqrt */
#include <string.h> /* memcpy */

/* LINPACK uses the same function (dqrsl) to find derived results from
 * the QR decomposition. It uses decimal coding to define the kind of
 * output with following alternatives (although we will not use all of
 * these): */

#define FIT 1
#define RESID 10
#define COEF 100
#define QTY 1000
#define QY 10000

/* Function called form getF to evaluate the sum of all eigenvalues */

static double getEV(double *x, int nr, int nc, int isDB)
{
    int i;
    double sumev;
    if (isDB)
	for(i = 0, sumev = 0; i < nr; i++)
	    sumev += x[i * nr + i];
    else
	for(i = 0, sumev = 0; i < nr * nc; i++)
	    sumev += x[i] * x[i];

    return sumev;
}

/* LAPACK function dgesdd for SVD. Returns first singular value */

static double svdfirst(double *x, int nr, int nc)
{
    char jobz[2] = "N";
    int minrc = (nr < nc) ? nr : nc;
    int i, len = nr*nc, info, lwork;
    double dummy = 0, query;

    /* copy data: dgesdd will destroy the original */
    double *xwork = (double *) R_alloc(len, sizeof(double));
    memcpy(xwork, x, len * sizeof(double));

    /* singular values */
    double *sigma = (double *) R_alloc(minrc, sizeof(double));

    int *iwork = (int *) R_alloc(8 * minrc, sizeof(int));
    /* query and set optimal work array */
    info = 0;
    lwork = -1;
    F77_CALL(dgesdd)(jobz, &nr, &nc, xwork, &nr, sigma, &dummy,
		     &nr, &dummy, &nc, &query, &lwork, iwork, &info);
    if (info != 0)
	error("error %d from Lapack dgesdd", info);
    lwork = (int) query;
    double *work = (double *) R_alloc(lwork, sizeof(double));
    /* call svd */
    F77_CALL(dgesdd)(jobz, &nr, &nc, xwork, &nr, sigma, &dummy,
		     &nr, &dummy, &nc, work, &lwork, iwork, &info);
    if (info != 0)
	error("error %d from Lapack dgesdd, pos 2", info);
    return sigma[0];
}

/* LAPACK function to return the first eigenvalue of symmetric
   matrix */

static double eigenfirst(double *x, int nr)
{
    /* no eigenvectors (jobz), range of eigenvalus (range), lower
       triagle (uplo) */
    char jobz[2] = "N", range[2] = "I", uplo[2] = "L";
    double vl = 0.0, vu = 0.0; /* range of ev magnitudes, not used */
    double abstol = 0.0, dummy = 0.0;
    /* il, iu: we want only largest eigenvalue, and they come in
       *ascending* order */
    int il = nr, iu = nr, naxes = 1;
    double *eval = (double *) R_alloc(nr, sizeof(double));
    int i, len = nr*nr;

    /* work arrays, their sizes and info. */
    int *isuppz = (int *) R_alloc(2 * nr, sizeof(int));
    double *work, tmp;
    int lwork, *iwork, liwork, info, itmp;

    /* input will be destroyed: copy here */
    double *rx = (double *) R_alloc(len, sizeof(double));
    memcpy(rx, x, len * sizeof(double));

    /* query and set optimal work arrays */
    lwork = -1;
    liwork = -1;
    F77_CALL(dsyevr)(jobz, range, uplo, &nr, rx, &nr, &vl, &vu, &il, &iu,
		     &abstol, &naxes, eval, &dummy, &nr, isuppz,
		     &tmp, &lwork, &itmp, &liwork, &info);
    if (info != 0)
	error("error %d in work query in LAPACK routine dsyevr", info);
    lwork = (int) tmp;
    liwork = itmp;
    work = (double *) R_alloc(lwork, sizeof(double));
    iwork = (int *) R_alloc(liwork, sizeof(int));

    /* Finally run the eigenanalysis */
    F77_CALL(dsyevr)(jobz, range, uplo, &nr, rx, &nr, &vl, &vu, &il, &iu,
		     &abstol, &naxes, eval, &dummy, &nr, isuppz,
		     work, &lwork, iwork, &liwork, &info);
    if (info != 0)
	error("error %d in LAPACK routine dsyever", info);
    return eval[0];
}

/* function to test svdfirst & eigenfirst from R */
SEXP test_ev(SEXP x, SEXP svd)
{
    int KIND = asInteger(svd);
    int nr = nrows(x), nc = ncols(x);
    SEXP ans = PROTECT(allocVector(REALSXP, 1));
    if (KIND)
	REAL(ans)[0] = svdfirst(REAL(x), nr, nc);
    else
	REAL(ans)[0] = eigenfirst(REAL(x), nr);
    UNPROTECT(1);
    return ans;
}

/* transpose a matrix. */

static void transpose(double *x, double *tx, int nr, int nc)
{
    int ij, i, j;
    for (i = 0, ij = 0; i < nr; i++)
	for (j = 0; j < nc; j++)
	    tx[ij++] = x[j * nr + i];
}

/* test transpose from an R session */

SEXP test_trans(SEXP x)
{
    int nr = nrows(x), nc = ncols(x);
    SEXP tx = PROTECT(allocMatrix(REALSXP, nc, nr));
    transpose(REAL(x), REAL(tx), nr, nc);
    UNPROTECT(1);
    return tx;
}

/* Reconstruct data X from weighted QR decomposition. We need this
   only for CCA, and there the data are weighted by row sums, but we
   need the original unweighted data. So we do here both qrX and
   de-weighting. */

static void qrXw(double *qr, int rank, double *qraux, double *X, double *w,
		 int nr, int nc)
{
    int i, j, ij, len = nr*nc, info = 0, qrkind;
    double dummy = 0, wsqrt;
    /* Extract  R from qr into upper triangle of X */
    for(i = 0; i < len; i++)
	X[i] = 0;
    for(j = 0; j < nc; j++)
	for(i = 0; i <= j; i++) {
	    ij = i + nr*j;
	    X[ij] = qr[ij];
	}
    /* Find data as Qy: if y = R then X = QR. The data will over-write
       R. No pivoting, and aliased variables will be moved to last
       columns. Uses Linpack. */
    qrkind = QY;
    for(j = 0; j < nc; j++)
	F77_CALL(dqrsl)(qr, &nr, &nr, &rank, qraux, X + j*nr, X + j*nr,
			&dummy, &dummy, &dummy, &dummy, &qrkind, &info);

    /* de-weight X */
    for(i = 0; i < nr; i++) {
	wsqrt = sqrt(w[i]);
	for (j = 0; j < nc; j++)
	    X[i + nr*j] /= wsqrt;
    }
}

/* function to test qrX from R. Use with CCA model 'm' as
   .Call("test_qrXw", m$CCA$QR, weights(m)) */

SEXP test_qrXw(SEXP QR, SEXP w)
{
    int nc, nr;
    double *qr = REAL(VECTOR_ELT(QR, 0));
    int rank = asInteger(VECTOR_ELT(QR, 1));
    double *qraux = REAL(VECTOR_ELT(QR, 2));
    nr = nrows(VECTOR_ELT(QR, 0));
    nc = ncols(VECTOR_ELT(QR, 0));
    SEXP X = PROTECT(allocMatrix(REALSXP, nr, nc));
    qrXw(qr, rank, qraux, REAL(X), REAL(w), nr, nc);
    UNPROTECT(1);
    return X;
}

/* Function do_getF is modelled after R function getF embedded in
 * permutest.cca. The do_getF provides a drop-in replacement to the R
 * function, and is called directly the R function */

SEXP do_getF(SEXP perms, SEXP E, SEXP QR, SEXP QZ, SEXP first,
	     SEXP isPartial, SEXP isDB)
{
    int i, j, k, ki,
	nperm = nrows(perms), nr = nrows(E), nc = ncols(E),
	FIRST = asInteger(first), PARTIAL = asInteger(isPartial),
	DISTBASED = asInteger(isDB);
    double ev1;
    SEXP ans = PROTECT(allocMatrix(REALSXP, nperm, 2));
    double *rans = REAL(ans);
    SEXP Y = PROTECT(duplicate(E));
    double *rY = REAL(Y);

    /* pointers and new objects to the QR decomposition */

    double *qr = REAL(VECTOR_ELT(QR, 0));
    int qrank = asInteger(VECTOR_ELT(QR, 1));
    double *qraux = REAL(VECTOR_ELT(QR, 2));
    double *Zqr, *Zqraux;
    int Zqrank;
    if (PARTIAL) {
	Zqr = REAL(VECTOR_ELT(QZ, 0));
	Zqrank = asInteger(VECTOR_ELT(QZ, 1));
	Zqraux = REAL(VECTOR_ELT(QZ, 2));
    }

    double *fitted = (double *) R_alloc(nr * nc, sizeof(double));
    /* separate resid needed only in some cases */
    double *resid;
    if (PARTIAL || FIRST)
	resid = (double *) R_alloc(nr * nc, sizeof(double));
    /* work array and variables for QR decomposition */
    double *qty = (double *) R_alloc(nr, sizeof(double));
    double dummy;
    int info, qrkind;

    /* distance-based methods need to transpose data */
    double *transY;
    if (DISTBASED)
	transY = (double *) R_alloc(nr * nr, sizeof(double));

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
		if (DISTBASED)    /* shuffle rows & cols symmetrically */
		    rY[i + nr*j] = REAL(E)[ki + nr * iperm[k + nperm*j]];
		else   /* shuffle rows */
		    rY[i + nr*j] = REAL(E)[ki + nr*j];
	    }
	}

	/* Partial model: qr.resid(QZ, Y) with LINPACK */
	if (PARTIAL) {
	    qrkind = RESID;
	    for(i = 0; i < nc; i++)
		F77_CALL(dqrsl)(Zqr, &nr, &nr, &Zqrank, Zqraux, rY + i*nr,
				&dummy, qty, &dummy, rY + i*nr, &dummy,
				&qrkind, &info);
	    /* distances need symmetric residuals */
	    if (DISTBASED) {
		transpose(rY, transY, nr, nr);
		qrkind = RESID;
		for(i = 0; i < nc; i++)
		    F77_CALL(dqrsl)(Zqr, &nr, &nr, &Zqrank, Zqraux,
				    transY + i*nr, &dummy, qty, &dummy,
				    rY + i*nr, &dummy, &qrkind, &info);
	    }
	}

	/* qr.fitted(QR, Y) + qr.resid(QR, Y) with LINPACK */
	if (PARTIAL || FIRST)
	    qrkind = FIT + RESID;
	else
	    qrkind = FIT;
	for (i = 0; i < nc; i++)
	    F77_CALL(dqrsl)(qr, &nr, &nr, &qrank, qraux, rY + i*nr, &dummy,
			    qty, &dummy, resid + i*nr, fitted + i*nr,
			    &qrkind, &info);

	/* Eigenvalues: either sum of all or the first If the sum of
	 * all eigenvalues does not change, we have only ev of CCA
	 * component in the first column, and the second column is
	 * rubbish that should be filled in the calling R function
	 * with the correct value. */

	if (FIRST) {
	    if (DISTBASED) { /* needs symmetric matrix */
		transpose(fitted, transY, nr, nr);
		qrkind = FIT;
		for(i = 0; i < nc; i++)
		    F77_CALL(dqrsl)(qr, &nr, &nr, &qrank, qraux,
				    transY + i*nr, &dummy, qty, &dummy,
				    &dummy, fitted + i*nr, &qrkind, &info);
		ev1 = eigenfirst(fitted, nr);
	    } else {
		ev1 = svdfirst(fitted, nr, nc);
		ev1 = ev1 * ev1;
	    }
	    rans[k] = ev1;
	} else {
	    rans[k] = getEV(fitted, nr, nc, DISTBASED);
	}
	if (PARTIAL || FIRST)
	    rans[k + nperm] = getEV(resid, nr, nc, DISTBASED);

    } /* end permutation loop */

    UNPROTECT(3);
    return ans;
}
