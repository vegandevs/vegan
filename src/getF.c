/* Function to evaluate F-value in permutest.cca. For instance, in R
 * CMD check in Macbook Air this function uses 1/4 of computing time,
 * and in applications with constrained ordination this is the major
 * function. Even small speed-up in this function will have a
 * considerable impact in running time.
 */

/* handle passing strings to Fortran from C */
#define USE_FC_LEN_T       /* from WRExt section 6.6.1 */

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Linpack.h> /* QR */
#include <R_ext/Lapack.h>  /* SVD, eigen */
#include <R_ext/Applic.h> /* R version of QR decomposition dqrdc2: not
			     in public API */

/* handle passing strings to Fortran from C */
#ifndef FCONE              /* from Writing R Extensions section 6.6.1 */
# define FCONE
#endif

#include <math.h> /* sqrt */
#include <string.h> /* memcpy, memset */

/* The following file is in goffactor.c file in vegan */
extern
void wcentre(double *, double *, double *, int *, int *);

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
    int len = nr*nc, info, lwork;
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
		     &nr, &dummy, &nc, &query, &lwork, iwork, &info FCONE);
    if (info != 0)
	error("error %d from Lapack dgesdd", info);
    lwork = (int) query;
    double *work = (double *) R_alloc(lwork, sizeof(double));
    /* call svd */
    F77_CALL(dgesdd)(jobz, &nr, &nc, xwork, &nr, sigma, &dummy,
		     &nr, &dummy, &nc, work, &lwork, iwork, &info FCONE);
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
    int len = nr*nr;

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
		     &tmp, &lwork, &itmp, &liwork, &info FCONE FCONE FCONE);
    if (info != 0)
	error("error %d in work query in LAPACK routine dsyevr", info);
    lwork = (int) tmp;
    liwork = itmp;
    work = (double *) R_alloc(lwork, sizeof(double));
    iwork = (int *) R_alloc(liwork, sizeof(int));

    /* Finally run the eigenanalysis */
    F77_CALL(dsyevr)(jobz, range, uplo, &nr, rx, &nr, &vl, &vu, &il, &iu,
		     &abstol, &naxes, eval, &dummy, &nr, isuppz,
		     work, &lwork, iwork, &liwork, &info FCONE FCONE FCONE);
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

static void qrXw(double *qr, int rank, double *qraux, int *pivot, double *X,
    double *w, int nr, int nc, int discard)
{
    int i, j, ij, len = nr*nc, info = 0, qrkind;
    double dummy[1] = {0.0}, wsqrt;
    double *xwork = (double *) R_alloc(len, sizeof(double));
    /* Extract  R from qr into upper triangle of X */
    for(i = 0; i < len; i++)
	xwork[i] = 0;
    for(j = 0; j < nc; j++)
	for(i = 0; i <= j; i++) {
	    ij = i + nr*j;
	    xwork[ij] = qr[ij];
	}
    /* pivot to zero-base with option to discard first columns */
    for(j = 0; j < nc; j++)
        pivot[j] = pivot[j] - 1 - discard;
    /* Find data as Qy: if y = R then X = QR. The data will over-write
       R. No pivoting, and aliased variables will be moved to last
       columns. Uses Linpack. */
    qrkind = QY;
    /* fill X in the order of the pivot */
    for(j = 0; j < nc; j++) {
	if (pivot[j] >= 0)
	    F77_CALL(dqrsl)(qr, &nr, &nr, &rank, qraux, xwork + j*nr,
			    X + pivot[j]*nr, dummy, dummy, dummy, dummy,
			    &qrkind, &info);
    }

    /* de-weight X */
    for(i = 0; i < nr; i++) {
	wsqrt = sqrt(w[i]);
	for (j = 0; j < nc; j++)
	    X[i + nr*j] /= wsqrt;
    }
}

/* function to test qrX from R. Use with CCA model 'm' as
   .Call("test_qrXw", m$CCA$QR, weights(m)) */

SEXP test_qrXw(SEXP QR, SEXP w, SEXP discard)
{
    int nc, nr;
    QR = PROTECT(duplicate(QR));
    double *qr = REAL(VECTOR_ELT(QR, 0));
    int rank = asInteger(VECTOR_ELT(QR, 1));
    double *qraux = REAL(VECTOR_ELT(QR, 2));
    int *pivot = INTEGER(VECTOR_ELT(QR, 3));
    nr = nrows(VECTOR_ELT(QR, 0));
    nc = ncols(VECTOR_ELT(QR, 0));
    SEXP X = PROTECT(allocMatrix(REALSXP, nr, nc));
    memset(REAL(X), 0, nr * nc * sizeof(double));
    qrXw(qr, rank, qraux, pivot, REAL(X), REAL(w), nr, nc, asInteger(discard));
    UNPROTECT(2);
    return X;
}

/* QR decomposition. Actually R does not use this function, but an
   edited version. From the user point of view, the main difference is
   that the R function returns rank, but this does not do so
   directly. So we test here if we can make this compatible with
   R. This function can be called as .Call("do_QR", x), where x is a
   (centred) matrix. The function returns an object of class "qr" and
   similar to R::qr() result object, expect that the order of columns
   (pivoting) is different than in R.  */

SEXP do_QR(SEXP x)
{
    /* set up */
    int k;
    int nr = nrows(x), nx = ncols(x);
    double TOL = 1e-7;
    SEXP qraux = PROTECT(allocVector(REALSXP, nx));
    SEXP pivot = PROTECT(allocVector(INTSXP, nx));
    memset(INTEGER(pivot), 0, nx * sizeof(int));
    double *work = (double *) R_alloc(nx, sizeof(double));
    int job = 1;
    x = PROTECT(duplicate(x));

    /* QR decomposition with Linpack */
    F77_CALL(dqrdc)(REAL(x), &nr, &nr, &nx, REAL(qraux),
		    INTEGER(pivot), work, &job);
    /* get rank */
    for (k = 1; k < nx; k++) {
	if (fabs(REAL(x)[nr * k + k]) < fabs(TOL * REAL(x)[0]))
	    break;
    }
	
    /* pack up */
    SEXP qr = PROTECT(allocVector(VECSXP, 4));
    SEXP labs = PROTECT(allocVector(STRSXP, 4));
    SET_STRING_ELT(labs, 0, mkChar("qr"));
    SET_STRING_ELT(labs, 1, mkChar("rank")); 
    SET_STRING_ELT(labs, 2, mkChar("qraux"));
    SET_STRING_ELT(labs, 3, mkChar("pivot"));
    setAttrib(qr, R_NamesSymbol, labs);
    SEXP cl = PROTECT(allocVector(STRSXP, 1));
    SET_STRING_ELT(cl, 0, mkChar("qr"));
    classgets(qr, cl);
    UNPROTECT(2); /* cl, labs */
    SET_VECTOR_ELT(qr, 0, x);
    SET_VECTOR_ELT(qr, 1, ScalarInteger(k));
    SET_VECTOR_ELT(qr, 2, qraux);
    SET_VECTOR_ELT(qr, 3, pivot);
    UNPROTECT(4); /* qr, x, pivot, qraux */
    return qr;
}

/* Function do_getF is modelled after R function getF embedded in
 * permutest.cca. The do_getF provides a drop-in replacement to the R
 * function, and is called directly the R function */

SEXP do_getF(SEXP perms, SEXP E, SEXP QR, SEXP QZ,  SEXP effects,
             SEXP w, SEXP first, SEXP isPartial, SEXP isCCA, SEXP isDB)
{
    int i, j, k, ki, p, nterms = length(effects),
	nperm = nrows(perms), nr = nrows(E), nc = ncols(E),
	FIRST = asInteger(first), PARTIAL = asInteger(isPartial),
	DISTBASED = asInteger(isDB), WEIGHTED = asInteger(isCCA);
    /* check that we got terms */
    if (nterms == 0)
        error("model has no terms to test");
    /* check that permutations matrix has correct number of
     * observations */
    if (ncols(perms) != nr)
	error("\'permutations\' matrix should have %d columns, but it has %d",
	      nr, ncols(perms));
    double ev1, ev0, ev;
    SEXP ans = PROTECT(allocMatrix(REALSXP, nperm, nterms + 1));
    double *rans = REAL(ans);
    memset(rans, 0, nperm * (nterms + 1) * sizeof(double));
    SEXP Y = PROTECT(duplicate(E));
    QR = PROTECT(duplicate(QR));
    QZ = PROTECT(duplicate(QZ));
    double *rY = REAL(Y);
    if (TYPEOF(effects) != INTSXP)
	effects = coerceVector(effects, INTSXP);
    PROTECT(effects);
    int *term = INTEGER(effects);
    
    /* pointers and new objects to the QR decomposition */
    
    double *qr = REAL(VECTOR_ELT(QR, 0));
    int qrank = asInteger(VECTOR_ELT(QR, 1));
    double *qraux = REAL(VECTOR_ELT(QR, 2));
    int *pivot = INTEGER(VECTOR_ELT(QR, 3));
    double *Zqr, *Zqraux;
    int Zqrank;
    if (PARTIAL) {
	Zqr = REAL(VECTOR_ELT(QZ, 0)); 
	Zqrank = asInteger(VECTOR_ELT(QZ, 1));
	Zqraux = REAL(VECTOR_ELT(QZ, 2));
    }
    
    double *fitted = (double *) R_alloc(nr * nc, sizeof(double));
    double *resid = (double *) R_alloc(nr * nc, sizeof(double));
    /* work array and variables for QR decomposition */
    double *qty = (double *) R_alloc(nr, sizeof(double));
    double dummy[1] = {0.0};
    int info, qrkind;
    /* Weighted methods currently need re-evaluation of QR
       decomposition (probably changed in the future, but now for the
       compatibility with the current code). For this we need to
       reconstruct constraints and conditions. */
    int nx = ncols(VECTOR_ELT(QR, 0)), nz = 0, *zpivot;
    double *wperm, *Xorig, *Zorig, *Zperm, *qrwork, *zqrwork, qrtol=1e-7;
    if (WEIGHTED) {
	if (PARTIAL) {
	    nz = ncols(VECTOR_ELT(QZ, 0));
	    /* want to have weighted Z: set all weights = 1 */
	    double *w1 = (double *) R_alloc(nr, sizeof(double));
	    for(i = 0; i < nr; i++)
		w1[i] = 1;
	    zpivot = INTEGER(VECTOR_ELT(QZ, 3));
	    Zorig = (double *) R_alloc(nr * nz, sizeof(double));
	    memset(Zorig, 0, nr * nz * sizeof(double));
	    Zperm = (double *) R_alloc(nr * nz, sizeof(double));
	    qrXw(Zqr, Zqrank, Zqraux, zpivot, Zorig, w1, nr, nz, 0);
	    zqrwork = (double *) R_alloc(2 * nz, sizeof(double));
	}
	wperm = (double *) R_alloc(nr, sizeof(double));
	Xorig = (double *) R_alloc(nr * nx, sizeof(double));
	memset(Xorig, 0, nr * nx * sizeof(double));
	qrXw(qr, qrank, qraux, pivot, Xorig, REAL(w), nr, nx, nz);
	qrwork = (double *) R_alloc(2 * nx, sizeof(double));
    }

    /* distance-based methods need to transpose data */
    double *transY;
    if (DISTBASED)
	transY = (double *) R_alloc(nr * nr, sizeof(double));

    /* permutation matrix must be duplicated */
    if (TYPEOF(perms) != INTSXP)
	perms = coerceVector(perms, INTSXP);
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
	    if (PARTIAL)
		for(j = 0; j < nz; j++)
		    Zperm[i + nr*j] = Zorig[ki + nr*j];
	    if (WEIGHTED) {
		wperm[i] = REAL(w)[ki];
	    }
	}

	/* Partial model: qr.resid(QZ, Y) with LINPACK */
	if (PARTIAL) {
	    /* Re-do QR decomposition with changed weights */
	    if (WEIGHTED) {
	        memcpy(Zqr, Zperm, nr * nz * sizeof(double));
		/* dqrdc2 is not in R API */
		for (i = 0; i < nz; i++)
		    zpivot[i] = i + 1;
		F77_CALL(dqrdc2)(Zqr, &nr, &nr, &nz, &qrtol, &Zqrank,
	                         Zqraux, zpivot, zqrwork);
	    }
	    qrkind = RESID;
	    for(i = 0; i < nc; i++)
		F77_CALL(dqrsl)(Zqr, &nr, &nr, &Zqrank, Zqraux, rY + i*nr,
	                        dummy, qty, dummy, rY + i*nr, dummy,
				&qrkind, &info);
	    /* distances need symmetric residuals */
	    if (DISTBASED) {
		transpose(rY, transY, nr, nr);
		qrkind = RESID;
		for(i = 0; i < nc; i++)
		    F77_CALL(dqrsl)(Zqr, &nr, &nr, &Zqrank, Zqraux,
				    transY + i*nr, dummy, qty, dummy,
				    rY + i*nr, dummy, &qrkind, &info);
	    }

	}

	/* CONSTRAINED COMPONENT */

	/* Re-weight & residualize constraints and re-do QR */
	if (WEIGHTED) {
	    memcpy(qr, Xorig, nr * nx * sizeof(double));
	    wcentre(Xorig, qr, wperm, &nr, &nx);
	    if (PARTIAL) {
		qrkind = RESID;
		for (i = 0; i < nx; i++)
	            F77_CALL(dqrsl)(Zqr, &nr, &nr, &Zqrank, Zqraux,
	                qr + i*nr, dummy, qty, dummy, qr + i*nr,
	                dummy, &qrkind, &info);
	    }
	    for(i = 0; i < nx; i++)
		pivot[i] = i + 1;
	    F77_CALL(dqrdc2)(qr, &nr, &nr, &nx, &qrtol, &qrank, 
			     qraux, pivot, qrwork);  
	}
	
	/* qr.fitted(QR, Y) + qr.resid(QR, Y) with LINPACK */

	/* If there are effects, we go for all but the full rank first */
	ev0 = 0; /* must be set for later use outside the loop */
	if (nterms > 1) {
	    qrkind = FIT;
	    for (p = 0; p < (nterms - 1); p++) {
		for (i = 0; i < nc; i++)
		    F77_CALL(dqrsl)(qr, &nr, &nr, term + p, qraux, rY + i*nr,
				    dummy, qty, dummy, dummy, fitted + i*nr,
				    &qrkind, &info);
		ev = getEV(fitted, nr, nc, DISTBASED);
		rans[k + p*nperm] = ev - ev0;
		ev0 = ev;
	    }
	}
	/* Evaluate full-rank model */
	if (PARTIAL || FIRST)
	    qrkind = FIT + RESID;
	else
	    qrkind = FIT;
	for (i = 0; i < nc; i++)
	    F77_CALL(dqrsl)(qr, &nr, &nr, &qrank, qraux, rY + i*nr, dummy,
			    qty, dummy, resid + i*nr, fitted + i*nr,
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
				    transY + i*nr, dummy, qty, dummy,
				    dummy, fitted + i*nr, &qrkind, &info);
		ev1 = eigenfirst(fitted, nr);
	    } else {
		ev1 = svdfirst(fitted, nr, nc);
		ev1 = ev1 * ev1;
	    }
	    rans[k] = ev1;
	} else {
	    rans[k + (nterms - 1) * nperm] =
		getEV(fitted, nr, nc, DISTBASED) - ev0;
	}
	if (PARTIAL || FIRST)
	    rans[k + nterms * nperm] = getEV(resid, nr, nc, DISTBASED);

    } /* end permutation loop */

    UNPROTECT(6);
    return ans;
}
