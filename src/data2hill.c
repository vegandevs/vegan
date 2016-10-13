/*********************************************************
 Changes a data matrix into a form used internally in Mark
 Hill's programs DECORANA and TWINSPAN.
 There was an R version, but profiling showed it was a bottleneck.
 Using this C version saved 30& running time in all matrix sizes
 tested (although it matters only for large data sets).
 The original R function is at the end of this function as a comment.
**********************************************************/

#include <R.h>

void data2hill(double *x,
	       int *mi, int *n, int *nid, int *ibegin, int *iend,
	       int *idat, double *qidat)
{
     int nr, nc, i, j, ij, now;
     
     nr = *mi;
     nc = *n;
     if (nr <= 0 || nc <= 0)
	  error("zero extent dimensions");
     
     now=0;
     for (i=0; i<nr; i++) {
	  for (j=0; j<nc; j++) {
	       ij = i+nr*j;
	       if (x[ij] > 0.0) {
		    idat[now] = j+1;
		    qidat[now] = x[ij];
		    now++;
	       }
	  }
	  iend[i] = now;
     }
     ibegin[0] = 1;
     for (i=1; i<nr; i++)
	  ibegin[i] = iend[i-1] + 1;
     *mi = nr;
     *n = nc;
     *nid = now; 
}

/*  "data2hill" <- */
/*    function(x) */
/*  { */
/*    x <- as.matrix(x) */
/*    nc <- ncol(x) */
/*    nr <- nrow(x) */
/*    rcount <- apply(x>0, 1, sum) */
/*    iend <- cumsum(rcount) */
/*    ibegin <- c(0,iend[-nr])+1 */
/*    nz <- iend[nr] */
/*    idat <- rep(NA, nz) */
/*    qidat <- rep(NA, nz) */
/*    sp.index <- 1:nc */
/*    for (j in 1:nr) { */
/*      hit <- x[j,] > 0 */
/*      idat[ibegin[j]:iend[j]] <- sp.index[hit] */
/*      qidat[ibegin[j]:iend[j]] <- x[j, hit] */
/*    } */
/*    list(mi=nr, n=nc, nid=nz, ibegin=ibegin, iend=iend, idat=idat, qidat=qidat) */
/*  } */

/* C interface to hide the ugly calls to decorana.f */

#include <R.h>
#include <Rinternals.h>

/* Fortran routines called from decorana.f */

void F77_NAME(eigy)(double*, double*, double*, int*, int*, int*, double*,
		    int*, int*, int*, int*, int*, int*, int*, double*,
		    double*, double*, double*, double*, double*, double*,
		    double*, int*, int*, int*, double*, double*);
void F77_NAME(cutup)(double*, int*, int*, int*);
void F77_NAME(yxmult)(double*, double*, int*, int*, int*, int*, int*,
		      int*, double*);

SEXP do_decorana(SEXP veg, SEXP ira, SEXP iresc, SEXP rshort, SEXP imk,
		 SEXP aidot, SEXP adotj)
{
    /* decorana CONSTANTS */
    int NAXES = 4, MKPAD = 4;
    double ZEROEIG = 1e-7;
    /* input parameters */
    int ra = asInteger(ira), resc = asInteger(iresc),
	mk = asInteger(imk) + MKPAD;
    double xshort = asReal(rshort);
    /* internal parameters */
    int nr = nrows(veg), nc = ncols(veg), nid;
    int i, j;

    /* PART 1: R data matrix to CEP condensed format */

    /* check type of veg */
    if (TYPEOF(veg) != REALSXP)
	veg = coerceVector(veg, REALSXP);
    PROTECT(veg);
    double *xveg = REAL(veg);
    /* No. of non-zero items in veg */
    for (i = 0, nid = 0; i < nr*nc; i++)
	if (xveg[i] > 0)
	    nid++;
    /* allocate vectors for the CEP format */
    int *ibegin = (int *) R_alloc(nr, sizeof(int));
    int *iend = (int *) R_alloc(nr, sizeof(int));
    int *idat = (int *) R_alloc(nid, sizeof(int));
    double *qidat = (double *) R_alloc(nid, sizeof(double));
    /* data to internal CEP format */
    data2hill(xveg, &nr, &nc, &nid, ibegin, iend, idat, qidat);
    UNPROTECT(1); /* veg */

    /* PART 2: Call decorana Fortran functions */

    /* return objects */
    SEXP xeig = PROTECT(allocMatrix(REALSXP, nr, NAXES));
    SEXP yeig = PROTECT(allocMatrix(REALSXP, nc, NAXES));
    SEXP eig = PROTECT(allocVector(REALSXP, NAXES));
    double *rxeig = REAL(xeig);
    double *ryeig = REAL(yeig);
    double *reig = REAL(eig);

    /* internal vectors for decorana */
    int *ix = (int *) R_alloc(3 * nr, sizeof(int));
    double *ywork = (double *) R_alloc(4 * nc, sizeof(double));

    /* Call decorana.f functions */
    for (i = 0; i < NAXES; i++) {
	F77_CALL(eigy)(rxeig + i*nr, ryeig + i*nc, reig + i, &i, &ra, &resc,
		       &xshort, &nr, &mk, &nc, &nid, ibegin, iend, idat, qidat,
		       ywork, ywork + nc, ywork + 2*nc, ywork + 3*nc,
		       rxeig, rxeig + nr, rxeig + 2*nr,
		       ix, ix + nr, ix + 2*nr,
		       REAL(aidot), REAL(adotj));
	if (reig[i] < ZEROEIG) {
	    for(j = 0; j < nr; j++)
		rxeig[i*nr + j] = 0;
	    for(j = 0; j < nc; j++)
		ryeig[i*nc + j] = 0;
	    reig[i] = 0;
	}
	if (!ra && i != NAXES - 1)
	    F77_CALL(cutup)(rxeig + i*nr, ix + i*nr, &nr, &mk);
    }
    for (i = 0; i < NAXES; i++) {
	F77_CALL(yxmult)(ryeig + i*nc, rxeig + i*nr, &nr, &nc, &nid,
			 ibegin, iend, idat, qidat);
	for (j = 0; j < nr; j++)
	    rxeig[i*nr + j] /= REAL(aidot)[j];
    }

    /* decorana object */
    SEXP result = PROTECT(allocVector(VECSXP, 9));
    SEXP names = PROTECT(allocVector(STRSXP, 9));
    SEXP class = PROTECT(allocVector(STRSXP, 1));
    SET_STRING_ELT(names, 0, mkChar("evals"));
    SET_STRING_ELT(names, 1, mkChar("rproj"));
    SET_STRING_ELT(names, 2, mkChar("cproj"));
    SET_STRING_ELT(names, 3, mkChar("adotj"));
    SET_STRING_ELT(names, 4, mkChar("aidot"));
    SET_STRING_ELT(names, 5, mkChar("ira"));
    SET_STRING_ELT(names, 6, mkChar("iresc"));
    SET_STRING_ELT(names, 7, mkChar("short"));
    SET_STRING_ELT(names, 8, mkChar("mk"));
    setAttrib(result, R_NamesSymbol, names);
    SET_STRING_ELT(class, 0, mkChar("decorana"));
    classgets(result, class);
    UNPROTECT(2); /* class, names */
    SET_VECTOR_ELT(result, 0, eig);
    SET_VECTOR_ELT(result, 1, xeig);
    SET_VECTOR_ELT(result, 2, yeig);
    SET_VECTOR_ELT(result, 3, adotj);
    SET_VECTOR_ELT(result, 4, aidot);
    SET_VECTOR_ELT(result, 5, ira);
    SET_VECTOR_ELT(result, 6, iresc);
    SET_VECTOR_ELT(result, 7, rshort);
    SET_VECTOR_ELT(result, 8, imk);

    UNPROTECT(4); /* xeig, yeig, eig, result */
    return result;
}
