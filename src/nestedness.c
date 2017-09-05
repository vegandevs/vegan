/* C functions for null model simulation.

   These functions are intended to be called via R functions. The main
   vehicle is make.commsim.R that defines null models.

   The static functions are only visible to other C functions in this
   file and cannot be called from R or from C code in other files
   ("translation units"). The void functions (if there are any) can be
   called from R using .C() interface. The SEXP functions can be
   called from R using .Call() interface. Most actual null models are
   static void and are intended to be called via SEXP functions.

   Note on RNG: static void functions are called by other C functions
   (SEXP functions with .Call interface in R), and we should
   GetRNGstate and PutRNGstate only once in that calling function
   instead of getting and putting RNGstate at every function call. The
   non-static functions can be called directly in R code (with .C
   interface) and they must get and put RNGstate within their code.

*/

#include <R.h>
#include <Rmath.h>
#include <R_ext/Utils.h> /* check user interrupts */

/* Utility functions as macros */

/* Random integer 0..imax */

#define IRAND(imax) (int) (((double) (imax + 1)) * unif_rand())

/* 2 different random integers */

#define I2RAND(vec, m) vec[0] = IRAND(m); \
    do {vec[1] = IRAND(m) ;} while(vec[1] == vec[0])

/*
 * Quasiswap or sum-of-squares reducing swap of Miklos & Podani. A quasiswap
 * step takes a random 2x2 submatrix and adds (-1,+1,+1,-1). If the submatrix
 * was (1,0,0,1) it is swapped to (0,1,1,0), but if it was, say, (2,0,0,1) it
 * is swapped to (1,1,1,0) which reduces sums-of-squares. We start with a
 * random matrix with given marginal totals (from R r2dtable) but possibly
 * some values >1. Then we perform quasiswaps on random 2x2 submatrices so
 * long that only 1 and 0 are left.  The function only does the quasiswaps,
 * and it assumes that input matrix 'm' (dimensions 'nr', 'nc') was produced
 * by r2dtable or some other function with given marginal totals, but some
 * values possibly > 1.
 */

/* row & col indices to a vector index */

#define INDX(i, j, nr) (i) + (nr)*(j)

static void quasiswap(int *m, int *nr, int *nc, int *thin)
{
    int i, n, mtot, ss, row[2], col[2], nr1, nc1, a, b, c, d;
    size_t intcheck;

    nr1 = (*nr) - 1;
    nc1 = (*nc) - 1;

    /* Get matrix total 'mtot' and sum-of-squares 'ss' */

    n = (*nr) * (*nc);
    for (i = 0, mtot = 0, ss = 0; i < n; i++) {
	mtot += m[i];
	ss += m[i] * m[i];
    }

    /* Get R RNG in the calling C function */
    /* GetRNGstate(); */

    /* Quasiswap while there are entries > 1 */

    intcheck  = 0; /* check interrupts */
    while (ss > mtot) {
	for (i = 0; i < *thin; i++) {
	    I2RAND(row, nr1);
	    I2RAND(col, nc1);
	    /* a,b,c,d notation for a 2x2 table */
	    a = INDX(row[0], col[0], *nr);
	    b = INDX(row[0], col[1], *nr);
	    c = INDX(row[1], col[0], *nr);
	    d = INDX(row[1], col[1], *nr);
	    if (m[a] > 0 && m[d] > 0 && m[a] + m[d] - m[b] - m[c] >= 2) {
		ss -= 2 * (m[a] + m[d] - m[b] - m[c] - 2);
		m[a]--;
		m[d]--;
		m[b]++;
		m[c]++;
	    } else if (m[b] > 0 && m[c] > 0 &&
		       m[b] + m[c] - m[a] - m[d] >= 2) {
		ss -= 2 * (m[b] + m[c] - m[a] - m[d] - 2);
		m[a]++;
		m[d]++;
		m[b]--;
		m[c]--;
	    }
	}
	/* interrupt? */
	if (intcheck % 10000 == 9999)
	    R_CheckUserInterrupt();
	intcheck++;
    }

    /* Set R RNG in the calling function */
    /* PutRNGstate(); */
}

/* Trial swap: try 'thin' times and swap when you can. This gives zero
 * to many swaps for one call.
 */

static void trialswap(int *m, int *nr, int *nc, int *thin)
{
    int i, a, b, c, d, row[2], col[2];

    /* Get and Set RNG in calling C function */
    /* GetRNGstate(); */

    for (i=0; i < *thin; i++) {
	I2RAND(row, (*nr) - 1);
	I2RAND(col, (*nc) - 1);
	a = INDX(row[0], col[0], *nr);
	b = INDX(row[0], col[1], *nr);
	c = INDX(row[1], col[0], *nr);
	d = INDX(row[1], col[1], *nr);
        /* there are 16 possible matrices, but only two can be
	 * swapped. Find signature of each matrix with bitwise shift
	 * and OR. */
	switch(m[a] | m[b] << 1 | m[c] << 2 | m[d] << 3) {
	case 6: /* 0110 -> 1001 */
	    m[a] = 1;
	    m[b] = 0;
	    m[c] = 0;
	    m[d] = 1;
	    break;
	case 9: /* 1001 -> 0110 */
	    m[a] = 0;
	    m[b] = 1;
	    m[c] = 1;
	    m[d] = 0;
	    break;
	default:
	    break;
	}
    }

    /* PutRNGstate(); */
}

/* Ordinary swap: swap if you can, stop after you swapped, or repeat
 * thin times. The data matrix 'm' must be binary: this is not
 * checked.
 */

static void swap(int *m, int *nr, int *nc, int *thin)
{

    int i, a, b, c, d, row[2], col[2];
    size_t intcheck;

    /* Get and Put RNG in calling C function */
    /* GetRNGstate(); */

    for (i=0, intcheck=0; i < *thin; i++) {
	for(;;) {
	    if (intcheck % 10000 == 9999)
		R_CheckUserInterrupt();
	    intcheck++;
	    I2RAND(row, (*nr) - 1);
	    I2RAND(col, (*nc) - 1);
	    a = INDX(row[0], col[0], *nr);
	    b = INDX(row[0], col[1], *nr);
	    c = INDX(row[1], col[0], *nr);
	    d = INDX(row[1], col[1], *nr);
	    if(m[a] + m[b] + m[c] + m[d] != 2)
		continue;
	    if (m[a] == 1 && m[d] == 1) {
		m[a] = 0;
		m[d] = 0;
		m[b] = 1;
		m[c] = 1;
		break;
	    } 
	    if (m[c] == 1 && m[b] == 1) {
		m[a] = 1;
		m[d] = 1;
		m[b] = 0;
		m[c] = 0;
		break;
	    }
	}
    }
    /* PutRNGstate(); */
}

/* Strona et al. 2014 (NATURE COMMUNICATIONS | 5:4114 |
 * DOI:10.1038/ncomms5114 | www.nature.com/naturecommunications)
 * suggested a boosted sequential binary swap method. Instead of
 * looking for random 2x2 submatrices, they look for 2 rows and
 * collect a list of unique species that occur only in one row, and
 * allocate these randomly to rows preserving counts.
 */

/* uniq is a work vector to hold indices of unique species (occurring
 * only in one of two random rows). uniq must be allocated in the
 * calling function, with safe size 2 * (max. number of species) or
 * with belt and suspenders 2 * (*nc). */

static void curveball(int *m, int *nr, int *nc, int *thin, int *uniq)
{
    int row[2], i, j, jind, ind, nsp1, nsp2, itmp, tmp;

    /* Set RNG in calling C code */
    /* GetRNGstate(); */

    for (i = 0; i < *thin; i++) {
	/* Random sites */
	I2RAND(row, (*nr)-1);
	/* uniq is a vector of unique species for a random pair of
	   rows, It need not be zeroed between thin loops because ind
	   keeps track of used elements. */
	for (j = 0, ind = -1, nsp1 = 0, nsp2 = 0; j < (*nc); j++) {
	    jind = j * (*nr);
	    if (m[row[0] + jind] > 0 && m[row[1] + jind] == 0) {
		uniq[++ind] = j;
		nsp1++;
	    }
	    if (m[row[1] + jind] > 0 && m[row[0] + jind] == 0) {
		uniq[++ind] = j;
		nsp2++;
	    }
	}
	/* uniq contains indices of unique species: shuffle these and
	 * allocate nsp1 first to row[0] and the rest to row[1] */
	if (nsp1 > 0 && nsp2 > 0) { /* something to swap? */
	    for (j = ind; j >= nsp1; j--) {
		tmp = uniq[j];
		itmp = IRAND(j);
		uniq[j] = uniq[itmp];
		uniq[itmp] = tmp;
	    }
	    for (j = 0; j < nsp1; j++) {
		m[INDX(row[0], uniq[j], *nr)] = 1;
		m[INDX(row[1], uniq[j], *nr)] = 0;
	    }
	    for (j = nsp1; j <= ind; j++) {
		m[INDX(row[0], uniq[j], *nr)] = 0;
		m[INDX(row[1], uniq[j], *nr)] = 1;
	    }
	}
    }

    /* PutRNGstate(); */
}

/* 'swapcount' is a C translation of Peter Solymos's R code. It is
 * similar to 'swap', but can swap > 1 values and so works for
 * quantitative (count) data.
 */


/* 'isDiag' is a utility function for 'swapcount' to find the largest
 * value that can be swapped and whether in diagonal or antidiagonal
 * way. The input is a 2x2 submatrix 'sm'.
*/

static int isDiag(int *sm, int *change)
{
    int i, sX;
    int retval;

    /* sX: number of non-zero cells */
    for (i = 0, sX = 0; i < 4; i++)
	    if (sm[i] > 0)
		    sX++;

    /* default values */
    retval = 0;
    *change = 0;
    switch (sX) {
    case 0:
    case 1:
	    /* nothing to swap*/
	    break;
    case 2:
	    /* diagonal and antidiagonal swappable */
	    if (sm[1] > 0 && sm[2] > 0) {
		    retval = (sm[1] < sm[2]) ? sm[1] : sm[2];
		    if (sm[1] != sm[2])
			    *change = 1;
	    }
	    else if (sm[0] > 0 && sm[3] > 0) { 
		    retval = (sm[0] < sm[3]) ? -sm[0] : -sm[3];
		    if (sm[0] != sm[3])
			    *change = 1;
	    } 
	    break;
    case 3:
	    /* always swappable: case depends on the empty corner */
	    if (sm[0] == 0 || sm[3] == 0) {
		    retval = (sm[1] < sm[2]) ? sm[1] : sm[2];
		    if (sm[1] == sm[2])
			    *change = -1;
	    } else {
		    retval = (sm[0] < sm[3]) ? -sm[0] : -sm[3];
		    if (sm[0] == sm[3])
			    *change = -1;
	    }
	    break;
    case 4:
	    /* always swappable: return diagonal case */
	    retval = (sm[1] < sm[2]) ? sm[1] : sm[2];
	    if (sm[1] == sm[2])
		    *change = -2;
	    else
		    *change = -1;
	    break;
    }
    return retval;
}


/* isDiagFill: Largest swappable element and swap policies for
 * fill-neutral swapping
 */

static int isDiagFill(int *sm)
{
    int i, sX;
    int retval;

    /* sX: number of non-zero cells */
    for (i = 0, sX = 0; i < 4; i++)
	    if (sm[i] > 0)
		    sX++;

    retval = 0;
    switch (sX) {
    case 0:
    case 1:
	    /* nothing to swap*/
	    break;
    case 2:
	    /* equal diagonal and antidiagonal fill-neutrally
	     * swappable */
	    if ((sm[0] == sm[3]) && (sm[1] == sm[2])) {
		    if (sm[1] > 0)
			    retval = (sm[1] < sm[2]) ? sm[1] : sm[2];
		    else
			    retval = (sm[0] < sm[3]) ? -sm[0] : -sm[3];
	    }
	    break;
    case 3:
	    /* fill-neutrally swappable if diagonal & antidiagonal
	     * unequal */
	    if ((sm[0] != sm[3]) && (sm[1] != sm[2])) {
		    if (sm[0] == 0 || sm[3] == 0) {
			    retval = (sm[1] < sm[2]) ? sm[1] : sm[2];
		    } else {
			    retval = (sm[0] < sm[3]) ? -sm[0] : -sm[3];
		    }
	    }
	    break;
    case 4:
	    /* never swappable (minelement-1 always swappable) */
	    break;
    }
    return retval;
}

static void swapcount(int *m, int *nr, int *nc, int *thin)
{
    int row[2], col[2], k, ij[4], changed,
	pm[4] = {1, -1, -1, 1} ;
    int sm[4], ev;
    size_t intcheck;

    /* GetRNGstate(); */

    changed = 0;
    intcheck = 0;
    while (changed < *thin) {
	/* Select a random 2x2 matrix*/
	I2RAND(row, *nr - 1);
	I2RAND(col, *nc - 1);
	ij[0] = INDX(row[0], col[0], *nr);
	ij[1] = INDX(row[1], col[0], *nr);
	ij[2] = INDX(row[0], col[1], *nr);
	ij[3] = INDX(row[1], col[1], *nr);
	for (k = 0; k < 4; k ++)
	    sm[k] = m[ij[k]];
	/* The largest value that can be swapped */
	ev = isDiagFill(sm);
 	if (ev != 0) { 
		for (k = 0; k < 4; k++)
			m[ij[k]] += pm[k]*ev;
		changed++;
	}
	if (intcheck % 10000 == 9999)
	    R_CheckUserInterrupt();
	intcheck++;
    }

    /* PutRNGstate(); */
}

/* rswapcount for "reducing swap of count data" is a minor variant of
 * swapcount, but it tries to change the fill: you first make a matrix
 * with correct marginal totals, but possibly wrong fill and then you
 * run this to change the fill while maintaining the totals. The idea
 * is similar as quasiswap for presence/absence data.
 */

static void rswapcount(int *m, int *nr, int *nc, int *mfill)
{
    int row[2], col[2], i, k, ij[4], n, change, cfill,
       pm[4] = {1, -1, -1, 1} ;
    int sm[4], ev;
    size_t intcheck;

    /* Get the current fill 'cfill' */
    n = (*nr) * (*nc);
    for (i = 0, cfill=0; i < n; i++) {
	if (m[i] > 0) 
	    cfill++;
    }
    /* GetRNGstate in calling C code */
    /* GetRNGstate(); */

    /* Loop while fills differ */
    intcheck = 0;
    while (cfill != *mfill) {
	/* Select a random 2x2 matrix*/
	I2RAND(row, *nr - 1);
	I2RAND(col, *nc - 1);
	ij[0] = INDX(row[0], col[0], *nr);
	ij[1] = INDX(row[1], col[0], *nr);
	ij[2] = INDX(row[0], col[1], *nr);
	ij[3] = INDX(row[1], col[1], *nr);
	for (k = 0; k < 4; k ++)
	    sm[k] = m[ij[k]];
	/* The largest value that can be swapped */
	ev = isDiag(sm, &change);
	if (ev != 0) {
	    /* Fill does not change, but swap to bail out from
	     * non-swappable configurations */
	    if (change == 0) {
		for (k = 0; k < 4; k++)
		    m[ij[k]] += pm[k]*ev;
	    } 
	    else if ((change < 0 && *mfill < cfill) ||
		     (change > 0 && *mfill > cfill)) {
		for (k = 0; k < 4; k++)
		    m[ij[k]] += pm[k]*ev;
		cfill += change;
	    } 
	}
	if (intcheck % 10000 == 9999)
	    R_CheckUserInterrupt();
	intcheck++;
    }
    /* PutRNGstate(); */
}

/* 'isDiagSimple' needed for 'abuswap' */

static int isDiagSimple(double *sm)
{
    int i, sX;
    int retval = 0;

    /* sX: number of non-zero cells */
    for (i = 0, sX = 0; i < 4; i++)
	if (sm[i] > 0)
	    sX++;
    
    switch(sX) {
    case 0:
    case 1:
	    /* never swappable */
	    retval = 0;
	    break;
    case 2:
	    /* diagonal and antidiagonal swappable */
	    if ((sm[1] > 0 && sm[2] > 0) || (sm[0] > 0 && sm[3] > 0))
		    retval = 1;
	    else
		    retval = 0;
	    break;
    case 3:
	    /* never swappable */
	    retval = 0;
	    break;
    case 4:
	    /* always swappable */
	    retval = 1;
	    break;
    }
    return retval;
}

/* 'abuswap' to do Hardy 2008 J Ecol 96: 914-926 */

static void abuswap(double *m, int *nr, int *nc, int *thin, int *direct)
{
    int row[2], col[2], k, ij[4], changed, ev;
    size_t intcheck;
    double sm[4];

    /* GetRNGstate in calling C code */
    /* GetRNGstate(); */

    changed = 0;
    intcheck = 0;
    while (changed < *thin) {
	/* Select a random 2x2 matrix*/
	 I2RAND(row, *nr - 1);
	 I2RAND(col, *nc - 1);
	 ij[0] = INDX(row[0], col[0], *nr);
	 ij[1] = INDX(row[1], col[0], *nr);
	 ij[2] = INDX(row[0], col[1], *nr);
	 ij[3] = INDX(row[1], col[1], *nr);
	 for (k = 0; k < 4; k++)
	      sm[k] = m[ij[k]];
	 ev = isDiagSimple(sm);
	 /* Swap */
	 if (ev == 1) {
	      /* fixed column sums */
	      if (*direct == 0) {
		   m[ij[0]] = sm[1];
		   m[ij[1]] = sm[0];
		   m[ij[2]] = sm[3];
		   m[ij[3]] = sm[2];
	      }
	      /* fixed row sums */
	      else {
		   m[ij[0]] = sm[2];
		   m[ij[1]] = sm[3];
		   m[ij[2]] = sm[0];
		   m[ij[3]] = sm[1];
	      }
	      changed++;
	 }
	 if (intcheck % 10000 == 9999)
	     R_CheckUserInterrupt();
	 intcheck++;
    }
    
    /* PutRNGstate(); */
}


/* .Call wrappers to nestedness functions for make.commsim.R */

#include <Rinternals.h>
#include <string.h>

/* Sequential methods:

void trialswap(int *m,    int *nr, int *nc, int *thin)
void      swap(int *m,    int *nr, int *nc, int *thin)
void swapcount(int *m,    int *nr, int *nc, int *thin)
void curveball(int *m,    int *nr, int *nc, int *thin, int *uniq)
void   abuswap(double *m, int *nr, int *nc, int *thin, int *direct)

*/


/*  Sequential swap models: static void functions trialswap, swap and
 *  swapcount have identical signatures and all called via do_swap
 *  with .Call() interface in the R code.
*/

static void (*swap_fun)(int*, int*, int*, int*);

SEXP do_swap(SEXP x, SEXP nsim, SEXP thin, SEXP method)
{
    int nr = nrows(x), nc = ncols(x), ny = asInteger(nsim),
	ithin = asInteger(thin);
    size_t i, N = nr*nc;

    /* trialswap, swap and swapcount have identical function signature */
    const char *cmethod = CHAR(STRING_ELT(method, 0));
    if (strcmp("trialswap", cmethod) == 0)
	swap_fun = trialswap;
    else if (strcmp("swap", cmethod) == 0)
	swap_fun = swap;
    else if (strcmp("swapcount", cmethod) == 0)
	swap_fun = swapcount;
    else
	error("unknown sequential null model \"%s\"", cmethod);

    SEXP out = PROTECT(alloc3DArray(INTSXP, nr, nc, ny));
    int *iout = INTEGER(out);
    if(TYPEOF(x) != INTSXP)
	x = coerceVector(x, INTSXP);
    PROTECT(x);

    int *ix = (int *) R_alloc(N, sizeof(int));

    /* sequential trialswap of ix and save result to the iout
       array */
    memcpy(ix, INTEGER(x), N * sizeof(int));
    GetRNGstate();
    for(i = 0; i < ny; i++) {
	swap_fun(ix, &nr, &nc, &ithin);
	memcpy(iout + i * N, ix, N * sizeof(int));
    }
    PutRNGstate();
    UNPROTECT(2);
    return out;
}

/* curveball has five arguments and needs a work vector, but otherwise
 * the code below mostly duplicate do_swap. Curveball could be
 * combined to do_swap with a couple of if's. */

SEXP do_curveball(SEXP x, SEXP nsim, SEXP thin)
{
    int nr = nrows(x), nc = ncols(x), ny = asInteger(nsim),
	ithin = asInteger(thin);
    size_t i, N = nr*nc;

    SEXP out = PROTECT(alloc3DArray(INTSXP, nr, nc, ny));
    int *iout = INTEGER(out);
    if(TYPEOF(x) != INTSXP)
	x = coerceVector(x, INTSXP);
    PROTECT(x);

    /* difference to do_swap: need a work vector */
    int *iwork = (int *) R_alloc(2*nc, sizeof(int));
    int *ix = (int *) R_alloc(N, sizeof(int));

    /* sequential trialswap of ix and save result to the iout
       array */
    memcpy(ix, INTEGER(x), N * sizeof(int));
    GetRNGstate();
    for(i = 0; i < ny; i++) {
	/* different call than in do_swap */
	curveball(ix, &nr, &nc, &ithin, iwork);
	memcpy(iout + i * N, ix, N * sizeof(int));
    }
    PutRNGstate();
    UNPROTECT(2);
    return out;
}

/* Mostly similar to do_swap, but works on REALSXP instead of INTSXP,
 * and has five arguments (and the last is input, unlike in
 * curveball) */

SEXP do_abuswap(SEXP x, SEXP nsim, SEXP thin, SEXP direct)
{
    int nr = nrows(x), nc = ncols(x), ny = asInteger(nsim),
	ithin = asInteger(thin), idirect = asInteger(direct);
    size_t i, N = nr*nc;

    SEXP out = PROTECT(alloc3DArray(REALSXP, nr, nc, ny));
    double *rout = REAL(out);
    if(TYPEOF(x) != REALSXP)
	x = coerceVector(x, REALSXP);
    PROTECT(x);

    double *rx = (double *) R_alloc(N, sizeof(double));

    /* sequential swap as in do_swap */
    memcpy(rx, REAL(x), N * sizeof(double));
    GetRNGstate();
    for(i = 0; i < ny; i++) {
	abuswap(rx, &nr, &nc, &ithin, &idirect);
	memcpy(rout + i * N, rx, N * sizeof(double));
    }
    PutRNGstate();
    UNPROTECT(2);
    return out;
}

/* Non-sequential methods:

void  quasiswap(int *m, int *nr, int *nc, int *thin)
void rswapcount(int *m, int *nr, int *nc, int *mfill)

*/


/* SEXP x should be 3D array from r2dtable. This will be changed in
   situ, and the original data will be overwritten with quasiswapped
   data. The function does not duplicate its argument, and input x
   will be overwritten.
*/

/* NOTE: Current input 'x' is a 3D array from R function r2dtable
 *  (which is pretty fast). The underlying C function 'rcont2' was
 *  made user callable in R-devel subversion commit 71765 ( hornik |
 *  2016-12-09 17:58:47 +0200 (Fri, 09 Dec 2016) ) and we may consider
 *  calling it directly here in the future.
 */

static void (*qswap_fun)(int *, int *, int *, int *);

SEXP do_qswap(SEXP x, SEXP nsim, SEXP arg4, SEXP method)
{
    /* arg4 is thin for quasiswap and fill for rswapcount */
    int nr = nrows(x), nc = ncols(x), ny = asInteger(nsim),
	iarg4 = asInteger(arg4);
    size_t i, N = nr*nc;

    /* quasiswap and rswapcount have identical function signatures */
    const char *cmethod = CHAR(STRING_ELT(method, 0));
    if (strcmp("quasiswap", cmethod) == 0)
	qswap_fun = quasiswap;
    else if (strcmp("rswapcount", cmethod) == 0)
	qswap_fun = rswapcount;
    else
	error("unknown null model \"%s\"", cmethod);

    /* we must check that input x is integer: some null models set
     * storage.mode "double". */
    if (TYPEOF(x) != INTSXP)
	x = coerceVector(x, INTSXP);
    PROTECT(x);
    int *ix = INTEGER(x);

    GetRNGstate();
    for(i = 0; i < ny; i++) {
	qswap_fun(ix + i * N, &nr, &nc, &iarg4);
    }
    PutRNGstate();
    UNPROTECT(1);
    return x;
}

/* Fill an array of matrices with ones honouring row and column
 * totals. This is similar as Podani's original suggestion for initial
 * filling of a matrix for quasiswap. We have used r2dtable() for
 * quasiswap, but Podani originally studied a method where rows and
 * columns are picked with equal probabilities as long as they still
 * can be filled. */

SEXP do_rcfill(SEXP n, SEXP rs, SEXP cs)
{
    int nrow = length(rs), ncol = length(cs), nmat = asInteger(n);
    int *rfill, *cfill, *rind, *cind;
    int rlen, clen, i, j, k, offset;

    if(TYPEOF(rs) != INTSXP)
	rs = coerceVector(rs, INTSXP);
    PROTECT(rs);
    if(TYPEOF(cs) != INTSXP)
	cs = coerceVector(cs, INTSXP);
    PROTECT(cs);
    int *rowsum = INTEGER(rs);
    int *colsum = INTEGER(cs);
    
    rfill = (int *) R_alloc(nrow, sizeof(int));
    cfill = (int *) R_alloc(ncol, sizeof(int));
    rind = (int *) R_alloc(nrow, sizeof(int));
    cind = (int *) R_alloc(ncol, sizeof(int));
    SEXP out = PROTECT(alloc3DArray(INTSXP, nrow, ncol, nmat));
    int *x = INTEGER(out);
    memset(x, 0, nmat * nrow * ncol * sizeof(int));

    /* collect matrices */
    GetRNGstate();
    for (k = 0; k < nmat; k++) {
	offset = k * nrow * ncol;
	/* initialize fills (0) and indices (0..n) */
	for (i = 0, rlen = -1; i < nrow; i++) {
	    if (rowsum[i] > 0) /* skip empty rows */
		rind[++rlen] = i;
	    rfill[i] = 0;
	}
	for (j = 0, clen = -1; j < ncol; j++) {
	    if (colsum[j] > 0)
		cind[++clen] = j;
	    cfill[j] = 0;
	}
	/* items 0..rlen/clen of rind/cind have the indices of
	 * rows/columns which still can be filled. When the row/column
	 * gets full, replace its index with the last index and reduce
	 * rlen/clen by one */
	while (rlen >= 0) {
	    i = IRAND(rlen);
	    j = IRAND(clen);
	    x[offset + rind[i] + nrow * cind[j]]++;
	    if (++rfill[rind[i]] >= rowsum[rind[i]])
		rind[i] = rind[rlen--];
	    if (++cfill[cind[j]] >= colsum[cind[j]])
		cind[j] = cind[clen--];
	}
    }
    PutRNGstate();
    UNPROTECT(3);
    return out;
}

/* backtracking is a brute force method to fill a matrix with 1's
 * honouring margin totals: do something and if it fails, do something
 * else (Sedgewick). The approach here may not be the fastest, but it
 * is fun to do. We have a vector of indices 'ind' with three
 * compartment: up to index 'ielig' we have indices of eligible zeros
 * that can be picked, then up to 'izero' we have indices of zeros
 * that cannot be picked because their row or column sums are filled,
 * and after 'izero' we have 'npick' indices that we have picked and
 * that will be 1. We fill as long as there are eligible indices, and
 * if there are none but we need to pick more, we "backtrack" or
 * remove picked items, update marginal sums and eligible points. This
 * means of lot of swapping. */

#define EMPTY (-1)
#define SWAP(a,b) tmp=a;a=b;b=tmp

static void backtrack(int *out, int *rowsum, int *colsum, int fill,
		      int nr, int nc, int *rfill, int *cfill, int *ind)
{
    int tmp, i, ir, ic;
    int izero = nr * nc - 1, ielig = nr * nc - 1, npick = 0;

    /* initialize */
    for(i = 0; i < nr * nc; i++)
	ind[i] = i;
    memset(rfill, 0, nr * sizeof(int));
    memset(cfill, 0, nc * sizeof(int));
    
    /* Start working */
    while(npick < fill) { /* outernmost loop (placeholder) */
	/* fill */
	while(ielig > EMPTY) {
	    i = IRAND(ielig); /* eligible: always succeed */
	    ir = ind[i] % nr; /* row */
	    ic = ind[i] / nr; /* column */
	    npick++;
	    SWAP(ind[i], ind[ielig]); /* move after izero */
	    if (ielig < izero) {
		SWAP(ind[ielig], ind[izero]);
	    }
	    ielig--;
	    izero--;
	    /* update fills and move from eligible if marginal sum reached */
	    if (++rfill[ir] == rowsum[ir])
		for(i = ielig; i > EMPTY; i--)
		    if (ind[i] % nr == ir) {
			SWAP(ind[i], ind[ielig]);
			ielig--;
		    }
	    if (++cfill[ic] == colsum[ic])
		for (i = ielig; i > EMPTY; i--)
		    if (ind[i] / nr == ic) {
			SWAP(ind[i], ind[ielig]);
			ielig--;
		    }
	}
	/* get out */
	R_CheckUserInterrupt();
	if (npick == fill)
	    break;
	/* backtrack: remove a picked item and update marginal totals
	 * and see if any items become eligible */
	i = IRAND(npick-1) + izero + 1;
	rfill[ind[i] % nr]--;
	cfill[ind[i] / nr]--;
	npick--;
	SWAP(ind[i], ind[izero + 1]);
	SWAP(ind[izero + 1], ind[0]); /* elig is EMPTY: move to ind[0] */
	ielig = 0;
	izero++;
	/* see what can be moved to eligible */
	for (i = izero; i > ielig; i--) {
	    ir = ind[i] % nr;
	    ic = ind[i] / nr;
	    if (rfill[ir] < rowsum[ir] && cfill[ic] < colsum[ic]) {
		SWAP(ind[i], ind[ielig+1]);
		ielig++;
	    }
	}
    }

    /* output */
    memset(out, 0, nr * nc * sizeof(int));
    /* put 1's in their places */
    for (i = izero+1; i < nr*nc; i++)
	out[ind[i]] = 1;
}

/* .Call interface to backtrack. Input arguments are n (number of
 * matrices), rs (rowsums) and cs (colsums). */

SEXP do_backtrack(SEXP n, SEXP rs, SEXP cs)
{
    int i, fill, nr = length(rs), nc = length(cs), nmat = asInteger(n);
    int N = nr * nc;

    /* check & cast */
    if(TYPEOF(rs) != INTSXP)
	rs = coerceVector(rs, INTSXP);
    PROTECT(rs);
    if(TYPEOF(cs) != INTSXP)
	cs = coerceVector(cs, INTSXP);
    PROTECT(cs);
    int *rowsum = INTEGER(rs);
    int *colsum = INTEGER(cs);

    /* initialize work arrays for backtrack()*/
    int *ind = (int *) R_alloc(nr * nc, sizeof(int));
    int *rfill = (int *) R_alloc(nr, sizeof(int));
    int * cfill = (int *) R_alloc(nc, sizeof(int));
    for (i = 0, fill = 0; i < nr; i++)
	fill += rowsum[i];
    int *x = (int *) R_alloc(nr * nc, sizeof(int));

    SEXP out =  PROTECT(alloc3DArray(INTSXP, nr, nc, nmat));
    int *iout = INTEGER(out);

    GetRNGstate();
    /* Call static C function */
    for(i = 0; i < nmat; i++) {
	backtrack(x, rowsum, colsum, fill, nr, nc, rfill, cfill, ind);
	memcpy(iout + i * N, x, N * sizeof(int));
    }
    PutRNGstate();

    UNPROTECT(3);
    return out;
}
	
#undef EMPTY
#undef SWAP
/* undef: do_backtrack */

#undef IRAND
#undef INDX
/* undef: all of nestedness.c */
