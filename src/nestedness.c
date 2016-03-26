#include <R.h>
#include <Rmath.h>

/* Utility functions */

/* Random integer 0..imax */

#define IRAND(imax) (int) (((double) (imax + 1)) * unif_rand())

/* 2 different random integers */

void i2rand(int *vec, int imax)
{
    vec[0] = IRAND(imax);
    do {
	vec[1] = IRAND(imax);
    } while (vec[1] == vec[0]);
}


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

void quasiswap(int *m, int *nr, int *nc)
{
    int i, n, mtot, ss, row[2], col[2], nr1, nc1, a, b, c, d;

    nr1 = (*nr) - 1;
    nc1 = (*nc) - 1;

    /* Get matrix total 'mtot' and sum-of-squares 'ss' */

    n = (*nr) * (*nc);
    for (i = 0, mtot = 0, ss = 0; i < n; i++) {
	mtot += m[i];
	ss += m[i] * m[i];
    }

    /* Get R RNG */
    GetRNGstate();

    /* Quasiswap while there are entries > 1 */

    while (ss > mtot) {
	i2rand(row, nr1);
	i2rand(col, nc1);
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

    /* Set R RNG */
    PutRNGstate();
}

/* Trial swap: try 'thin' times and swap when you can. This gives zero
 * to many swaps for one call.
 */

void trialswap(int *m, int *nr, int *nc, int *thin)
{

    int i, a, b, c, d, row[2], col[2], sX;

    GetRNGstate();

    for (i=0; i < *thin; i++) {
	i2rand(row, (*nr) - 1);
	i2rand(col, (*nc) - 1);
	a = INDX(row[0], col[0], *nr);
	b = INDX(row[0], col[1], *nr);
	c = INDX(row[1], col[0], *nr);
	d = INDX(row[1], col[1], *nr);
        /* only two filled items can be swapped */
	sX = m[a] + m[b] + m[c] + m[d];
	if (sX != 2)
	    continue;
	if (m[a] == 1 && m[d] == 1) {
	    m[a] = 0;
	    m[d] = 0;
	    m[b] = 1;
	    m[c] = 1;
	} else if (m[c] == 1 && m[b] == 1) {
	    m[a] = 1;
	    m[d] = 1;
	    m[b] = 0;
	    m[c] = 0;
	}
    }

    PutRNGstate();
}

/* Ordinary swap: swap if you can, stop after you swapped, or repeat
 * thin times. The data matrix 'm' must be binary: this is not
 * checked.
 */

void swap(int *m, int *nr, int *nc, int *thin)
{

    int i, a, b, c, d, row[2], col[2], sX;

    GetRNGstate();

    for (i=0; i < *thin; i++) {
	for(;;) {
	    i2rand(row, (*nr) - 1);
	    i2rand(col, (*nc) - 1);
	    a = INDX(row[0], col[0], *nr);
	    b = INDX(row[0], col[1], *nr);
	    c = INDX(row[1], col[0], *nr);
	    d = INDX(row[1], col[1], *nr);
	    sX = m[a] + m[b] + m[c] + m[d];
	    if (sX != 2)
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
    PutRNGstate();
}

/* Strona et al. 2014 (NATURE COMMUNICATIONS | 5:4114 |
 * DOI:10.1038/ncomms5114 | www.nature.com/naturecommunications)
 * suggested a boosted sequential binary swap method. Instead of
 * looking for random 2x2 submatrices, they look for 2 rows and swap
 * every swappable item in one sweep. Swappable item is a presense in
 * only one of the compared rows, and the maximum number of swappable
 * items is the smaller of unique presences in compared sites. They
 * used a bizarre name ("curveball") and narrative for their method,
 * and we also use this name here.
 */

void curveball(int *m, int *nr, int *nc, int *thin, int *wrk1, int *wrk2)
{
    int row[2], i, j, jind, ind1, ind2, itmp, tmp;

    /* Set RNG */
    GetRNGstate();

    for (i = 0; i < *thin; i++) {
	/* Random sites */
	i2rand(row, (*nr)-1);
	/* Vectors of unique species for each random species */
	for (j = 0, ind1 = -1, ind2 = -1; j < (*nc); j++) {
	    jind <- j * (*nr)
	    if (m[row[0] + jind] > 0 && m[row[1] + jind] == 0)
		wrk1[++ind1] = j;
	    if (m[row[1] + jind] > 0 && m[row[0] + jind] == 0)
		wrk2[++ind2] = j;
	}
	/* If both rows had unique species, we swap the smaller number */
	if (ind1 > -1 && ind2 > -1) {
	    /* All items of smaller set are swapped, but we need a random
	     * subset for the longer one */
	    if (ind1 > ind2)
		for (j = ind1; j >= ind2; j--) {
		    tmp = wrk1[j];
		    itmp = IRAND(j);
		    wrk1[itmp] = wrk1[j];
		    wrk1[j] = tmp;
		}
	    if (ind2 > ind1)
		for (j = ind2; j >= ind1; j--)  {
		    tmp = wrk2[j];
		    itmp = IRAND(j);
		    wrk2[itmp] = wrk2[j];
		    wrk2[j] = tmp;
		}
	}
	/* Do swaps */
	if (ind1 > ind2)
	    ind1 = ind2;
	for (j = 0; j < ind1; j++) {
	    m[INDX(row[0], wrk1[j], *nr)] = 0;
	    m[INDX(row[0], wrk2[j], *nr)] = 1;
	    m[INDX(row[1], wrk2[j], *nr)] = 0;
	    m[INDX(row[1], wrk1[j], *nr)] = 1;
	}
    }

    PutRNGstate();
}

/* 'swapcount' is a C translation of Peter Solymos's R code. It is
 * similar to 'swap', but can swap > 1 values and so works for
 * quantitative (count) data.
 */


/* 'isDiag' is a utility function for 'swapcount' to find the largest
 * value that can be swapped and whether in diagonal or antidiagonal
 * way. The input is a 2x2 submatrix 'sm'.
*/

int isDiag(int *sm, int *change)
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

int isDiagFill(int *sm)
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

void swapcount(int *m, int *nr, int *nc, int *thin)
{
    int row[2], col[2], k, ij[4], changed, 
	pm[4] = {1, -1, -1, 1} ;
    int sm[4], ev;

    GetRNGstate();

    changed = 0;
    while (changed < *thin) {
	/* Select a random 2x2 matrix*/
	i2rand(row, *nr - 1);
	i2rand(col, *nc - 1);
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
    }

    PutRNGstate();
}

/* rswapcount for "reducing swap of count data" is a minor variant of
 * swapcount, but it tries to change the fill: you first make a matrix
 * with correct marginal totals, but possibly wrong fill and then you
 * run this to change the fill while maintaining the totals. The idea
 * is similar as quasiswap for presence/absence data.
 */

void rswapcount(int *m, int *nr, int *nc, int *mfill)
{
    int row[2], col[2], i, k, ij[4], n, change, cfill,
       pm[4] = {1, -1, -1, 1} ;
    int sm[4], ev;

    /* Get the current fill 'cfill' */
    n = (*nr) * (*nc);
    for (i = 0, cfill=0; i < n; i++) {
	if (m[i] > 0) 
	    cfill++;
    }
 
    GetRNGstate();

    /* Loop while fills differ */
    while (cfill != *mfill) {
	/* Select a random 2x2 matrix*/
	i2rand(row, *nr - 1);
	i2rand(col, *nc - 1);
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
    }
    PutRNGstate();
}

/* 'isDiagSimple' needed for 'abuswap' */

int isDiagSimple(double *sm)
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

void abuswap(double *m, int *nr, int *nc, int *thin, int *direct)
{
    int row[2], col[2], k, ij[4], changed, ev;
    double sm[4];

    GetRNGstate();

    changed = 0;
    while (changed < *thin) {
	/* Select a random 2x2 matrix*/
	 i2rand(row, *nr - 1);
	 i2rand(col, *nc - 1);
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
    }
    
    PutRNGstate();
}

#undef IRAND
#undef INDX
