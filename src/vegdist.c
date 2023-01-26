/*
 * Distance measures for community ecologists.  The measures here were
 * recommended by Peter Minchin, since they have a good rank-order
 * relation with gradient distance.  The standard distances are found
 * in standard R library mva in function dist (distance.c).
 *
 * The calling program defines the index by a number (internally and
 * factually). Extra care is needed to get the numbers right there
 * above when calling the code.
 *
 * Number "99" is an extra case: It is not for vegdist.R, but for
 * something simpler.
 *
 * (C) 2001-2005, Jari Oksanen. You are free to use this code if you
 * accept GPL2.
 *
 * Oct 2003: Added Morisita, Horn-Morisita, "Jaccard", and Mountford.
 * May 2005: Added Raup-Crick.
 */


/* Standard R headers */
#include <R_ext/Arith.h>

#include <R.h>
#include <Rmath.h>
#include <R_ext/Utils.h> /* interrupt */
#include <float.h>
#include <string.h> /* memset */

/* Indices: The numbers here MUST match the order of indices in vegdist.R */

#define MANHATTAN 1
#define EUCLIDEAN 2
#define CANBERRA 3
#define BRAY 4
#define KULCZYNSKI 5
#define GOWER 6
#define MORISITA 7
#define HORN 8
#define MOUNTFORD 9
#define JACCARD 10
#define RAUP 11
#define MILLAR 12
#define CHAO 13
#define GOWERDZ 14
#define CAO 15
#define MAHALANOBIS 16
#define CLARK 17
#define CHISQUARE 18
#define CHORD 19
#define HELLINGER 20
#define AITCHISON 21
#define RAITCHISON 22
#define MATCHING 50
#define NOSHARED 99

/* Distance functions */

/* Manhattan distance: duplicates base R */

static double veg_manhattan(double *x, int nr, int nc, int i1, int i2)
{
     double dist;
     int count, j;
  
     dist = 0.0;
     count = 0;
     for (j=0; j<nc; j++) {
	  if (!ISNAN(x[i1]) && !ISNAN(x[i2])) {
	       dist += fabs( x[i1] - x[i2] );
	       count++;
	  }
	  i1 += nr;
	  i2 += nr;
     }
     if (count == 0) dist = NA_REAL;
     return dist;
}

/* Gower is like Manhattan, but data were standardized to range 0..1
 * for rows before call and dist is divided by the number of non-zero
 * pairs. There is an alternative implementation in cluster package.
 * Some extra manipulations are needed in the calling R function.
 */

static double veg_gower(double *x, int nr, int nc, int i1, int i2)
{
     double dist;
     int count, j;
  
     dist = 0.0;
     count = 0;
     for (j=0; j<nc; j++) {
	  if (!ISNAN(x[i1]) && !ISNAN(x[i2])) {
	       dist += fabs( x[i1] - x[i2] );
	       count++;
	  }
	  i1 += nr;
	  i2 += nr;
     }
     if (count == 0) dist = NA_REAL;
     dist /= (double) count;
     return dist;
}

/* Identical to veg_gower except skipping count update for double
 * zeros.  Gower 1971 proposed skipping double zeros for
 * presence/absence data (but not for quantitative data), and many
 * people think this should always be done in Gower distance. With
 * presence/absence data this gives Jaccard of binary data.
 */

static double veg_gowerDZ(double *x, int nr, int nc, int i1, int i2)
{
     double dist;
     int count, j;
  
     dist = 0.0;
     count = 0;
     for (j=0; j<nc; j++) {
	  if (!ISNAN(x[i1]) && !ISNAN(x[i2])) {
	      if (x[i1] > 0 || x[i2] > 0) {
		  dist += fabs( x[i1] - x[i2] );
		  count++;
	      }
	  }
	  i1 += nr;
	  i2 += nr;
     }
     if (count == 0) dist = NA_REAL;
     dist /= (double) count;
     return dist;
}

/* Euclidean distance: duplicates base R. This function can be (and will be)
* used for many named distance indices with transforming & standardizing input
* data in the calling R code (vegdist): Chord, Hellinger, Chi-square,
* Aitchison, Mahalanobis -- see veg_distance() below for actual list at the
* moment.
*/

static double veg_euclidean(double *x, int nr, int nc, int i1, int i2)
{
     double dist, dev;
     int count, j;

     count = 0;
     dist = 0.0;
     for (j=0; j<nc; j++) {
	  if (!ISNAN(x[i1]) && !ISNAN(x[i2])) {
	       dev = x[i1] - x[i2];
	       dist += dev*dev;
	       count++;
	  }
	  i1 += nr;
	  i2 += nr;
     }
     if (count == 0) return NA_REAL;
     return sqrt(dist);
}

/* Chord & Hellinger distances can be calculated after appropriate
* standardization in veg_euclidean, but it is also possible to use
* this function for non-standardized (Chord) data or square root
* transformed (Hellinger) data. It may be that veg_euclidean() is
* numerically more stable: the current function collects sum of
* squares and cross-products which can have loss of precision. See
* driver routine veg_distance() to see which function currently is
* used for Chord and Hellinger distances.
*/

static double veg_chord(double *x, int nr, int nc, int i1, int i2)
{
    double dist, cp = 0.0, ss1 = 0.0, ss2 = 0.0;
    int count = 0, j;

    for (j = 0; j < nc; j++) {
	if (!ISNAN(x[i1]) && !ISNAN(x[i2])) {
	    count++;
	    cp += x[i1] * x[i2];  /* cross products */
	    ss1 += x[i1] * x[i1]; /* sum of squares */
	    ss2 += x[i2] * x[i2];
	}
	i1 += nr;
	i2 += nr;
    }

    if (count == 0) return NA_REAL;
    dist = 2.0 * (1.0 - cp/sqrt(ss1 * ss2));
    return sqrt(dist);
}

/* It may not make sense to use direct formula to Chord distance as
* this involves getting sums of squares which may be numerically
* unstable. It makes even less sense to find Hellinger distance via
* Chord of sqrt-transformed data (sqrt(2)*sqrt(2) != 2), but it makes
* perfect sense to calculate Hellinger distance directly: this is more
* stable than working via sqrt(x) data.
*/

static double veg_hellinger(double *x, int nr, int nc, int i1, int i2)
{
    double dist, cp = 0.0, ss1 = 0.0, ss2 = 0.0;
    int count = 0, j;

    for (j = 0; j < nc; j++) {
	if (!ISNAN(x[i1]) && !ISNAN(x[i2])) {
	    count++;
	    cp += sqrt(x[i1] * x[i2]);
	    ss1 += x[i1];
	    ss2 += x[i2];
	}
	i1 += nr;
	i2 += nr;
    }

    if (count == 0) return NA_REAL;
    dist = 2.0 * (1.0 - cp/sqrt(ss1 * ss2));
    return sqrt(dist);
}

/* Canberra distance: duplicates R base, but is scaled into range
 * 0...1  
*/

static double veg_canberra(double *x, int nr, int nc, int i1, int i2)
{
     double numer, denom, dist;
     int count, j;

     count = 0;
     dist = 0.0;
     for (j=0; j<nc; j++) {
	  if (!ISNAN(x[i1]) && !ISNAN(x[i2])) {
	       if (x[i1] != 0 || x[i2] != 0) {
		    count++;
		    denom = fabs(x[i1]) + fabs(x[i2]);
		    if (denom > 0.0) {
			 numer = fabs(x[i1] - x[i2]);
			 dist += numer/denom;
		    }
		    else {
			 dist += R_PosInf;
		    }
	       }
	  }
	  i1 += nr;
	  i2 += nr;
     }
     if (count == 0) return NA_REAL;
     dist /= (double)count;
     return dist;
}

/* Clark's (1952) coefficient of divergence is similar to Canberra
 * index, but uses squared terms instead of minimum terms. This is one
 * of the indices discussed by Legendre & De Caceres, Ecol Lett 16,
 * 951-963 (2012) and therefore implemented here. Discussed as index
 * D11 in Legendre & Legendre 2012, Numerical Ecology.
 */

static double veg_clark(double *x, int nr, int nc, int i1, int i2)
{
    double t1, denom, dist;
    int count, j;

    count = 0;
    dist = 0.0;
    for (j=0; j<nc; j++) {
	if (!ISNAN(x[i1]) && !ISNAN(x[i2])) {
	    if (x[i1] != 0 || x[i2] != 0) {
		count++;
		denom = x[i1] + x[i2];
		if (denom > 0.0) {
		    t1 = (x[i1] - x[i2])/denom;
		    dist += t1 * t1;
		}
		else {
		    dist += R_PosInf;
		}
	    }
	}
	i1 += nr;
	i2 += nr;
    }
    if (count == 0) return NA_REAL;
    dist /= (double)count;
    return sqrt(dist);
}

/*  Bray-Curtis and Jaccard indices:
 *
 * Jaccard = (2 * Bray)/(1 + Bray). If Jaccard is requested, Bray is
 * calculated in this function and it is left as the task of the
 * caller to translate this into Jaccard. Actually, Jaccard is
 * redundant, but since people ask for Jaccard, they get it.
 */

static double veg_bray(double *x, int nr, int nc, int i1, int i2)
{
     double dist, total;
     int count, j;
  
     total = 0.0;
     count = 0;
     dist = 0;
     for (j=0; j<nc; j++) {
	  if (!ISNAN(x[i1]) && !ISNAN(x[i2])) {
	       dist += fabs(x[i1] - x[i2]);
	       total += x[i1] + x[i2];
	       count++;
	  }
	  i1 += nr;
	  i2 += nr;
     }
     if (count==0) return NA_REAL;
     dist /= total;
     return dist;
}

/* Kulczynski index */

static double veg_kulczynski(double *x, int nr, int nc, int i1, int i2)
{
     double sim, dist, t1, t2;
     int count, j;

     t1 = 0.0;
     t2 = 0.0;
     count = 0;
     sim = 0.0;
     for (j=0; j<nc; j++) {
	  if (!ISNAN(x[i1]) && !ISNAN(x[i2])) {
	       sim += (x[i1] < x[i2]) ? x[i1] : x[i2] ;
	       t1 += x[i1];
	       t2 += x[i2];
	       count++;
	  }
	  i1 += nr;
	  i2 += nr;
     }
     if (count==0) return NA_REAL;
     dist = 1 - sim/t1/2 - sim/t2/2;
     if (dist < 0)
	  dist = 0;
     return dist;
}

/* Morisita index.  Can only be used with integer data, and may still
 * fail with unfortunate pairs of species occurring only once.
 */

static double veg_morisita(double *x, int nr, int nc, int i1, int i2)
{
     double sim, dist, t1, t2, tlam1, tlam2;
     int count, j;

     t1 = 0.0;
     t2 = 0.0;
     count = 0;
     sim = 0.0;
     tlam1 = 0.0;
     tlam2 = 0.0;
     for (j=0; j<nc; j++) {
	  if (!ISNAN(x[i1]) && !ISNAN(x[i2])) {
	       sim += x[i1]*x[i2];
	       t1 += x[i1];
	       t2 += x[i2];
	       tlam1 += x[i1]*(x[i1] - 1);
	       tlam2 += x[i2]*(x[i2] - 1);
	       count++;
	  }
	  i1 += nr;
	  i2 += nr;
     }
     if (count==0) return NA_REAL;
     dist = 1 - 2*sim/(tlam1/t1/(t1-1) + tlam2/t2/(t2-1))/t1/t2;
     if (dist < 0)
	  dist = 0;
     return dist;
}

/* Horn-Morisita index */

static double veg_horn(double *x, int nr, int nc, int i1, int i2)
{
     double sim, dist,  t1, t2, sq1, sq2;
     int count, j;

     t1 = 0.0;
     t2 = 0.0;
     count = 0;
     sim = 0.0;
     sq1 = 0.0;
     sq2 = 0.0;
     for (j=0; j<nc; j++) {
	  if (!ISNAN(x[i1]) && !ISNAN(x[i2])) {
	       sim += x[i1]*x[i2];
	       t1 += x[i1];
	       t2 += x[i2];
	       sq1 += x[i1]*x[i1];
	       sq2 += x[i2]*x[i2];
	       count++;
	  }
	  i1 += nr;
	  i2 += nr;
     }
     if (count==0) return NA_REAL;
     dist = 1 - 2*sim/(sq1/t1/t1 + sq2/t2/t2)/t1/t2;
     if (dist < 0)
	 dist = 0;
     return dist;
}


/* Mountford index theta is defined as the root of an exponential
 * equation (function mount_fun below). The value of theta is found
 * using Newton method (derivatives in mount_der) in veg_mountford.
 * The result is divided by log(2) to put dissimilarities into
 * conventional range 0...1.
 */

#define MAXIT 20
#define EPS 1e-12
#define TOL 1e-5

static double mount_fun(double theta, double j, double a, double b)
{
     return(exp(theta*a) + exp(theta*b) - exp(theta*(a+b-j)) - 1);
}

static double mount_der(double theta, double j, double a, double b)
{
     return(a*exp(theta*a) + b*exp(theta*b) - (a+b-j)*exp(theta*(a+b-j)));
}

static double veg_mountford(double *x, int nr, int nc, int i1, int i2)
{
     double dist, oldist, A,  B, J;
     int sim, t1, t2, j, count;
     
     sim = 0;
     t1 = 0;
     t2 = 0;
     count = 0;
     for (j = 0; j < nc; j++) {
	  if (!ISNAN(x[i1]) && !ISNAN(x[i2])) {
	       if (x[i1] > 0.0 && x[i2] > 0.0)
		    sim++;
	       if (x[i1] > 0)
		    t1++;
	       if (x[i2] > 0)
		    t2++;
	       count++;
	  }
	  i1 += nr;
	  i2 += nr;
     }
     if (count == 0) return NA_REAL;
     if (t1 == 0 || t2 == 0)
	  dist = NA_REAL;
     else if (sim == 0)
	  dist = 0;
     else if (sim == t1 || sim == t2)
	  dist = M_LN2;
     else {
	  J = (double)(sim);
	  A = (double)(t1);
	  B = (double)(t2);
	  dist = 2*J/(2*A*B - (A+B)*J);
	  for (j = 0; j < MAXIT; j++) {
	       oldist = dist;
	       dist -= mount_fun(dist, J, A, B)/mount_der(dist, J, A, B);
	       if(fabs(oldist - dist)/oldist < TOL || fabs(oldist - dist) < EPS) 
		    break;
	  }
     }
     /* index in range 0..log(2) */
     return M_LN2 - dist;
}

#undef MAXIT
#undef EPS
#undef TOL

/* Raup-Crick dissimilarity: R code supplied by 
 * Michael.Bedward@environment.nsw.gov.au.  
 *
 * Here his original comments:
 *
 * "Attached is a little function to calculate the probabilistic
 * Raup-Crick dissimilarity metric for presence-absence data.  Rather
 * than the permutation procedure used in the original reference (Raup
 * & Crick 1979 Paleontology, as related by Legendre and Legendre in
 * Numerical Ecology), this function uses phyper() for a faster and
 * more precise calculcation.  I subsequently found that the same
 * (obvious) idea is in the literature under a variety of other names
 * (or sometimes no name).
 *
 * Compared to other metrics for p/a data, Raup-Crick seems to be very
 * robust for small samples."
 *
 * This is a direct port from Bedward's R to C (Jari Oksanen, May 2005).
 */

static double veg_raup(double *x, int nr, int nc, int i1, int i2)
{
	double dist, J, A, B;
	int sim, t1, t2, j, count;
	
	sim = 0;
	t1 = 0;
	t2 = 0;
	count = 0;
	for (j = 0; j < nc; j++) {
		if (!ISNAN(x[i1]) && !ISNAN(x[i2])) {
			if (x[i1] > 0.0 && x[i2] > 0.0)
				sim++;
			if (x[i1] > 0)
				t1++;
			if (x[i2] > 0)
				t2++;
			count++;
		}
		i1 += nr;
		i2 += nr;
	}
	if (count == 0) return NA_REAL; 
	J = (double) (sim - 1);
	A = (t1 < t2) ? (double) t1 : (double) t2;
	B = (t1 < t2) ? (double) t2 : (double) t1;
	dist = 1 - phyper(J, A, (double) count - A, B, 1, 0);
	return dist;
}

/* "Millar dissimilarity" is unpublished.  I found this in the lecture
 * notes of Marti Anderson over the internet, and she attributes this
 * idea to her colleague Russell Millar.  The index is basically
 * binomial deviance under H0 that species are equally common in the
 * two compared communities.  This could be easily generalized over
 * to, say, Poisson case.
 */

static double veg_millar(double *x, int nr, int nc, int i1, int i2)
{
     double dist, t1, t2, nk, lognk;
     int count, j;
  
     count = 0;
     dist = 0;
     for (j=0; j<nc; j++, i1 += nr, i2 += nr) {
	  if (!ISNAN(x[i1]) && !ISNAN(x[i2])) {
	       nk = x[i1] + x[i2];
	       if (nk == 0) continue;
	       lognk = log(nk);
	       t1 = (x[i1] > 0) ? x[i1] * (log(x[i1]) - lognk) : 0;
	       t2 = (x[i2] > 0) ? x[i2] * (log(x[i2]) - lognk) : 0;
	       dist += (t1 + t2 + nk * M_LN2)/nk;
	       count++;
	  }
     }
     if (count==0) return NA_REAL;
     if (dist < 0)
	 dist = 0;
     return dist;
}

/* Chao's index (Ecol. Lett. 8, 148-159; 2005) tries to take into
 * account the number of unseen shared species using Chao's method for
 * estimating the number of unseen species. June 2006.
 */

static void chaoterms(double *x, int nr, int nc, int i1, int i2,
		      double *U, double *V)
{
    double ionce = 0, itwice = 0, jonce = 0, jtwice = 0, itot = 0, jtot = 0,
	ishare = 0, jshare = 0, ishar1 = 0, jshar1 = 0;
    int count = 0, j;
  
    for (j=0; j<nc; j++) {
	if (!ISNAN(x[i1]) && !ISNAN(x[i2])) {
	    count++;
	    itot += x[i1];
	    jtot += x[i2];
	    if (x[i1] > 0 && x[i2] > 0) {
		ishare += x[i1];
		jshare += x[i2];
		if (fabs(x[i2] - 1) < 0.01) {
		    ishar1 += x[i1];
		    jonce += 1;
		} else if (fabs(x[i2] - 2) < 0.01) {
		    jtwice += 1;
		}
		if (fabs(x[i1] - 1) < 0.01) { 
		    jshar1 += x[i2];
		    ionce += 1;
		} else if (fabs(x[i1] - 2) < 0.01) {
		    itwice += 1;
		}
	    }
	}
	i1 += nr;
	i2 += nr;
    }
    if (count==0) {
	*U = NA_REAL;
	*V = NA_REAL;
	return;
    }
    *U = ishare/itot;
    if (ishar1 > 0) {
	if (jonce < 1) jonce = 1; /* Never true if got here? */
	if (jtwice < 1) jtwice = 1;
	*U += (jtot - 1)/jtot * jonce/jtwice/2.0 * ishar1/itot;
    }
    if (*U > 1) *U = 1;
    *V = jshare/jtot;
    if (jshar1 > 0) {
	if (ionce < 1) ionce = 1; /* This never true? */
	if (itwice < 1) itwice = 1;
	*V += (itot - 1)/itot * ionce/itwice/2.0 * jshar1/jtot;
    }
    if (*V > 1) *V = 1;
    return;
}

static double veg_chaojaccard(double *x, int nr, int nc, int i1, int i2)
{
    double dist, U, V;
    chaoterms(x, nr, nc, i1, i2, &U, &V);
    if (U <= 0 || V <= 0)
	dist = 1;
    else
	dist = 1 - U*V/(U + V - U*V);
    if (dist < 0)
	dist = 0;
    return dist;
}

/* veg_cao implements Cao index (CYd) of Cao Y, Williams WP, Bark AW:
 *   Water Envir Res 69, 95-106; 1997. Anderson MJ & Thompson AA: Ecol
 *   Appl 14, 1921-1935; 2004 use different but equal formulation.
 */

static double veg_cao(double *x, int nr, int nc, int i1, int i2)
{
     double dist, x1, x2, t1, t2;
     int count, j;
  
     count = 0;
     dist = 0;
     for (j=0; j<nc; j++, i1 += nr, i2 += nr) {
	  if (!ISNAN(x[i1]) && !ISNAN(x[i2])) {
	       /* skip the rest of the loop if both species are
		  absent */
	       if (x[i1] == 0 && x[i2] == 0) continue;
	       /* Cao uses arbitrary value of 0.1 for zeros to avoid
		  log(0). Obviously this indicates the use of counts
		  (integer), but we accept non-integer data (with a
		  warning in R) and put the truncation to the same 0.1
		  to avoid discontinuities with non-integer data */
	       x1 = (x[i1] < 0.1) ? 0.1 : x[i1];
	       x2 = (x[i2] < 0.1) ? 0.1 : x[i2];
	       t1 = x1 + x2;
	       /* Cao et al. used log10, but we do not and so our
		  results are log(10) = 2.302585 times higher */
	       t2 = x1 * log(x2) + x2 * log(x1);
	       dist += log(t1) - M_LN2 - t2/t1;
	       count++;
	  }
     }
     if (count==0) return NA_REAL;
     if (dist < 0)
	 dist = 0;
     dist /= (double)count;
     return dist;
}


/* veg_noshared is not a proper dissimilarity index, but a pretty
 * useless helper function. It returns 1 when there are no shared
 * species, and 0 if two sites have at least one shared species, and
 * it stops looping after finding the first shared species (hence it
 * should be fast).
 */

static double veg_noshared(double *x, int nr, int nc, int i1, int i2)
{
     double dist;
     int j, count;
     dist = 1;
     count = 0;
     for (j = 0; j<nc; j++) {
	  if (!ISNAN(x[i1]) && !ISNAN(x[i2])) {
	       count++;
	       if (x[i1] > 0 && x[i2] > 0) {
		    dist = 0;
		    break;
	       }
	  }
	  i1 += nr;
	  i2 += nr;
     }
     if (count == 0) return NA_REAL;
     return(dist);
}

/* Simple matching coefficient.  This is not to be called from
 * vegdist, but must be called separately.
 */

static double veg_matching(double *x, int nr, int nc, int i1, int i2)
{
     double dist;
     int j, count, matches;
     matches = 0;
     count = 0;
     for (j = 0; j<nc; j++) {
	  if (!ISNAN(x[i1]) && !ISNAN(x[i2])) {
	       count++;
	       if (x[i1] == x[i2])
		    matches++;
	  }
	  i1 += nr;
	  i2 += nr; 
     }
     if (count == 0) return NA_REAL;
     dist = 1.0 - (double) matches / (double) count;
     return(dist);
}

/* Driver */

static double (*distfun)(double*, int, int, int, int);

static void veg_distance(double *x, int *nr, int *nc, double *d, int *diag,
			 int *method)
{
    int dc, i, j;
    size_t ij; /* can be a long vector exceeding int maximum */
    switch(*method) {
    case MANHATTAN:
	distfun = veg_manhattan;
	break;
    case EUCLIDEAN:
    case MAHALANOBIS:
    case CHISQUARE:
    case AITCHISON:
    case RAITCHISON:
	distfun = veg_euclidean;
	break;
    case CHORD:
	distfun = veg_chord;
	break;
    case HELLINGER:
	distfun = veg_hellinger;
	break;
    case CANBERRA:
	distfun = veg_canberra;
	break;
    case BRAY:
    case JACCARD:
	distfun = veg_bray;
	break;
    case KULCZYNSKI:
	distfun = veg_kulczynski;
	break;
    case GOWER:
	distfun = veg_gower;
	break;
    case MORISITA:
	distfun = veg_morisita;
	break;
    case HORN:
	distfun = veg_horn;
	break;
    case MOUNTFORD:
	distfun = veg_mountford;
	break;
    case RAUP:
	distfun = veg_raup;
	break;
    case MILLAR:
	distfun = veg_millar;
	break;
    case CHAO:
	distfun = veg_chaojaccard;
	break;
    case GOWERDZ:
	distfun = veg_gowerDZ;
	break;
    case CAO:
        distfun = veg_cao;
        break;
    case CLARK:
	distfun = veg_clark;
	break;
    case MATCHING:
	distfun = veg_matching;
	break;
    case NOSHARED:
	distfun = veg_noshared;
	break;	 
    default:
	error("Unknown distance in the internal C function");
    }

    dc = (*diag) ? 0 : 1;
    ij = 0;
    for (j=0; j <= *nr; j++) {
	if (j % 200 == 199)
	    R_CheckUserInterrupt();
	for (i=j+dc; i < *nr; i++) {
	    d[ij++] = distfun(x, *nr, *nc, i, j);
	}
    }
}


/* .Call interface to veg_distance */

#include <Rinternals.h>

SEXP do_vegdist(SEXP x, SEXP method)
{
    SEXP dist;
    int nr, nc, imethod = asInteger(method), diag=0;
    R_xlen_t ndis;

    nr = nrows(x);
    nc = ncols(x);
    ndis = (R_xlen_t) nr * (nr-1)/2;

    PROTECT(dist = allocVector(REALSXP, ndis));
    if(TYPEOF(x) != REALSXP)
	x = coerceVector(x, REALSXP);
    PROTECT(x);

    veg_distance(REAL(x), &nr, &nc, REAL(dist), &diag, &imethod);

    UNPROTECT(2);
    return dist;
}

/* Minimum terms for designdist. Returns a matrix where the diagonal
 * contains row sums and below-diagonal elements the sums of pairwise
 * minima for rows. Input x must be a matrix.
 */

SEXP do_minterms(SEXP x)
{
    int nr = nrows(x), nc = ncols(x), i, j, k;
    double t1, t2, sum;

    SEXP terms = PROTECT(allocMatrix(REALSXP, nr, nr));
    double *rterm = REAL(terms);
    memset(rterm, 0, nr * nr * sizeof(double));
    if(TYPEOF(x) != REALSXP)
	x = coerceVector(x, REALSXP);
    PROTECT(x);
    double *rx = REAL(x);

    /* end R house keeping, start actual work */
    for(i = 0; i < nr; i++) {
	for(j = i; j < nr; j++) {
	    sum = 0.0;
	    /* allow for NA result and do not check ISNAN */
	    for(k = 0; k < nc; k++) {
		t1 = rx[i + nr*k];
		t2 = rx[j + nr*k];
		sum += (t1 < t2) ? t1 : t2;
	    }
	    /* fill the lower triangle and the diagonal */
	    rterm[j + nr*i] = sum;
	}
    }

    /* back to R house keeping: set dimnames */
    SEXP rnames = getAttrib(x, R_DimNamesSymbol);
    if (!isNull(rnames) && !isNull(VECTOR_ELT(rnames, 0))) {
	SEXP dimnames = PROTECT(allocVector(VECSXP, 2));
	SET_VECTOR_ELT(dimnames, 0, duplicate(VECTOR_ELT(rnames, 0)));
	SET_VECTOR_ELT(dimnames, 1, duplicate(VECTOR_ELT(rnames, 0)));
	setAttrib(terms, R_DimNamesSymbol, dimnames);
	UNPROTECT(1);
    }

    UNPROTECT(2);
    return terms;
}

/* Extract Chao terms U & V for designdist */

SEXP do_chaoterms(SEXP x)
{
    int nr = nrows(x), nc = ncols(x);
    R_xlen_t nterms, i, j, ij;
    nterms = (R_xlen_t) nr * (nr-1)/2;

    if (TYPEOF(x) != REALSXP)
	x = coerceVector(x, REALSXP);
    PROTECT(x);

    SEXP U = PROTECT(allocVector(REALSXP, nterms));
    SEXP V = PROTECT(allocVector(REALSXP, nterms));
    double *ru = REAL(U);
    double *rv = REAL(V);

    /* collect U & V */

    ij = 0;
    for (j = 0; j < nr; j++)
	for (i = j + 1; i < nr; i++) {
	    chaoterms(REAL(x), nr, nc, i, j, ru + ij, rv + ij);
	    ij++;
	}

    /* out */

    SEXP out = PROTECT(allocVector(VECSXP, 2));
    SEXP names = PROTECT(allocVector(STRSXP, 2));
    SET_STRING_ELT(names, 0, mkChar("U"));
    SET_STRING_ELT(names, 1, mkChar("V"));
    setAttrib(out, R_NamesSymbol, names);
    UNPROTECT(1); /* names */
    SET_VECTOR_ELT(out, 0, U);
    SET_VECTOR_ELT(out, 1, V);

    UNPROTECT(4);
    return out;
}
