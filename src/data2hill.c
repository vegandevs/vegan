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


