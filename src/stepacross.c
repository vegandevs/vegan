/* stepacross: Extended dissimilarities according to 
 *
 *         G. De'ath, Plant Ecology 144, 191-199; 1999 
 * 
 * 1. Find all dissimilarities above a threshold and make them NA.
 * 2. See if you can find a stepacross  point k so that d[ij] = d[ki] +
 *    d[kj] is not NA. Take the shortest of these dissimilarities as new d[ij].
 * 3. Repeat step 2 with updated dissimilarity matrix until no NA are
 *    left, or you cannot find any new non-NA steps. The latter case is
 *    an error: Disconnected data set.
 *
 * The function inspects only dissimilarities originally marked as
 * NA. It does not update paths involving sites i and j once d[ij] was
 * changed. Neither does it try to find shortest path for
 * dissimilarities below the threshold, although such could be found
 * in semi-metrics such as Bray-Curtis index.
 *
 * (C) 2003 Jari Oksanen. You are allowed to use this software under GPL2.
 */

#include <R.h>

#define EPS (1e-6)
#define IND(N,a,b) (N)*(a) - (a)*((a)+1)/2 + (b) - (a) - (1)

void C_stepacross(double *dist, int *n, double *toolong, int *trace)
{
     int i, j, k, nacount, oldcount, ind, ki, kj, inew, *newind, ndist;
     double stepdis, steptry, *newdist, limit;

     limit = *toolong - EPS;
     ndist = (*n) * ((*n) - 1) / 2;
     if (*toolong > 0)
	     for (i = 0; i < ndist; i++)
		     if (dist[i] >= limit) 
			     dist[i] = NA_REAL;

     for (i = 0, nacount = 0; i < ndist; i++)
	     if (ISNA(dist[i]))
		     nacount++;

     if (*trace)
	  Rprintf("Too long or NA distances: %d out of %d (%.1f%%)\n", 
		  nacount, ndist, 100.0*nacount/ndist);

     newdist = (double *) R_alloc(nacount, sizeof(double));
     newind = (int *) R_alloc(nacount, sizeof(int));

     while (nacount > 0) {
	  if (*trace)
	       Rprintf("Stepping across %d dissimilarities...\n", nacount);
	  oldcount = nacount;
	  inew = 0;
	  for (i = 0, ind = 0; i < *n; i++) {
	       for (j = i + 1; j < *n; j++) {
		    if (ISNA(dist[ind])) {
			 stepdis = DOUBLE_XMAX;
			 for (k = 0; k < *n; k++) {
			      if (k == i || k == j) continue;
			      ki = (k > i) ? IND(*n, i, k) : IND(*n, k, i);
			      if (ISNA(dist[ki])) continue;
			      kj = (k > j) ? IND(*n, j, k) : IND(*n, k, j);
			      if (ISNA(dist[kj])) continue;
			      steptry = dist[ki] + dist[kj];
			      if (steptry < stepdis)
				   stepdis = steptry;
			 }
			 if (stepdis < DOUBLE_XMAX) {
			      newdist[inew] = stepdis;
			      newind[inew] = ind;
			      nacount--;
			      inew++;
			 }
		    }
		    ind++;
	       }
	  }
	  if (oldcount == nacount) {
	       warning("Disconnected data: Result will contain NAs");
	       break;
	  }
	  for (k = 0; k < inew; k++)
	       dist[newind[k]] = newdist[k];
     }
}

/* stepabyss and visitabyss implement depth-first search for
 * connectivity of dissimilarity matrix with threshold
 * limit-EPS. Algorithm is directly taken from Sedgewick (1990),
 * "Algorithms in C," pages 423-427. Macro IND and constant EPS were
 * defined above.
 */

static void visitabyss(int k, int id, int *val, int n, double *dist)
{
     int t, ki;

     val[k] = id;
     for (t = 0; t < n; t++) {
	  if (k == t) continue;
	  ki = (t > k) ? IND(n, k, t) : IND(n, t, k);
	  if (!ISNA(dist[ki]))
	       if (val[t] == 0)
		    visitabyss(t, id, val, n, dist);
     }
}

void stepabyss(double *dist, int *n, double *toolong, int *val)
{
     int k, id, ndist;
     double limit;

     limit = *toolong - EPS;
     ndist = (*n) * ((*n) - 1)/2;
     if (*toolong > 0)
	     for (k = 0; k < ndist; k++)
		     if (dist[k] >= limit)
			     dist[k] = NA_REAL;

     for (k = 0; k < *n; k++)
	  val[k] = 0;

     for (k = 0, id = 0; k < *n; k++)
	  if (val[k] == 0) 
	       visitabyss(k, ++id, val, *n, dist);
}

/* Function dykstrapath impelements Dijkstra's shortest path algorithm
 * for graph traversal. Return matrix outval contains shortest path
 * distances.  Function matrixpfs() in Sedgewick, 1990, page 466
 * (priority first search for dense graphs).
 */

#define UNSEEN 1e8

void dykstrapath(double *dist, int *n, double *toolong, int *trace, 
		 double *outval)
{
     int i, k, t, min = 0, ki, ndist, nacount = 0;
     double priority, *val, limit;

     /* val[*n] is a sentinel: length needs be n+1
      */
     val = (double *) R_alloc((*n) + 1, sizeof(double));

     limit = *toolong - EPS;
     ndist = (*n) * ((*n) - 1)/2;
     if (*toolong > 0)
	     for(i = 0; i < ndist; i++)
		     if (dist[i] >= limit)  
			     dist[i] = NA_REAL;
     
     if (*trace) {
	     for (i = 0, nacount = 0; i < ndist; i++)
		     if (ISNA(dist[i]))
			     nacount++;
	     Rprintf("Too long or NA distances: %d out of %d (%.1f%%)\n", 
		     nacount, ndist, 100.0 * nacount/ndist); 
	     Rprintf("Stepping across %d dissimilarities...\n", ndist);
     }

     for (i = 0; i < *n; i++) {
	  for (k = 0; k <= *n; k++)
	       val[k] = -UNSEEN;
	  
	  val[*n] = -(UNSEEN + 1);
	  
	  for (k = i; k != *n; k = min, min = *n) {
	       val[k] = -val[k];
	       if (val[k] == UNSEEN) 
		    val[k] = 0;
	       for (t = 0; t < *n; t++)
		    if (val[t] < 0) {
			 ki = (t > k) ? IND(*n, k, t) : IND(*n, t, k);
			 priority = val[k] + dist[ki];
			 if (!ISNA(priority) && (val[t] < -priority)) {
			      val[t] = -priority;
			 }
			 if (val[t] > val[min])
			      min = t;
		    }
	  }
	  ki = IND(*n, i, i+1);
	  for (t = i + 1; t < *n; t++)
	       outval[ki + t - i - 1] = val[t];
     }
     for (i = 0, nacount = 0; i < ndist; i++)
	  if (ISNA(dist[i]) && outval[i] == 0) {
	       outval[i] = NA_REAL;
	       nacount++;
	  }
     if (nacount)
	  warning("Disconnected data: Result will contain NAs"); 
}

/* Function primtree finds *a* minimum spanning tree for dist. Vector
 * dad returns the vertex joined to point i (the first item is 0: it
 * is the link from 0 to 0) and val the length of corresponding
 * edge. If the distance matrix is disconnected at level toolong,
 * returns a minimum spanning forest where the missing links are
 * marked as NA.  Uses Prim's method (priority-first search for a
 * dense graph) of Sedgewick, 1990, page 466, function
 * matrixpfs(). --- I didn't intend to have minimum spanning tree,
 * because I think it is pretty useless, but this just appeared when I
 * tried to implement Dijkstra's shortest path algorithm above.
 */

void primtree(double *dist, double *toolong, int *n,  double *val, int *dad)
{
     int k, t, min = 0, ki, ndist;
     double priority;

     ndist = (*n) * ((*n) - 1)/2;
     if (*toolong > 0)
	  for (k = 0; k < ndist; k++)
	       if (dist[k] >= *toolong - EPS)
		    dist[k] = NA_REAL;

     for (k = 0; k <= *n; k++) {
	  dad[k] = NA_INTEGER;
	  val[k] = -UNSEEN;
     }

     val[*n] = -(UNSEEN + 1);

     for (k = 0; k != *n; k = min, min = *n) {
	  val[k] = -val[k];
	  if (val[k] == UNSEEN) 
	       val[k] = 0;
	  for (t = 0; t < *n; t++)
	       if (val[t] < 0) {
		    if (t == k) continue;
		    ki = (t > k) ? IND(*n, k, t) : IND(*n, t, k);
		    priority = dist[ki];
		    if (!ISNA(priority) && (val[t] < -priority)) {
			 val[t] = -priority;
			 dad[t] = k;
		    }
		    if (val[t] > val[min])
			 min = t;
	       }
     }
}
#undef UNSEEN
