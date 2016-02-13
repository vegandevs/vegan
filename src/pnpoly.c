/**********************************************************'
 * pnpoly.c: Find Points in a Polygon
 * ********************************************************/

/* Function pnpoly is lightly edited from the function of Wm. Randolph
 * Franklin by Jari Oksanen to work in R and for a vector of points.
 * The source of this code was at
 * http://www.ecse.rpi.edu/Homepages/wrf/Research/Short_Notes/pnpoly.html,
 * and the original copyright notice is below.  This version is Friday, 
 * October 13, 2006. */

/* Motivation for yet another point in polygon function: There was
 * none in base R, and I wanted to liberate ordisurf from dependence
 * on non-standard R packages. I think there should be such a function
 * connected with chull (in.chull) to make this unnecessary in
 * vegan. */

/* License to Use */

/* Copyright (c) 1970-2003, Wm. Randolph Franklin */

/* Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use, copy,
 * modify, merge, publish, distribute, sublicense, and/or sell copies
 * of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions: */

/*    Redistributions of source code must retain the above copyright
         notice, this list of conditions and the following
         disclaimers. */
/*    Redistributions in binary form must reproduce the above
         copyright notice in the documentation and/or other materials
         provided with the distribution. */
/*    The name of W. Randolph Franklin may not be used to endorse or
         promote products derived from this Software without specific
         prior written permission.  */

/* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
 * BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
 * ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
 * CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.  */

void pnpoly(int *npol, double *xp, double *yp, int *np, double *x, double *y, 
	    int *c)
{
     int i, j, k;
 
     for (k = 0; k < *np; k++)
	  c[k] = 0;
     
     for (k = 0; k < *np; k++) {
	  for (i = 0, j = *npol-1; i < *npol; j = i++) {
	       if ((((yp[i] <= y[k]) && (y[k] < yp[j])) ||
		    ((yp[j] <= y[k]) && (y[k] < yp[i]))) &&
		   (x[k] < (xp[j] - xp[i]) * (y[k] - yp[i]) / (yp[j] - yp[i]) 
		    + xp[i]))
		    c[k] = !c[k];
	  }
     }
}
