/* ============================================================ */
/*                                                              */
/*    pchip: Piecewice Cubic Hermite Interpolating Polynomial   */
/*                                                              */
/* ============================================================ */
/*
*/
/*
** phi relevance function.
**  - Cubic Hermite Spline
** Rita P. Ribeiro
*/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "Boolean.h"
#include "pchip.h"

int findInterval2(double *xt, int n, double x,
		  Rboolean rightmost_closed, Rboolean all_inside,
		  Rboolean left_open, // <- new in findInterval2()
		  int ilo, int *mflag)
{
    int istep, middle, ihi;

/* computes  `left' := max( i ; 1 <= i <= n   &&  xt[i] <= x )  .
 ******  i n p u t  ******
  xt	numeric vector of length  n , assumed to be nondecreasing
  n	length(xt)
  x	the point whose location with respect to the sequence  xt  is
	to be determined.
  rightmost_closed {logical} indicating if the rightmost xt[] interval
	should be closed, i.e. result := n-1 if x == x[n]
	(when left_open, the *leftmost* interval should be closed.)
  all_inside {logical} indicating if result should be coerced
		to lie inside {1, n-1}
  left_open  {logical} use intervals (s, t] instead of [s, t)
  ilo   typically the result of the last call to findInterval(.)
	`ilo' used to be a static variable (in Fortran) which is not
	desirable in R anymore (threads!).
	Instead, you *should* use a reasonable value, in the first call.
 ******  o u t p u t  ******
  left, mflag  both integers, whose value is
   0     -1      if            x <  xt[1]
   i      0      if  xt[i]  <= x <  xt[i+1]
   n      1      if  xt[n]  <= x
	in particular,  mflag = 0 is the 'usual' case.  mflag != 0
	indicates that  x  lies outside the halfopen interval
	xt[1] <= y < xt[n] . the asymmetric treatment of the
	interval is due to the decision to make all pp functions cont-
	inuous from the right.
   Note that if all_inside, left is 1 instead of 0 and n-1 instead of n;
   and if rightmost_closed and x == xt[n],  left is    n-1 instead of n.
 ******  m e t h o d  ******
  the program is designed to be efficient in the common situation that
  it is called repeatedly, with  x  taken from an increasing or decreasing
  sequence. this will happen, e.g., when a pp function is to be graphed.
  The first guess for  left  is therefore taken to be the value returned at
  the previous call and stored in the  l o c a l   variable  ilo .
  a first check ascertains that  ilo < n (this is necessary since the
  present call may have nothing to do with the previous call).
  then, if  xt[ilo] <= x < xt[ilo+1], we set  left = ilo
  and are done after just three comparisons.
  otherwise, we repeatedly double the difference  istep = ihi - ilo
  while also moving  ilo  and  ihi  in the direction of  x , until
		      xt[ilo] <= x < xt[ihi] ,
  after which we use bisection to get, in addition, ilo+1 = ihi .
  left = ilo  is then returned.
*/

#define left_boundary  { *mflag = -1; \
	return((all_inside || (rightmost_closed && x == xt[1])) ? 1 : 0); }

#define right_boundary { *mflag = +1;					\
	return((all_inside || (rightmost_closed && x == xt[n]))		\
		? (n - 1) : n); }

#define X_grtr(XT_v) x > XT_v || (!left_open && x >= XT_v)
#define X_smlr(XT_v) x < XT_v ||  (left_open && x <= XT_v)

    if(n == 0) { *mflag = 0 ; return 0; }
    /* 1-indexing : */
    --xt;

    if(ilo <= 0) {
	if (X_smlr(xt[1]))		left_boundary;
	ilo = 1;
    }
    ihi = ilo + 1;
    if (ihi >= n) {
	if (X_grtr(xt[n]))		right_boundary;
	if (n <= 1) /* x < xt[1] */	left_boundary;
	ilo = n - 1;
	ihi = n;
    }

    if (X_smlr(xt[ihi])) {
	if (X_grtr(xt[ilo])) {
	    /* `lucky': same interval as last time */
	    *mflag = 0;	   return ilo;
	}
	/* **** now x < xt[ilo] .	decrease  ilo  to capture  x */
	if(!left_open) for(istep = 1; ; istep *= 2) {
	    ihi = ilo;
	    ilo = ihi - istep;
	    if (ilo <= 1)
		break;
	    if (x >= xt[ilo])		goto L50;
	} else for(istep = 1; ; istep *= 2) {
	    ihi = ilo;
	    ilo = ihi - istep;
	    if (ilo <= 1)
		break;
	    if (x > xt[ilo])		goto L51;
	}
	ilo = 1;
	if (X_smlr(xt[1]))		left_boundary;
    }
    else {
	/* **** now x >= xt[ihi] .	increase  ihi  to capture  x */
	if(!left_open) for(istep = 1; ; istep *= 2) {
	    ilo = ihi;
	    ihi = ilo + istep;
	    if (ihi >= n)
		break;
	    if (x < xt[ihi])		goto L50;
	}
	else for(istep = 1; ; istep *= 2) {
	    ilo = ihi;
	    ihi = ilo + istep;
	    if (ihi >= n)
		break;
	    if (x <= xt[ihi])		goto L51;
	}
	if (X_grtr(xt[n]))		right_boundary;
	ihi = n;
    }

    if (left_open) goto L51; /* There _is_ a path to here, avoiding return and goto */

L50: // ! left_open
    /* **** now xt[ilo] <= x < xt[ihi] . narrow the interval. */
    for(;;) {
	middle = (ilo + ihi) / 2;
	if (middle == ilo) {
	    *mflag = 0;	   return ilo;
	}
	/* note. it is assumed that middle = ilo in case ihi = ilo+1 . */
	if (x >= xt[middle])
	    ilo = middle;
	else
	    ihi = middle;
    }

L51: // left_open
    /* **** now xt[ilo] < x <= xt[ihi] . narrow the interval. */
    for(;;) {
	middle = (ilo + ihi) / 2;
	if (middle == ilo) {
	    *mflag = 0;	   return ilo;
	}
	/* note. it is assumed that middle = ilo in case ihi = ilo+1 . */
	if (x > xt[middle])
	    ilo = middle;
	else
	    ihi = middle;
    }
} /* findInterval2 */

// has been in API -- keep for compatibility:
int findInterval(double *xt, int n, double x,
		 Rboolean rightmost_closed,  Rboolean all_inside, int ilo,
		 int *mflag)
{
    return findInterval2(xt, n, x, rightmost_closed, all_inside, FALSE, ilo, mflag);
}

hermiteSpl *pchip_set(int n,
		      double *x, double *y, double *m) {


  int i;
  double *h, *delta, *new_m;
  hermiteSpl *H;

  if((H = (hermiteSpl *) calloc(1,sizeof(hermiteSpl))) == NULL) perror("pchip.c: memory allocation error");

  // fill-in the struct -- not sure of malloc

  H->npts = n;

  if((H->x = (double *)calloc(n,sizeof(double))) == NULL) perror("pchip.c: memory allocation error");
  if((H->a = (double *)calloc(n,sizeof(double))) == NULL) perror("pchip.c: memory allocation error");
  if((H->b = (double *)calloc(n,sizeof(double))) == NULL) perror("pchip.c: memory allocation error");
  if((H->c = (double *)calloc(n,sizeof(double))) == NULL) perror("pchip.c: memory allocation error");
  if((H->d = (double *)calloc(n,sizeof(double))) == NULL) perror("pchip.c: memory allocation error");

  if((h = (double *)calloc(n,sizeof(double))) == NULL) perror("pchip.c: memory allocation error");
  if((delta = (double *)calloc(n,sizeof(double))) == NULL) perror("pchip.c: memory allocation error");

  //n +1
  memcpy(H->x,x,n*sizeof(double));
  memcpy(H->a,y,n*sizeof(double));

  // auxiliary vectors
  for(i = 0;i < n-1; i++) {
    h[i] = x[i+1] - x[i];
    delta[i] = (y[i+1] - y[i])/ h[i];
  }

  new_m = pchip_slope_monoFC(n, m, delta);

  memcpy(H->b,new_m,n*sizeof(double));

  for(i = 0;i < n-1; i++) {
    H->c[i] = (3 * delta[i] - 2 * new_m[i] - new_m[i+1]) /
      h[i];
    H->d[i] = (new_m[i] - 2 * delta[i] + new_m[i+1]) /
      (h[i] *  h[i]);
  }

  return H;
}

/*
Slopes for shape-preserving Hermite cubic polynomials
 */


/**
 * Modify the slopes  m_k := s'(x_k) using Fritsch & Carlson (1980)'s algorithm
 *
 * @param m  numeric vector of length n, the preliminary desired slopes s'(x_i), i = 1:n
 * @param S the divided differences (y_{i+1} - y_i) / (x_{i+1} - x_i);        i = 1:(n-1)
 * @return m*: the modified m[]'s: Note that m[] is modified in place
 * @author Martin Maechler, Date: 19 Apr 2010
 */
// adapted
double *pchip_slope_monoFC(int n, double *m, double *delta) {

  for(int k = 0; k < n - 1; k++) {
    /* modify both (m[k] & m[k+1]) if needed : */
    double Sk = delta[k];


    int k1 = k + 1;

    if(fabs(Sk) == 0) {
      m[k] = m[k1] = 0.;

    } else {

      double
	alpha = m[k ] / Sk,
	beta  = m[k1] / Sk, a2b3, ab23;

      if(fabs(m[k]) !=0 && alpha < 0) {
	m[k] = -m[k];
	alpha = m[k] / Sk;
      }

      if(fabs(m[k1]) !=0 && beta < 0) {
	m[k1] = -m[k1];
	beta = m[k1] / Sk;
      }

      a2b3 = 2*alpha + beta - 3;
      ab23 = alpha + 2*beta - 3;

      if(a2b3 > 0 && ab23 > 0 &&
	 alpha * (a2b3 + ab23) < a2b3*a2b3) {
	/* we are outside the monotonocity region ==> fix slopes */
	double tauS = 3*Sk / sqrt(alpha*alpha + beta*beta);
	m[k ] = tauS * alpha;
	m[k1] = tauS * beta;

      }
    }
  } /* end for */

  return m;
}


//  Evaluate the cubic polynomial.
//  Find, from the left, the interval that contains or is nearest to xval.
// Check for linear extrapolation
// use cubic Hermite polynomials, even for extrapolation
void  pchip_val(hermiteSpl *H, double xval, int extrapol,
		double *yval) {

  int i = 1, rightmost_closed = 0, all_inside = 0, mfl = 0;
  double s;

  i = findInterval(H->x,H->npts,
		   xval,
		   rightmost_closed,all_inside,i,&mfl);


  // if extrapol is linear
  if(extrapol == 0 && (i == 0 || i == H->npts)) {

    if(i == H->npts) i--;

    *yval = H->a[i] + H->b[i] * (xval - H->x[i]);

    return;
  }


  i--;

  s = (xval - H->x[i]);
  *yval = H->a[i] + s * (H->b[i] +
		     s * (H->c[i] +
			  s * H->d[i]));

}