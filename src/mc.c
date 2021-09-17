/*
    Algorithm for the skewness estimator medcouple (MC)
    --------------------------------------------------
    ( originally matlabmc.c and also  mc/mcrsoft/spmc.c )
*/

#include <stdlib.h>
#include <math.h>
#include "R.h"
#include <inttypes.h>
// -> int64_t

#include "Rmath.h"
/* -> fmax2(.,.) */
#include "Utils.h"

#define _i_whimed_
#define _d_whimed_



#ifdef _d_whimed_

# define _WHIMED_	whimed
# define _WGT_TYPE_	double
# define _WGT_SUM_TYPE_ double
# undef _d_whimed_

#elif defined (_i_whimed_)


# define _WGT_TYPE_	int
# define _WGT_SUM_TYPE_ int64_t
# undef _i_whimed_

#else
# error "must define correct  whimed_  macro !"
#endif


/* Interface routines to be called via .C() and those from API : */
#include "robustbase.h"
/*
  including

 whimed_i(a,iw,n): the weighted high median of an array a of length n,
		   using the positive integer weights iw[].

 * which is in ./wgt_himed.c_templ
 *		 ~~~~~~~~~~~~~~~~~
*/

/* Includes the auxiliary function

   h_kern(a,b, ai,bi,ab, eps):	the values h(a,b)  needed to compute the mc
*/
static
double h_kern(double a, double b, int ai, int bi, int ab, double eps, Rboolean do_scale);

// Called via .C() :
__declspec(dllexport) void mc_C(double *z, size_t in, double *eps, int *iter, double *out, int *scale)
{
    *out = mc_C_d(z, in, eps, iter, *scale);
    return;
}

/* MM:	The tolerance  'eps1' and 'eps2' can now be passed from R;
 *	the original code had only one 'eps' for both and hardcoded
 *	   eps =  0.0000000000001;  (== 1e-13 )
 *
 * MK:  eps1: for (relative) "equality" checks
 *      eps2: used to check for over- and underflow, respectively
 *      therefore I suggest eps1 = DBL_EPS and eps2 = DBL_MIN
 */

double whimed_i(double *a, int *w, int n,
		double* a_cand, double *a_srt, int* w_cand)
{

/*
  Algorithm to compute the weighted high median in O(n) time.

  The whimed is defined as the smallest a[j] such that the sum
  of the weights of all a[i] <= a[j] is strictly greater than
  half of the total weight.

  Arguments:

  a: double array containing the observations
  n: number of observations
  w: array of (int/double) weights of the observations.
*/
    int i;
    /* sum of weights: `int' do overflow when  n ~>= 1e5 */
    _WGT_SUM_TYPE_ wleft, wmid, wright, w_tot, wrest;
    double trial;

    w_tot = wrest = 0;
    for (i = 0; i < n; ++i)
	w_tot += w[i];


    if(n == 0) return 0;

/* REPEAT : */
    do {
	int n2 = n/2;/* =^= n/2 +1 with 0-indexing */
	for (i = 0; i < n; ++i)
	    a_srt[i] = a[i];
	rPsort(a_srt, n, n2);
	trial = a_srt[n2];

	wleft = 0;    wmid  = 0;    wright= 0;
	for (i = 0; i < n; ++i) {
	    if (a[i] < trial)
		wleft += w[i];
	    else if (a[i] > trial)
		wright += w[i];
	    else
		wmid += w[i];
	}
	/* wleft = sum_{i; a[i]	 < trial}  w[i]
	 * wmid	 = sum_{i; a[i] == trial}  w[i] at least one 'i' since trial is one a[]!
	 * wright= sum_{i; a[i]	 > trial}  w[i]
	 */
#ifdef DEBUG_whimed
	REprintf(" trial=%-g; w(left|mid|right) = (%g,%g,%g); ", trial,
		 (double)wleft, (double)wmid, (double)wright);
#endif

	int kcand = 0;
	if (2 * (wrest + wleft) > w_tot) {
	    for (i = 0; i < n; ++i) {
		if (a[i] < trial) {
		    a_cand[kcand] = a[i];
		    w_cand[kcand] = w[i];	++kcand;
		}
	    }
	}
	else if (2 * (wrest + wleft + wmid) <= w_tot) {
	    for (i = 0; i < n; ++i) {
		if (a[i] > trial) {
		    a_cand[kcand] = a[i];
		    w_cand[kcand] = w[i];	++kcand;
		}
	    }
	    wrest += wleft + wmid;
#ifdef DEBUG_whimed
	    REprintf(" new wrest = %g; ", (double)wrest);
#endif
	}
	else {
#ifdef DEBUG_whimed
	    REprintf(" -> found! return trial\n");
#endif
	    return trial;
	    /*==========*/
	}
	n = kcand;
#ifdef DEBUG_whimed
	REprintf("  ... and try again with  n:= kcand=%d\n", n);
#endif
	for (i = 0; i < n; ++i) {
	    a[i] = a_cand[i];
	    w[i] = w_cand[i];
	}
    } while(1);

}



// from Rmath

double fmax2(double x, double y)
{
#ifdef IEEE_754
	if (ISNAN(x) || ISNAN(y))
		return x + y;
#endif
	return (x < y) ? y : x;
}


//  Form /main/sort.c


static int rcmp(double x, double y, Rboolean nalast)
{
    int nax = ISNAN(x), nay = ISNAN(y);
    if (nax && nay)	return 0;
    if (nax)		return nalast ? 1 : -1;
    if (nay)		return nalast ? -1 : 1;
    if (x < y)		return -1;
    if (x > y)		return 1;
    return 0;
}
#define sort_body(TYPE_CMP, TYPE_PROT, TYPE_UNPROT)	\
    Rboolean nalast=TRUE;				\
    int i, j, h;					\
							\
    for (h = 1; h <= n / 9; h = 3 * h + 1);		\
    for (; h > 0; h /= 3)				\
	for (i = h; i < n; i++) {			\
	    v = TYPE_PROT(x[i]);			\
	    j = i;					\
	    while (j >= h && TYPE_CMP(x[j - h], v, nalast) > 0)	\
		 { x[j] = x[j - h]; j -= h; }		\
	    x[j] = v;					\
	    TYPE_UNPROT;				\
	}


void R_rsort(double *x, int n)
{
    double v;
    sort_body(rcmp,,)
}



#define psort_body						\
    Rboolean nalast=TRUE;					\
    R_xlen_t L, R, i, j;					\
								\
    for (L = lo, R = hi; L < R; ) {				\
	v = x[k];						\
	for(i = L, j = R; i <= j;) {				\
	    while (TYPE_CMP(x[i], v, nalast) < 0) i++;		\
	    while (TYPE_CMP(v, x[j], nalast) < 0) j--;		\
	    if (i <= j) { w = x[i]; x[i++] = x[j]; x[j--] = w; }\
	}							\
	if (j < k) L = i;					\
	if (k < i) R = j;					\
    }

static void rPsort2(double *x, R_xlen_t lo, R_xlen_t hi, R_xlen_t k)
{
    double v, w;
#define TYPE_CMP rcmp
    psort_body
#undef TYPE_CMP
}

void rPsort(double *x, int n, int k)
{
    rPsort2(x, 0, n-1, k);
}


//from Rmath/sign.c
double sign(double x)
{
    if (ISNAN(x))
	return x;
    return ((x > 0) ? 1 : ((x == 0)? 0 : -1));
}

double mc_C_d(const double z[], int n, const double eps[], int *iter, int scale)
{
/* NOTE:
    eps	  = c(eps1, eps2)
    iter := c(maxit, trace.lev)  as input
          = c(it, converged)     as output
*/
    int trace_lev = iter[1], it = 0;
    Rboolean converged = TRUE, do_scale = (Rboolean) scale;
    double medc; // "the" result
    static const double Large = DBL_MAX / 4.;

    if (n < 3) {
	medc = 0.; goto Finish;
    }
    /* copy data before sort()ing in place, also reflecting it -- dealing with +-Inf.
       NOTE: x[0] "empty" so we can use 1-indexing below */
    double *x  = (double *) calloc(n+1, sizeof(double));
    x[0] = 0;
    for (int i = 0; i < n; i++) {
	double zi = z[i];
	x[i+1] = - ((zi == INFINITY) ? Large :
		    (zi == -INFINITY ? -Large : zi));
    }

    R_rsort(&x[1], n); /* full sort */

    double xmed; // := median( x[1:n] ) = - median( z[0:(n-1)] ):
    if (n%2) { // n odd
	xmed = x[(n/2)+1];
    }
    else { // n even
	int ind = (n/2);
	xmed = (x[ind] + x[ind+1])/2;
    }

    double x_eps = eps[0] * (do_scale ? eps[0] + fabs(xmed) : fabs(xmed));
    if (fabs(x[1] - xmed) <= x_eps) {
	medc = -1.; goto Finish;
    } else if (fabs(x[n] - xmed) <= x_eps) {
	medc =	1.; goto Finish;
    }
    /* else : median is not at the border ------------------- */



    int i,j;
    /* center x[] wrt median --> such that then  median( x[1:n] ) == 0 */
    for (i = 1; i <= n; i++)
	x[i] -= xmed;

    if(do_scale) {
	/* MM: ==> This scaling is extremely outlier-dependent
	   --      it *kills*  equivariance when e.g. x[n] --> very large.
	   e.g., below '(eps[0] + fabs(xmed))' depends on rescaling

	   Should *NOT* be needed if everything else is *relative* instead of absolute
	   Consider replacing
	   1)  eps[0] * (eps[0] + fabs(xmed))   with  eps[0]*fabs(xmed)
	   2)        x[j] > x_eps               with     x[j] >= x_eps  (>= : for 0)
	*/

	/* Now scale to inside [-0.5, 0.5] and flip sign such that afterwards
	 *  x[1] >= x[2] >= ... >= x[n] */
	double xden = -2 * fmax2(-x[1], x[n]);
	for (i = 1; i <= n; i++)
	    x[i] /= xden;
	xmed /= xden;
	x_eps = eps[0] * (eps[0] + fabs(xmed));

    } else { // no re-scaling; still flipping signs :
	for (i = 1; i <= n; i++)
	    x[i] *= -1.;
    }

    j = 1;
    while (j <= n && x[j] >= x_eps) { /* test relative to xmed */
	/* x1[j] = x[j]; */
	j++;
    }

    i = 1;
    double *x2 = x+j-1; /* pointer -- corresponding to  x2[i] = x[j]; */
    while (j <= n && x[j] >= -x_eps) { /* test relative to xmed */
	/* x1[j] = x[j]; */
        /* x2[i] = x[j]; */
        j++;
        i++;
    }
    /* now  x1[] := {x | x_j > -eps}  also includes the median (0) */

    int h1 = j-1, /* == size of x1[] == the sum of those two sizes above */
    /* conceptually,  x2[] := {x | x_j <= eps}   (which includes the median 0) */
	h2 = i + (n-j);// == size of x2[] == maximal size of whimed() arrays


    /* work arrays for whimed_i() :  allocate *once* only !! */
    double *acand  = (double *) calloc(h2, sizeof(double)),
	   *a_srt  = (double *) calloc(h2, sizeof(double));

    int    *iw_cand= (int *)	calloc(h2, sizeof(int)),
    /* work arrays for the fast-median-of-table algorithm:
     *  currently still with  1-indexing */
	*left  = (int *) calloc((h2+1), sizeof(int)),
	*right = (int *) calloc((h2+1), sizeof(int)),
	*p     = (int *) calloc((h2+1), sizeof(int)),
	*q     = (int *) calloc((h2+1), sizeof(int));

    for (i = 1; i <= h2; i++) {
	left [i] = 1;
	right[i] = h1;
    }
    int64_t nr = ((int64_t) h1) * ((int64_t) h2), // <-- careful to *NOT* overflow
	knew = nr/2 +1;


    double trial = -2./* -Wall */;
    double *work= (double *) calloc(n, sizeof(double));
    int	   *iwt = (int *)    calloc(n, sizeof(int));
    Rboolean IsFound = FALSE;
    int nl = 0,
	neq = 0;
    /* MK:  'neq' counts the number of observations in the
     *      inside the tolerance range, i.e., where left > right + 1,
     *      since we would miss those when just using 'nl-nr'.
     *      This is to prevent index overflow in work[] later on.
     *      left might be larger than right + 1 since we are only
     *      testing with accuracy eps_trial and therefore there might
     *      be more than one observation in the `tolerance range`
     *      between < and <=.
     */
    while (!IsFound && (nr-nl+neq > n) && it < iter[0])
    {
	int64_t sum_p, sum_q;
	it++;
	j = 0;
	for (i = 1; i <= h2; i++)
	    if (left[i] <= right[i]) {
		iwt[j] = right[i] - left[i]+1;
		int k = left[i] + (iwt[j]/2);
		work[j] = h_kern(x[k], x2[i], k, i, h1+1, eps[1], do_scale);
		j++;
	    }

	trial = whimed_i(work, iwt, j, acand, a_srt, iw_cand);
	double eps_trial = eps[0] * (do_scale ? eps[0] + fabs(trial) : fabs(trial));

	j = 1;
	for (i = h2; i >= 1; i--) {
	    while (j <= h1 && h_kern(x[j],x2[i],j,i,h1+1,eps[1], do_scale) - trial > eps_trial) {
		// while (j <= h1 && h_kern(x[j],x2[i],j,i,h1+1,eps[1], do_scale) > trial) {

		j++;
	    }
/* 	    for(; j <= h1; j++) { */
/* 		register double h = h_kern(x[j],x2[i],j,i,h1+1,eps[1], do_scale); */
/* 		if(h > trial) break; */
/* 	    } */
	    p[i] = j-1;
	}
	j = h1;
	for (i = 1, sum_p=0, sum_q=0; i <= h2; i++) {
	    while (j >= 1 && trial - h_kern(x[j],x2[i],j,i,h1+1,eps[1], do_scale) > eps_trial)
		// while (j >= 1 && h_kern(x[j],x2[i],j,i,h1+1,eps[1], do_scale) < trial)
		j--;
	    q[i] = j+1;

	    sum_p += p[i];
	    sum_q += j;/* = q[i]-1 */
	}

	if(trace_lev >= 3) {
	    if (trace_lev == 3);
	    else { /* trace_lev >= 4 */
		Rboolean lrg = h2 >= 100;
		int i_m = lrg ? 95 : h2;




	    }
	}

	if (knew <= sum_p) {

	    for (i = 1, neq = 0; i <= h2; i++) {
		right[i] = p[i];
		if (left[i] > right[i]+1) neq += left[i]-right[i]-1;
	    }
	    nr = sum_p;
	}
	else { /* knew > sum_p */
	    IsFound = (knew <= sum_q); /* i.e. sum_p < knew <= sum_q */;


	    if(IsFound) {
		medc = trial;
	    } else { /*	 knew > sum_q */
	        for (i = 1; i <= h2; i++) {
		    left[i] = q[i];
		    if (left[i] > right[i]+1) neq += left[i]-right[i]-1;
		}
		nl = sum_q;
	    }
	}
	//R_CheckUserInterrupt();

    } /* end while loop */

    converged = IsFound || (nr-nl+neq <= n);
    if(!converged) {

	/* still: */
	medc = trial;
    }

    if (converged && !IsFound) { /* e.g., for  mc(1:4) : */
	j = 0;
	for (i = 1; i <= h2; i++) {
	    if (left[i] <= right[i]) {
		for (int k = left[i]; k <= right[i]; k++) {
		    work[j] = -h_kern(x[k],x2[i],k,i,h1+1,eps[1], do_scale);
		    j++;
		}
	    }
	}

	/* using rPsort(work, n,k), since we don't need work[] anymore:*/
	rPsort(work, /* n = */ j, /* k = */ knew-nl-1);
	medc = - work[knew-nl-1];
    }



Finish:
    iter[0] = it; /* to return */
    iter[1] = converged;

    return medc;

} /* end{ mc_C_d } */


/* h_kern() -- was called  calwork()  in original  rmc.c  code   and did
    if (fabs(a-b) < 2.0*eps) {
        if (ai+bi == ab) {
            return 0;
        }
        else {
            return (ai+bi < ab) ? 1 : -1 ;
        }
    }
    else {
        return (a+b)/(a-b);
    }
*/
static
double h_kern(double a, double b, int ai, int bi, int ab, double eps, Rboolean do_scale)
{
// eps := 'eps2' in R's mc()

/*     if (fabs(a-b) <= DBL_MIN) */
    /* check for zero division and positive b */
    // MK added a check '|| b > 0' ("or positive b"), but said "_and_ positive b" (r221)
    /* if (fabs(a-b) < 2.0*eps || b > 0) */
    // MM: don't see why (but it seems needed); the check should be *relative* to |a+b|
    if (b > 0 || fabs(a-b) <= eps*(do_scale ? 2. : fabs(a+b))) // '<=' since RHS maybe 0
	return sign((double)(ab - (ai+bi)));

    /* else */
    return (a+b)/(a-b);
}


/* Local variables section

 * Local variables:
 * mode: c
 * kept-old-versions: 12
 * kept-new-versions: 20
 * End:
 */
