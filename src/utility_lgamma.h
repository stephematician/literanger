/* This file is part of the C++ core of 'literanger'.
 *
 * literanger's C++ core was adapted from the C++ core of the 'ranger' package
 * for R Statistical Software <https://www.r-project.org>. The ranger C++ core
 * is Copyright (c) [2014-2018] [Marvin N. Wright] and distributed with MIT
 * license. literanger's C++ core is distributed with the same license, terms,
 * and permissions as ranger's C++ core.
 *
 * Copyright [2023] [Stephen Wade]
 *
 * This software may be modified and distributed under the terms of the MIT
 * license. You should have received a copy of the MIT license along with
 * literanger. If not, see <https://opensource.org/license/mit/>.
 *
 *
 * non-negative log gamma adapted from glibc
 * `sysdeps/ieee754/dbl-64/e_lgamma_r.c` with the following license:
 *
 * ====================================================
 * Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.
 *
 * Developed at SunPro, a Sun Microsystems, Inc. business.
 * Permission to use, copy, modify, and distribute this
 * software is freely granted, provided that this notice
 * is preserved.
 * ====================================================
 */
#ifndef LITERANGER_UTILITY_LGAMMA_H
#define LITERANGER_UTILITY_LGAMMA_H

#include <cmath>
#include <cstdint>


namespace literanger {

union double_parts {
    double x;
    std::int32_t words[2];
};


/* Documentation from `__ieee754_lgamma_r(x, signgamp)`
 *
 * Reentrant version of the logarithm of the Gamma function
 * with user provide pointer for the sign of Gamma(x).
 *
 * Method:
 *   1. Argument Reduction for 0 < x <= 8
 *	Since gamma(1+s)=s*gamma(s), for x in [0,8], we may
 *	reduce x to a number in [1.5,2.5] by
 *		lgamma(1+s) = log(s) + lgamma(s)
 *	for example,
 *		lgamma(7.3) = log(6.3) + lgamma(6.3)
 *			    = log(6.3*5.3) + lgamma(5.3)
 *			    = log(6.3*5.3*4.3*3.3*2.3) + lgamma(2.3)
 *   2. Polynomial approximation of lgamma around its
 *	minimum ymin=1.461632144968362245 to maintain monotonicity.
 *	On [ymin-0.23, ymin+0.27] (i.e., [1.23164,1.73163]), use
 *		Let z = x-ymin;
 *		lgamma(x) = -1.214862905358496078218 + z^2*poly(z)
 *	where
 *		poly(z) is a 14 degree polynomial.
 *   2. Rational approximation in the primary interval [2,3]
 *	We use the following approximation:
 *		s = x-2.0;
 *		lgamma(x) = 0.5*s + s*P(s)/Q(s)
 *	with accuracy
 *		|P/Q - (lgamma(x)-0.5s)| < 2**-61.71
 *	Our algorithms are based on the following observation
 *
 *                             zeta(2)-1    2    zeta(3)-1    3
 * lgamma(2+s) = s*(1-Euler) + --------- * s  -  --------- * s  + ...
 *                                 2                 3
 *
 *	where Euler = 0.5771... is the Euler constant, which is very
 *	close to 0.5.
 *
 *   3. For x>=8, we have
 *	lgamma(x)~(x-0.5)log(x)-x+0.5*log(2pi)+1/(12x)-1/(360x**3)+....
 *	(better formula:
 *	   lgamma(x)~(x-0.5)*(log(x)-1)-.5*(log(2pi)-1) + ...)
 *	Let z = 1/x, then we approximation
 *		f(z) = lgamma(x) - (x-0.5)(log(x)-1)
 *	by
 *				    3       5             11
 *		w = w0 + w1*z + w2*z  + w3*z  + ... + w6*z
 *	where
 *		|w - f(z)| < 2**-58.74
 *
 *   4. For negative x, since (G is gamma function)
 *		-x*G(-x)*G(x) = pi/sin(pi*x),
 *	we have
 *		G(x) = pi/(sin(pi*x)*(-x)*G(-x))
 *	since G(-x) is positive, sign(G(x)) = sign(sin(pi*x)) for x<0
 *	Hence, for x<0, signgam = sign(sin(pi*x)) and
 *		lgamma(x) = log(|Gamma(x)|)
 *			  = log(pi/(|x*sin(pi*x)|)) - lgamma(-x);
 *	Note: one should avoid compute pi*(-x) directly in the
 *	      computation of sin(pi*(-x)).
 *
 *   5. Special Cases
 *		lgamma(2+s) ~ s*(1-Euler) for tiny s
 *		lgamma(1)=lgamma(2)=0
 *		lgamma(x) ~ -log(x) for tiny x
 *		lgamma(0) = lgamma(inf) = inf
 *		lgamma(-integer) = +-inf
 */
inline double lgamma_nn(const double x) {

  /* Only support non-negative numbers in this adaptation. Do not store sign of
   * gamma function; simply return log |Gamma(x)| */

    if (std::isinf(x)) return INFINITY;
    if (std::isnan(x)) return NAN;

    const double_parts x_parts = { x };

    const std::int32_t & lx = x_parts.words[0];
    const std::int32_t & hx = x_parts.words[1];

    /* where are we? */
    if (hx < 0x3b900000) return -std::log(x);

  /* x == 1 and x == 2 */
	if(((hx - 0x3ff00000) | lx) == 0 || ((hx - 0x40000000) | lx) == 0) return 0;

  /* x < 2.0 */
    if(hx < 0x40000000) {
      /* [x_min,f_min] = local minimum of log gamma; unclear to me what `tt` is */
        constexpr double x_min =  1.46163214496836224576e+00,
                         f_min = -1.21486290535849611461e-01,
                       /* tt = -(tail of f_min) */
                         tt = -3.63867699703950536541e-18;
        double result, y;
        int i;
	    if(hx <= 0x3feccccc) {
          /* lgamma(x) = lgamma(x + 1) - log(x) */
		    result = -std::log(x);
		    if (hx >= 0x3FE76944) { y = 1.0 - x; i = 0; }
		    else if(hx >= 0x3FCDA661) { y = x - (x_min - 1.0); i = 1; }
		    else { y = x; i = 2; }
	    } else {
		    result = 0.0;
		    if (hx >= 0x3FFBB4C3) {y = 2.0 - x; i = 0; } /* [1.7316,2] */
		    else if(hx >= 0x3FF3B4C4) { y = x - x_min; i = 1; } /* [1.23,1.73] */
		    else { y = x - 1.0; i = 2; }
	    }
	    switch(i) {
	    case 0: {
            constexpr double a [12] = {
                7.72156649015328655494e-02, 3.22467033424113591611e-01,
                6.73523010531292681824e-02, 2.05808084325167332806e-02,
                7.38555086081402883957e-03, 2.89051383673415629091e-03,
                1.19270763183362067845e-03, 5.10069792153511336608e-04,
                2.20862790713908385557e-04, 1.08011567247583939954e-04,
                2.52144565451257326939e-05, 4.48640949618915160150e-05
            };
            const double z = y * y;
		    const double p1 = a[0] + z * (a[2] + z * (
                                  a[4] + z * (a[6] + z * (a[8] + z * a[10]))));
		    const double p2 = z * (a[1] + z * (a[3] + z * (
                                  a[5] + z * (a[7] + z * (a[9] + z * a[11])))));
		    const double p  = y * p1 + p2;
		    result += (p - 0.5 * y);
        } break;
	    case 1: {
            constexpr double a [15] {
                4.83836122723810047042e-01, -1.47587722994593911752e-01,
                6.46249402391333854778e-02, -3.27885410759859649565e-02,
                1.79706750811820387126e-02, -1.03142241298341437450e-02,
                6.10053870246291332635e-03, -3.68452016781138256760e-03,
                2.25964780900612472250e-03, -1.40346469989232843813e-03,
                8.81081882437654011382e-04, -5.38595305356740546715e-04,
                3.15632070903625950361e-04, -3.12754168375120860518e-04,
                3.35529192635519073543e-04
            };
		    const double z = std::pow(y, 3); // y * y;
		    const double p1 = a[0] + z * (
                                  a[3] + z * (a[6] + z * (a[9] + z * a[12])));
		    const double p2 = a[1] + z * (
                                  a[4] + z * (a[7] + z * (a[10] + z * a[13])));
		    const double p3 = a[2] + z * (
                                  a[5] + z * (a[8] + z * (a[11] + z * a[14])));
		    const double p  = (y * y) * p1 - (tt - z * (p2 + y * p3));
		    result += (f_min + p);
        } break;
	    case 2:
            constexpr double a [6] = {
               -7.72156649015328655494e-02, 6.32827064025093366517e-01,
                1.45492250137234768737e+00, 9.77717527963372745603e-01,
                2.28963728064692451092e-01, 1.33810918536787660377e-02
            };
            constexpr double b [6] = {
                1.00000000000000000000e+00, 2.45597793713041134822e+00,
                2.12848976379893395361e+00, 7.69285150456672783825e-01,
                1.04222645593369134254e-01, 3.21709242282423911810e-03
            };
		    const double p1 = y * (a[0] + y * (a[1] + y * (
                                  a[2] + y * (a[3] + y * (a[4] + y * a[5])))));
		    const double p2 = b[0] + y * (b[1] + y * (
                                  b[2] + y * (b[3] + y * (b[4] + y * b[5]))));
		    result -= (0.5 * y - p1 / p2);
	    }
        return result;

  /* 2.0 < x < 8.0 */
	} else if(hx < 0x40200000) {
        constexpr double a [7] = {
           -7.72156649015328655494e-02, 2.14982415960608852501e-01,
            3.25778796408930981787e-01, 1.46350472652464452805e-01,
            2.66422703033638609560e-02, 1.84028451407337715652e-03,
            3.19475326584100867617e-05
        };
       constexpr double b [7] = {
            1.00000000000000000000e+00, 1.39200533467621045958e+00,
            7.21935547567138069525e-01, 1.71933865632803078993e-01,
            1.86459191715652901344e-02, 7.77942496381893596434e-04,
            7.32668430744625636189e-06
        };
	    const int i = (int)x;
	    const double y = x - (double)i;
	    const double p = y * (a[0] + y * (a[1] + y * (a[2] + y * (
                                 a[3] + y * (a[4] + y * (a[5] + y * a[6]))))));
	    const double q = b[0] + y * (b[1] + y * (b[2] + y * (
                             b[3] + y * (b[4] + y * (b[5] + y * b[6])))));
	    double result = 0.5 * y + p / q;
	    double z = 1.0; /* lgamma(1 + s) = log(s) + lgamma(s) */
	    switch(i) {
	    case 7: z *= (y + 6.0);	/* FALLTHRU */
	    case 6: z *= (y + 5.0);	/* FALLTHRU */
	    case 5: z *= (y + 4.0);	/* FALLTHRU */
	    case 4: z *= (y + 3.0);	/* FALLTHRU */
	    case 3: z *= (y + 2.0);	/* FALLTHRU */
		    result += std::log(z); break;
	    }
        return result;

  /* 8.0 <= x < 2**58 */
	} else if (hx < 0x43900000) {
        constexpr double a [7] = {
            4.18938533204672725052e-01, 8.33333333333329678849e-02,
           -2.77777777728775536470e-03, 7.93650558643019558500e-04,
           -5.95187557450339963135e-04, 8.36339918996282139126e-04,
           -1.63092934096575273989e-03
        };
	    const double log_x = std::log(x);
	    const double z = 1.0 / x;
	    const double y = z * z;
	    const double w = a[0] + z * (a[1] + y * (a[2] + y * (
                                a[3] + y * (a[4] + y * (a[5] + y * a[6])))));
	    return (x - 0.5) * (log_x - 1.0) + w;
  /* 2**58 <= x <= inf */
	} else
        return x * (std::log(x) - 1.0);

}


} /* namespace literanger */


#endif /* LITERANGER_UTILITY_LGAMMA_H */

