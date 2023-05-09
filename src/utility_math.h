/* This file was adapted from the C++ core of the "ranger" package for R
 * Statistical Software.
 *
 * Adaptation was authored by Stephen Wade. The same license terms as the
 * original c++ core of the ranger package apply to the adaptation.
 *
 * License statement for C++ core of ranger:
 *
 * Copyright (c) [2014-2018] [Marvin N. Wright]
 *
 * This software may be modified and distributed under the terms of the MIT
 * license.
 *
 * Please note that the C++ core of ranger is distributed under MIT license and
 * the R package "ranger" under GPL3 license.
 */
#ifndef LITERANGER_UTILITY_MATH_H
#define LITERANGER_UTILITY_MATH_H

/* standard library headers */
#include <algorithm>
#include <cmath>
#include <limits>
#include <mutex>

/* general literanger headers */
#include "utility.h"


namespace literanger {

inline double beta_log_likelihood(double y, double mu, double nu) {

    static std::mutex lgamma_mutex;
    { /* Avoid 0 and 1 */
        const double eps = std::numeric_limits<double>::epsilon();
        y = std::max(eps, std::min(1 - eps, y));
        mu = std::max(eps, std::min(1 - eps, mu));
        nu = std::max(eps, std::min(1 - eps, nu));
    }
    double lgamma_result;
  /* FIXME: For now, we lock as lgamma isn't thread-safe */
    { std::unique_lock<std::mutex> lock(lgamma_mutex);
        lgamma_result = std::lgamma(nu) - std::lgamma(mu * nu) -
            std::lgamma((1 - mu) * nu); }

    return lgamma_result + (mu * nu - 1) * std::log(y) +
        ((1 - mu) * nu - 1) * std::log(1 - y);

}


inline double dZ(const double x) {
    return std::exp(-0.5 * x * x) / std::sqrt(2 * M_PI);
}


inline double pZ(const double x) {
    return 0.5 * (1 + std::erf(x / std::sqrt(2.0)));
}


inline double maxstat_p_value_Lausen92(const double b, const double lower) {
  /* assume upper = 1 - lower */
    if (b < 1) return 1;
    const double ln_ratio = 2 * std::log((1 - lower) / lower);
    const double db = dZ(b);
    const double p = (4 + (std::pow(b, 2) - 1) * ln_ratio) * (db / b);
  /* note: only valid for large b and may be g.t. one for small enough b,
   * therefore cap at one. */
    return std::max(0.0, std::min(1.0, p));
}


inline double maxstat_p_value_Lausen94(const double b, const size_t N,
                                       const std::vector<size_t> & m,
                                       const size_t n) {

    if (n < 2) return 2 * (1 - pZ(b));

    double D = 0;
    const double bsq = std::pow(b, 2);
    for (size_t j = 0; j != n - 1; ++j) {
        const double m_j = m[j];
        const double m_jp1 = m[j + 1];
        const double t = std::sqrt(
            1.0 - (m_j / m_jp1) * (N - m_jp1) / (N - m_j)
        );
        D += std::exp(-0.5 * bsq) * (t - (bsq / 4 - 1) * std::pow(t, 3) / 6);
    }

    return std::max(0.0, std::min(1.0, 2 * (1 - pZ(b)) + (D / M_PI)));

}


/** Compute adjusted p-values using Benjamini/Hochberg method */
inline std::vector<double> adjust_pvalues(const std::vector<double> & unadjusted) {

    const size_t n_pvalue = unadjusted.size();
    if (n_pvalue < 2) return unadjusted; /* short circuit */

  /* Order of p-values (decreasing) */
    const std::vector<size_t> index = order<true>(unadjusted);

    std::vector<double> adjusted(n_pvalue, 0);
    adjusted[index[0]] = unadjusted[index[0]];

    for (size_t j = 1; j != n_pvalue; ++j) {
        const size_t order_j = index[j];
        const size_t order_jm1 = index[j - 1];
        adjusted[order_j] = std::min(
            adjusted[order_jm1],
            n_pvalue / (double)(n_pvalue - j) * unadjusted[order_j]
        );
    }

    return adjusted;

}


} /* namespace literanger */


#endif /* LITERANGER_UTILITY_MATH_H */

