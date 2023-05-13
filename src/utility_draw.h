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
 */
#ifndef LITERANGER_UTILITY_DRAW_H
#define LITERANGER_UTILITY_DRAW_H

/* standard library headers */
#include <algorithm>
#include <cstddef>
#include <random>

/* general literanger headers */
#include "globals.h"


namespace literanger {

/* Declarations */

/** Draw with replacement from the sequence `[0, 1, ..., M]`.
 *
 * @param[in] n The number of values from the sequence to draw.
 * @param[in] limit The largest number in the sequence.
 * @param[in] gen A pseudo-random number generator.
 * @param[out] result The samples from the sequence drawn with replacement; must
 * be initially empty.
 * @param[out] inbag_counts A container with ones in elements that are
 * in the sample, and zeroes otherwise; must be (initially) all zeroes.
 */
void draw_replace(const size_t n, const size_t limit,
                  std::mt19937_64 & gen,
                  key_vector & result, count_vector & inbag_counts);


/** Draw with replacement from the sequence `[0, 1, ..., M]` with a specified
 * probability (or weights) for each number in the sequence.
 *
 * @param[in] n The number of values from the sequence to draw.
 * @param[in] weights A vector of weights for each number in the sequence
 * (this argument specifies the largest number in the sequence).
 * @param[in] gen A pseudo-random number generator.
 * @param[out] result The samples from the sequence drawn with replacement;
 * must be initially empty.
 * @param[out] inbag_counts A container with ones in the elements that are in
 * the sample, and zeroes otherwise; must be (initially) all zeros.
 */
void draw_replace_weighted(const size_t n, const dbl_vector & weights,
                           std::mt19937_64 & gen,
                           key_vector & result, count_vector & inbag_counts);


/**
 * @param[in] n The number of values from the sequence to draw.
 * @param[in] limit The largest number in the sequence.
 * @param[in] gen A pseudo-random number generator.
 * @param[out] result The samples from the sequence drawn without; must
 * be initially empty.
 * @param[out] inbag_counts A container with ones in elements that are
 * in the sample, and zeroes otherwise; must be (initially) all zeroes.
 */
void draw_no_replace(const size_t n,
                     const size_t limit,
                     const key_vector & skipped,
                     std::mt19937_64 & gen,
                     key_vector & result,
                     count_vector & inbag_counts);


/**
 * @param[out] result The samples from the sequence drawn without; must
 * be initially empty.
 * @param[out] inbag_counts A container with ones in elements that are
 * in the sample, and zeroes otherwise; must be (initially) all zeroes.
 */
void draw_no_replace_weighted(const size_t n,
                              const dbl_vector & weights,
                              std::mt19937_64 & gen,
                              key_vector & result,
                              count_vector & inbag_counts);


/* Definitions */

inline void draw_replace(const size_t n,
                         const size_t limit,
                         std::mt19937_64 & gen,
                         key_vector & result,
                         count_vector & inbag_counts) {

    if (result.size() != 0)
        throw std::invalid_argument("Require that output vector is initially "
            "empty");
    if (inbag_counts.size() != limit)
        throw std::invalid_argument("Require that output counts is initially "
            "zero and length equal to maximum drawn value.");

    std::uniform_int_distribution<size_t> U_rng(0, limit  - 1);
    result.reserve(n);

    for (size_t j = 0; j != n; ++j) {
        const size_t draw = U_rng(gen);
        result.push_back(draw);
        ++inbag_counts[draw];
    }

}


inline void draw_replace_weighted(const size_t n,
                                  const dbl_vector & weights,
                                  std::mt19937_64 & gen,
                                  key_vector & result,
                                  count_vector & inbag_counts) {

    if (result.size() != 0)
        throw std::invalid_argument("Require that output vector is initially "
            "empty");
    if (inbag_counts.size() != weights.size())
        throw std::invalid_argument("Require that output counts is initially "
            "zero and length equal to maximum drawn value.");

    std::discrete_distribution<> wtd_rng(weights.cbegin(), weights.cend());
    result.reserve(n);

    for (size_t j = 0; j != n; ++j) {
        const size_t draw = wtd_rng(gen);
        result.push_back(draw);
        ++inbag_counts[draw];
    }

}


inline void draw_no_replace(const size_t n, const size_t limit,
                            const key_vector & skipped,
                            std::mt19937_64 & gen,
                            key_vector & result, count_vector & inbag_counts) {

    if (result.size() != 0)
        throw std::invalid_argument("Require that output vector is initially "
            "empty");
    if (inbag_counts.size() != limit)
        throw std::invalid_argument("Require that output counts is initially "
            "zero and length equal to maximum drawn value.");
    const bool n_skipped = skipped.size();

    if (n < (limit / 10)) {
      /* draw simple */
        result.reserve(n);
        std::uniform_int_distribution<size_t> U_rng(0, limit - 1 - n_skipped);

        for (size_t j = 0; j != n; ++j) {
            size_t draw;
            do {
                draw = U_rng(gen);
                if (n_skipped > 0) {
                    for (size_t skip_value : skipped) {
                        if (draw >= skip_value) ++draw;
                    }
                }
            } while (inbag_counts[draw]);
            ++inbag_counts[draw];
            result.emplace_back(draw);
        }

    } else {

      /* draw fisher-yates - in theory the same as std::shuffle, but plausible
       * that this compiles faster */
        std::uniform_real_distribution<double> U_rng(0.0, 1.0);
        result.resize(limit);
        std::iota(result.begin(), result.end(), 0);

        if (n_skipped > 0) {
            /* FIXME: could swap (to back) and erase */
            for (key_vector::const_reverse_iterator key_it = skipped.crbegin();
                     key_it != skipped.crend(); ++key_it)
                result.erase(result.begin() + *key_it);
        }
        for (size_t j = 0; j != n; ++j) {
            size_t k = j + U_rng(gen) * (limit - n_skipped - j);
            std::swap(result[j], result[k]);
            ++inbag_counts[result[j]];
        }
        result.resize(n);

    }

}


inline void draw_no_replace_weighted(const size_t n, const dbl_vector & weights,
                                     std::mt19937_64 & gen,
                                     key_vector & result,
                                     count_vector & inbag_counts) {

    if (result.size() != 0)
        throw std::invalid_argument("Require that output vector is initially "
            "empty");
    if (inbag_counts.size() != weights.size())
        throw std::invalid_argument("Require that output counts is initially "
            "zero and length equal to maximum drawn value.");

    std::discrete_distribution<> wtd_rng(weights.cbegin(), weights.cend());
    result.reserve(n);
    for (size_t j = 0; j != n; ++j) {
        size_t draw;
        /* discard if already in-bag and redraw */
        do { draw = wtd_rng(gen); } while (inbag_counts[draw]);
        result.push_back(draw);
        ++inbag_counts[draw];
    }

}


} /* namespace literanger */


#endif /* LITERANGER_UTILITY_DRAW_H */

