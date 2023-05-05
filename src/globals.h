/* This file was adapted from the C++ core of the "ranger" package for R
 * Statistical Software.
 *
 * Adaptation was authored by Stephen Wade. The same license terms as the
 * original c++ core of the ranger package apply to the adaptation.
 *
 * License statement for c++ core of ranger:
 *
 * Copyright (c) [2014-2018] [Marvin N. Wright]
 *
 * This software may be modified and distributed under the terms of the MIT
 * license.
 *
 * Please note that the C++ core of ranger is distributed under MIT license and
 * the R package "ranger" under GPL3 license.
 */
#ifndef LITERANGER_GLOBALS_H
#define LITERANGER_GLOBALS_H

/* standard library headers */
#include <bitset>
#include <cstddef>
#include <limits>
#include <vector>


namespace literanger {

using count_vector = std::vector<size_t>;
using key_vector = std::vector<size_t>;
using dbl_vector = std::vector<double>;
using ull_bitenc = std::bitset<std::numeric_limits<size_t>::digits>;

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

const size_t DEFAULT_N_THREAD = 0;

const size_t DEFAULT_MIN_SPLIT_N_SAMPLE_CLASSIFICATION = 2;
const size_t DEFAULT_MIN_LEAF_N_SAMPLE_CLASSIFICATION = 1;
const size_t DEFAULT_MIN_SPLIT_N_SAMPLE_REGRESSION = 5;
const size_t DEFAULT_MIN_LEAF_N_SAMPLE_REGRESSION = 1;

const double STATUS_INTERVAL = 30.0;

// Threshold for q value split method switch
const double Q_THRESHOLD = 0.02;


} /* namespace literanger */


#endif /* LITERANGER_GLOBALS_H */

