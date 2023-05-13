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
#ifndef LITERANGER_DATA_DEFN_H
#define LITERANGER_DATA_DEFN_H

/* class declaration */
#include "Data.decl.h"

/* standard library headers */
#include <algorithm>
#include <iterator>
#include <numeric>
#include <random>
#include <stdexcept>


namespace literanger {

inline Data::Data(const size_t n_row, const size_t n_col) :
    n_row(n_row), n_col(n_col) { }


inline Data::~Data() { }


inline size_t Data::get_n_col() const { return n_col; }


inline size_t Data::get_n_row() const { return n_row; }


inline void Data::get_all_values(dbl_vector & all_values,
                                 const key_vector & sample_keys,
                                 const size_t predictor_key,
                                 const size_t start, const size_t end,
                                 const bool permute) const {

    if (start > end)
        throw std::invalid_argument("Start of interval must not be past end.");

    all_values.reserve(end - start);

    for (size_t j = start; j != end; ++j)
        all_values.emplace_back(get_x(sample_keys[j], predictor_key, permute));

    std::sort(all_values.begin(), all_values.end());
    all_values.erase(std::unique(all_values.begin(), all_values.end()),
                     all_values.end());

}


inline void Data::get_minmax_values(double & min, double & max,
                                    const key_vector & sample_keys,
                                    const size_t predictor_key,
                                    const size_t start, const size_t end,
                                    const bool permute) const {

    if (start > end)
        throw std::invalid_argument("Start of interval must not be past end.");

    if (sample_keys.size() > 0)
        max = min = get_x(sample_keys[start], predictor_key, permute);

    for (size_t j = start; j != end; ++j) {
        double value = get_x(sample_keys[j], predictor_key, permute);
        min = std::min(value, min);
        max = std::max(value, max);
    }

}


inline void Data::new_predictor_index() const {

    predictor_index.assign(n_col * n_row, 0);
    unique_predictor_values.clear();
    unique_predictor_values.reserve(n_col);
    max_n_unique_value = 0;

    for (size_t column = 0; column != n_col; ++column) {

      /* Get the unique values and update the index for this predictor */
        dbl_vector unique_values(n_row);
        for (size_t row = 0; row != n_row; ++row)
            unique_values[row] = get_x(row, column);

        std::sort(unique_values.begin(), unique_values.end());
        unique_values.erase(
            std::unique(unique_values.begin(), unique_values.end()),
            unique_values.end()
        );

        for (size_t row = 0; row < n_row; ++row)
            predictor_index[column * n_row + row] = std::distance(
                unique_values.begin(),
                std::lower_bound(
                    unique_values.begin(), unique_values.end(),
                    get_x(row, column)
                )
            );

      /* Save unique values for this predictor */
        unique_predictor_values.push_back(unique_values);
        max_n_unique_value = std::max(unique_values.size(), max_n_unique_value);

    }

}


inline void Data::finalise_predictor_index() const {
    predictor_index.clear();
    predictor_index.shrink_to_fit();
    unique_predictor_values.clear();
    unique_predictor_values.shrink_to_fit();
    max_n_unique_value = 0ull;
}


inline bool Data::has_predictor_index() const {
    return max_n_unique_value != 0;
}


inline size_t Data::get_index(const size_t sample_key,
                              const size_t predictor_key,
                              const bool permute) const {
    const size_t row = as_row_offset(sample_key, permute);
    if (predictor_key >= n_col)
        throw std::invalid_argument("Predictor key must be less than number of "
            "columns.");
    return predictor_index[predictor_key * n_row + row];
}


inline double Data::get_unique_predictor_value(const size_t predictor_key,
                                               const size_t offset) const {
    if (predictor_key >= n_col)
        throw std::invalid_argument("Predictor key must be less than number of "
            "columns.");
    return unique_predictor_values[predictor_key][offset];
}


inline size_t Data::get_n_unique_predictor_value(
    const size_t predictor_key
) const {
    if (predictor_key >= n_col)
        throw std::invalid_argument("Predictor key must be less than number of "
            "columns.");
    return unique_predictor_values[predictor_key].size();
}


inline size_t Data::get_max_n_unique_value() const {
    return std::max((size_t)3ull, max_n_unique_value);
    /* NOTE: unsure why lower bound of three */
}


inline dbl_vector Data::get_response_values() const {

    const size_t n_sample = get_n_row();
    dbl_vector result { };

    for (size_t sample_key = 0; sample_key != n_sample; ++sample_key) {
        const double value = get_y(sample_key, 0);
        const dbl_vector::const_iterator value_it = std::find(
            result.cbegin(), result.cend(), value
        );
        if (value_it == result.cend()) result.push_back(value);
    }

    return result;

}


inline void Data::new_response_index(const dbl_vector & response_values) const {

    response_index.clear();
    response_index.reserve(get_n_row());

    for (size_t j = 0; j != get_n_row(); ++j) {
        const size_t value_key = std::distance(
            response_values.cbegin(),
            std::find(response_values.cbegin(),
                      response_values.cend(),
                      get_y(j, 0))
        );
        if (value_key == response_values.size())
            throw std::invalid_argument("Response values does not contain all "
                "values observe in data.");
        response_index.emplace_back(value_key);
    }

}


inline void Data::finalise_response_index() const {
    response_index.clear();
    response_index.shrink_to_fit();
}


inline const key_vector & Data::get_response_index() const {
    return response_index;
}


inline void Data::new_sample_keys_by_response() const {

    sample_keys_by_response.assign(response_index.size(), key_vector());

    for (size_t j = 0; j != get_n_row(); ++j) {
        const size_t value_key = response_index[j];
        sample_keys_by_response[value_key].emplace_back(j);
    }

}


inline void Data::finalise_sample_keys_by_response() const {
    sample_keys_by_response.clear();
    sample_keys_by_response.shrink_to_fit();
}


inline const std::vector<key_vector> &
Data::get_sample_keys_by_response() const {
    return sample_keys_by_response;
}


inline void Data::new_permutation(const size_t seed) const {

    std::mt19937_64 gen;

    if (seed == 0) {
        std::random_device random_device;
        gen.seed(random_device());
    } else {
        gen.seed(seed);
    }

    permuted_sample_keys.assign(n_row, 0);
    std::iota(permuted_sample_keys.begin(), permuted_sample_keys.end(), 0);
    std::shuffle(permuted_sample_keys.begin(), permuted_sample_keys.end(), gen);

}


inline void Data::finalise_permutation() const {
    permuted_sample_keys.clear();
    permuted_sample_keys.shrink_to_fit();
}


inline size_t Data::as_row_offset(const size_t sample_key,
                                  const bool permute) const {
    return !permute ? sample_key : permuted_sample_keys[sample_key];
}


} /* namespace literanger */


#endif /* LITERANGER_DATA_DEFN_H */

