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
#ifndef LITERANGER_UTILITY_H
#define LITERANGER_UTILITY_H

/* standard library headers */
#include <algorithm>
#include <cstddef>
#include <iterator>
#include <memory>
#include <random>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>


namespace literanger {

/** Toggle-able printer */
struct toggle_print {

    /** R printf via forwarding
     * @param fmt Format specifier string e.g. "hello %f".
     * @param args Formatted arguments.
     * @tparam Types of formatted argument. */
    template <typename... ArgsT>
    void operator()(const char *fmt, ArgsT &&... args);

    /** Toggle (indicator) whether to print or not. */
    bool verbose;

    /** Function reference for a printf-like operation.
     * @param fmt Formatting string for a printf-like operation.
     * @param ... Additional arguments accepted by printf. */
    void (& print_fn)(const char * fmt, ...);

};


/** Split (or partition) a range of consecutive integers into equal parts.
 *
 * @param[out] result A vector of starting values of each partition and the
 *   value one past the end.
 * @param[in] start The first value in the entire range.
 * @param[in] end The last value in the entire range.
 * @param[in] n The number of partitions to split the range into (equally).
 */
void equal_split(std::vector<size_t> & result,
                 const size_t start, const size_t end,
                 const size_t n);


/**
 * Return the value with the greatest count if unique, or a randomly-selected
 * value amongst those with the greatest count
 *
 * @param[in] counts A map between values and the frequency of the value
 * @param[in] gen A random number generation engine
 * @param[in] order_invariant Indicator of whether ties are broken using a
 *   random draw that is _invariant_ to the order the values appear.
 */
template <typename ValueT, typename CountT>
ValueT most_frequent_value(const std::unordered_map<ValueT,CountT> & counts,
                           std::mt19937_64 & gen,
                           const bool order_invariant = true);


/** Get order of elements (decreasing) */
template <bool decreasing, typename ContainerT,
          typename std::enable_if<decreasing, std::nullptr_t>::type = nullptr>
std::vector<size_t> order(const ContainerT & x);


/** Get order of elements (increasing) */
template <bool decreasing, typename ContainerT,
          typename std::enable_if<!decreasing, std::nullptr_t>::type = nullptr>
std::vector<size_t> order(const ContainerT & x);


/** Sample ranks starting from 1. Ties are given the average rank.
 * @param x Values to rank
 * @returns Ranks of input values
 */
template <typename ContainerT>
std::vector<double> rank(const ContainerT & x);


/** Make a (shared resource) container for indicators of ordered predictors.
 * @param[in] predictor_names The names of the predictor variables in the order
 * they (will) appear in the data.
 * @param[in] names_of_unordered The names of the predictor variables
 * that are 'unordered' (and hence will be split by partition).
 * @tparam PtrT The type that manages a shared resource, defaults to a shared
 * pointer. */
template <template <typename...> class PtrT = std::shared_ptr>
PtrT<std::vector<bool>> make_is_ordered(
    const std::vector<std::string> & predictor_names,
    const std::vector<std::string> & names_of_unordered
);


/** Make a container of keys for the predictors that are always
 * candidates.
 * @param[in] predictor_names The names of the predictor variables in the order
 * they (will) appear in the data.
 * @param[in] names_of_always_draw The names of predictor that will always
 * be candidates for splitting.
 * @param[in] n_try The number of candidate predictors for each split.
 * @tparam PtrT The type that manages a shared resource, defaults to a shared
 * pointer. */
template <template <typename...> class PtrT = std::shared_ptr>
PtrT<std::vector<size_t>> make_draw_always_predictor_keys(
    const std::vector<std::string> & predictor_names,
    const std::vector<std::string> & names_of_always_draw,
    const size_t n_try
);


/** Get the column offset for a predictor by name.
 * @param[in] predictor_names The names of the predictor variables in the order
 * they (will) appear in the data.
 * @param[in] name The name of the predictor. */
size_t get_predictor_key(const std::vector<std::string> & predictor_names,
                         const std::string & name);


/* Definitions */

template <typename... ArgsT>
void toggle_print::operator()(const char *fmt, ArgsT &&... args) {
    if (verbose) print_fn(fmt, std::forward<ArgsT>(args)...);
}


inline void equal_split(std::vector<size_t> & result,
                        const size_t start, const size_t end,
                        const size_t n) {

    if (n == 0) throw std::domain_error("Cannot split into zero parts.");

    result.clear();
    result.reserve(n + 1);

    const size_t n_ = std::min(n, (end + 1) - start);
    const size_t dist = (end + 1) - start;
    const size_t chunk = dist / n_;
    size_t remainder = dist % n_;
    size_t next = start;

    for (size_t i = 0; i != n_; ++i) {
        result.emplace_back(next);
        next += chunk + (remainder ? (--remainder, 1) : 0);
    }
    result.emplace_back(end + 1);

}


template <typename ValueT, typename CountT>
ValueT most_frequent_value(const std::unordered_map<ValueT,CountT> & counts,
                           std::mt19937_64 & gen,
                           const bool order_invariant) {

    if (counts.empty())
        throw std::invalid_argument("Cannot find most frequent value for empty "
            "map.");

    std::vector<ValueT> major_values;
    major_values.reserve(counts.size());

  /* Get the value(s) with maximum count. */
    CountT max_count = 0;
    for (const auto & kv : counts) max_count = std::max(kv.second, max_count);

    for (const auto & kv : counts) {
        const CountT count = kv.second;
        const ValueT value = kv.first;
        if (count == max_count) major_values.emplace_back(value);
    }

  /* Short circuit if unique maximum. */
    if (major_values.size() == 1) return major_values[0];

    if (major_values.size() <= 1)
        throw std::runtime_error("Did not expect empty most frequent values.");

  /* random draw from the set of values with maximum count */
    std::uniform_int_distribution<size_t> U_rng(0, major_values.size() - 1);
    if (order_invariant) std::sort(major_values.begin(), major_values.end());
    return major_values[U_rng(gen)];

}


template <bool decreasing, typename ContainerT,
          typename std::enable_if<decreasing, std::nullptr_t>::type>
std::vector<size_t> order(const ContainerT & x) {
    std::vector<size_t> index(x.size());
    std::iota(index.begin(), index.end(), 0);
    std::sort(std::begin(index), std::end(index),
              [&](size_t lhs, size_t rhs){ return x[lhs] > x[rhs]; });
    return index;
}


template <bool decreasing, typename ContainerT,
          typename std::enable_if<!decreasing, std::nullptr_t>::type>
std::vector<size_t> order(const ContainerT & x) {
    std::vector<size_t> index(x.size());
    std::iota(index.begin(), index.end(), 0);
    std::sort(std::begin(index), std::end(index),
              [&](size_t lhs, size_t rhs){ return x[lhs] < x[rhs]; });
    return index;
}


template <typename ContainerT>
std::vector<double> rank(const ContainerT & x) {

    const size_t n_value = x.size();
    const std::vector<size_t> index = order<false>(x);

    std::vector<double> result(n_value);

    for (size_t end_j = 0; end_j != n_value;) {
        const size_t start_j = end_j;
      /* Get index of next non-equivalent element */
        while (end_j != n_value && x[index[start_j]] == x[index[end_j]])
            ++end_j;
      /* Assign (same) rank to all equivalent values */
        for (size_t k = start_j; k != end_j; ++k)
            result[index[k]] = (start_j + end_j - 1) / 2.0 + 1;
    }

    return result;

}


template <template <typename...> class PtrT>
PtrT<std::vector<bool>> make_is_ordered(
    const std::vector<std::string> & predictor_names,
    const std::vector<std::string> & names_of_unordered
) {
    using key_vector = std::vector<bool>;
    const size_t n_predictor = predictor_names.size();
    PtrT<key_vector> result(new key_vector(n_predictor, true));

    for (auto & name : names_of_unordered) {
        const size_t key = get_predictor_key(predictor_names, name);
        (*result)[key] = false;
    }

    return result;
}


template <template <typename...> class PtrT>
PtrT<std::vector<size_t>> make_draw_always_predictor_keys(
    const std::vector<std::string> & predictor_names,
    const std::vector<std::string> & names_of_always_draw,
    const size_t n_try
) {
    using key_vector = std::vector<size_t>;
    const size_t n_predictor = predictor_names.size();
    PtrT<key_vector> result( new key_vector() );

    if (names_of_always_draw.empty()) return result;
    result->reserve(n_predictor);

    for (auto & name: names_of_always_draw) {
        const size_t key = get_predictor_key(predictor_names, name);
        result->emplace_back(key);
    }

    if (result->size() + n_try > n_predictor)
        throw std::domain_error("Number of predictors to always consider "
            "splitting plus 'n_try' cannot be larger than total number of "
            "predictors (columns)");

    return result;
}


inline size_t get_predictor_key(
    const std::vector<std::string> & predictor_names,
    const std::string & name
) {
    auto it = std::find(predictor_names.cbegin(), predictor_names.cend(), name);
    if (it == std::end(predictor_names))
        throw std::invalid_argument("predictor `" + name + "` not found.");

    return (std::distance(predictor_names.cbegin(), it));
}


} /* namespace literanger */


#endif /* LITERANGER_UTILITY_H */

