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
#ifndef LITERANGER_TREE_DEFN_H
#define LITERANGER_TREE_DEFN_H

/* class declaration */
#include "Tree.decl.h"

/* standard library headers */
#include <algorithm>
#include <cmath>
#include <limits>
#include <random>
#include <stdexcept>
#include <utility>
#include <vector>

/* general literanger headers */
#include "utility_draw.h"
#include "utility_math.h"
/* required literanger class definitions */
#include "Data.defn.h"
#include "TreeBase.defn.h"


namespace literanger {

template <typename ImplT>
template <typename... ArgsT>
Tree<ImplT>::Tree(ArgsT &&... args) :
    TreeBase(std::forward<ArgsT>(args)...)
{ }


template <typename ImplT>
template <PredictionType prediction_type, typename result_type>
void Tree<ImplT>::predict(const std::shared_ptr<const Data> data,
                          const size_t sample_key,
                          result_type & result) {

    ImplT & tree_impl = *static_cast<ImplT *>(this);

    size_t node_key = 0, depth = 0;

    while (!max_depth || depth != max_depth) {

      /* Terminal node */
        if (left_children[node_key] == 0 && right_children[node_key] == 0)
            break;

        const size_t split_key = split_keys[node_key];
        const double value = data->get_x(sample_key, split_key);

        if ((*is_ordered)[split_key]) {
            if (value <= split_values[node_key]) {
                node_key = left_children[node_key];
            } else {
                node_key = right_children[node_key];
            }
        } else {
          /* NOTE: probably unsafe */
            const ull_bitenc split_enc = *((size_t *)(&split_values[node_key]));
            if (!split_enc.test(std::floor(value) - 1)) {
                node_key = left_children[node_key];
            } else {
                node_key = right_children[node_key];
            }
        }
        ++depth;

    }

    if (max_depth && depth == max_depth &&
            !(left_children[node_key] == 0 && right_children[node_key] == 0))
        throw std::runtime_error("Prediction failure tree does not obey "
            "maximum depth constraint.");

  /* Get a prediction from this tree */
    tree_impl.template predict_from_inbag<prediction_type>(node_key, result);

}


template <typename ImplT>
bool Tree<ImplT>::push_best_split(
    const size_t node_key,
    const std::shared_ptr<const Data> data,
    const key_vector & sample_keys,
    const key_vector & split_candidate_keys
) {

    double best_decrease = -INFINITY, best_value = 0;
    size_t best_split_key = 0;

  /* Prepares aggregates such as the count by response-value */
    new_node_aggregates(node_key, data, sample_keys);

    switch (split_rule) {
    case EXTRATREES: {
        for (const size_t split_key : split_candidate_keys) {
            if ((*is_ordered)[split_key]) {
                best_decrease_by_value_extratrees(
                    split_key, node_key, data, sample_keys,
                    best_decrease, best_split_key, best_value
                );
            } else {
                best_decrease_by_value_extratrees_unordered(
                    split_key, node_key, data, sample_keys,
                    best_decrease, best_split_key, best_value
                );
            }
        }
    } break;
    case LOGRANK: case BETA: case HELLINGER: {
        for (const size_t split_key : split_candidate_keys) {

            if ((*is_ordered)[split_key]) {
                const size_t n_sample_node = get_n_sample_node(node_key);
                const double q = n_sample_node / (
                    data->has_predictor_index() ?
                        (double)data->get_n_unique_predictor_value(split_key) :
                        INFINITY
                );
                if (!data->has_predictor_index() || q < Q_THRESHOLD) {
                    best_decrease_by_value_smallq(
                        split_key, node_key, data, sample_keys,
                        best_decrease, best_split_key, best_value
                    );
                } else {
                    best_decrease_by_value_largeq(
                        split_key, node_key, data, sample_keys,
                        best_decrease, best_split_key, best_value
                    );
                }
            } else {
                best_decrease_by_value_unordered(
                    split_key, node_key, data, sample_keys,
                    best_decrease, best_split_key, best_value
                );
            }
        }
    } break;
    case MAXSTAT: {
        double best_statistic = -INFINITY;
        std::vector<double> p_values, keys;
        p_values.reserve(split_candidate_keys.size());
        keys.reserve(split_candidate_keys.size());

        for (const size_t split_key : split_candidate_keys) {
            if (!(*is_ordered)[split_key])
                throw std::invalid_argument("Maximally selected rank "
                    "statistics metric not compatible with partition approach "
                    "to unordered predictors.");
            const double p_value = best_statistic_by_value(
                split_key, node_key, data, sample_keys,
                best_statistic, best_split_key, best_value
            );
            if (p_value >= 0) {
                p_values.emplace_back(p_value); keys.emplace_back(split_key);
            }
        }
        if (p_values.empty()) break;
        p_values = adjust_pvalues(p_values);
        best_decrease = std::max(
            best_decrease,
            -p_values[std::find(keys.cbegin(), keys.cend(), best_split_key) -
                          keys.cbegin()]
        );
    } break;
    default: { throw std::invalid_argument("Invalid split metric."); }
    } /* switch (split_rule) */

  /* Clean up aggregate data */
    finalise_node_aggregates();

    if (best_decrease < min_metric_decrease) return false;

    split_keys[node_key] = best_split_key;
    split_values[node_key] = best_value;

    return true;

}


template <typename ImplT>
void Tree<ImplT>::best_decrease_by_value_extratrees(
    const size_t split_key, const size_t node_key,
    const std::shared_ptr<const Data> data, const key_vector & sample_keys,
    double & best_decrease, size_t & best_split_key, double & best_value
) {

    ImplT & tree_impl = *static_cast<ImplT *>(this);

    const size_t n_sample_node = get_n_sample_node(node_key);
    dbl_vector candidate_values;

  /* Get the candidate values using extra-trees algorithm - i.e. randomly
   * sample between minimum and maximum value */
    double min = 0, max = 0;
    data->get_minmax_values(min, max, sample_keys, split_key,
                            start_pos[node_key], end_pos[node_key]);
    if (min == max) return;

  /* FIXME: should check n_random_split > 0? */
    candidate_values.reserve(n_random_split);
    std::uniform_real_distribution<double> U_rng(min, max);
    for (size_t j = 0; j != n_random_split; ++j)
        candidate_values.emplace_back(U_rng(gen));
    std::sort(candidate_values.begin(), candidate_values.end());
  /* Final element is never considered for splitting in worker. */
    candidate_values.emplace_back(INFINITY);

    const size_t n_candidate_value = candidate_values.size();
    if (n_candidate_value < 2) return;

  /* When a new best decrease is found, use the exact corresponding candidate
   * value as the new best value to split on. */
    auto update_best_value = [&](size_t j){ best_value = candidate_values[j]; };

    prepare_candidate_loop_via_value(split_key, node_key, data, sample_keys,
                                     candidate_values);

    tree_impl.best_decrease_by_real_value(
        split_key, n_sample_node, n_candidate_value,
        best_decrease, best_split_key, update_best_value
    );

    finalise_candidate_loop();

}


template <typename ImplT>
void Tree<ImplT>::best_decrease_by_value_extratrees_unordered(
    const size_t split_key, const size_t node_key,
    const std::shared_ptr<const Data> data,
    const key_vector & sample_keys,
    double & best_decrease, size_t & best_split_key, double & best_value
) {

    ImplT & tree_impl = *static_cast<ImplT *>(this);

    const size_t n_sample_node = get_n_sample_node(node_key);
    const size_t n_candidate_value =
        data->get_n_unique_predictor_value(split_key);

    ull_bitenc is_in_node, is_ex_node;
    for (size_t j = start_pos[node_key]; j != end_pos[node_key]; ++j) {
        const size_t value =
            std::floor(data->get_x(sample_keys[j], split_key)) - 1;
        is_in_node.set(value);
    }
    for (size_t k = 0; k != n_candidate_value; ++k)
        is_ex_node.set(!is_in_node.test(k));

  /* Number of unique partitions is equal to the number of different ways
   * we can select factor levels modulo the equivalence of negation (i.e.
   * swapping selected and not-selected). */
    const size_t n_partition = n_random_split;

  /* The output key is as per Tree::split_node (see comments for unordered); it
   * is drawn randomly from all available partitions that put at least one of
   * the observed levels to the right. */
    auto to_partition_key = [&](size_t j){
        ull_bitenc key = 0;
        { /* don't allow full or empty for splitting on present values */
            const size_t n_partition =
                (2ull << (is_in_node.count() - 1ull)) - 2ull;
            std::uniform_int_distribution<size_t> U_rng(1, n_partition);

            const ull_bitenc drawn_in_partition = U_rng(gen);
            size_t key_j = 0;
          /* For each in-node value, get the bit offset in the partition key,
           * and if the value was selected - set the bit */
            for (size_t k = 0; k != is_in_node.count(); ++k) {
                while (!is_in_node.test(key_j)) ++key_j;
                if (drawn_in_partition.test(k)) key.set(key_j);
            }
        }
        { /* allow full or empty for splitting on non-present values */
            const size_t n_partition =
                (2ull << (is_ex_node.count() - 1ull)) - 1ull;
            std::uniform_int_distribution<size_t> U_rng(0, n_partition);

            const ull_bitenc drawn_ex_partition = U_rng(gen);
            size_t key_j = 0;
          /* For each ex-node value, get the bit offset in the partition key,
           * and if the value was selected - set the bit */
            for (size_t k = 0; k != is_ex_node.count(); ++k) {
                while (!is_ex_node.test(key_j)) ++key_j;
                if (drawn_ex_partition.test(k)) key.set(key_j);
            }
        }
        return key;
    };

    tree_impl.best_decrease_by_partition(
        split_key, node_key, data, sample_keys,
        n_sample_node, n_partition, to_partition_key,
        best_decrease, best_split_key, best_value
    );

}


template <typename ImplT>
void Tree<ImplT>::best_decrease_by_value_smallq(
    const size_t split_key, const size_t node_key,
    const std::shared_ptr<const Data> data, const key_vector & sample_keys,
    double & best_decrease, size_t & best_split_key, double & best_value
) {

    ImplT & tree_impl = *static_cast<ImplT *>(this);

    const size_t n_sample_node = get_n_sample_node(node_key);
  /* All values of predictors that are present in this node are candidates for
   * splitting. */
    dbl_vector candidate_values;
    data->get_all_values(candidate_values, sample_keys, split_key,
                         start_pos[node_key], end_pos[node_key]);

    const size_t n_candidate_value = candidate_values.size();
  /* Break if pure or empty node. */
    if (n_candidate_value < 2) return;

  /* When a new best decrease is found, use the mid-point of the candidate
   * values to split on - except if the candidate values are too far apart and
   * the numerical difference between the right-hand value and the mid-point is
   * zero (then we'll use the left-hand value). */
    auto update_best_value = [&](size_t j){
        const double x0 = candidate_values[j];
        const double x1 = candidate_values[j + 1];
        best_value = (x0 + x1) / 2;
        if (best_value == x1) best_value = x0;
    };

    prepare_candidate_loop_via_value(split_key, node_key, data, sample_keys,
                                     candidate_values);

    tree_impl.best_decrease_by_real_value(
        split_key, n_sample_node, n_candidate_value,
        best_decrease, best_split_key, update_best_value
    );

    finalise_candidate_loop();

}


template <typename ImplT>
void Tree<ImplT>::best_decrease_by_value_largeq(
    const size_t split_key, const size_t node_key,
    const std::shared_ptr<const Data> data,
    const key_vector & sample_keys,
    double & best_decrease, size_t & best_split_key, double & best_value
) {

    ImplT & tree_impl = *static_cast<ImplT *>(this);

    prepare_candidate_loop_via_index(split_key, node_key, data, sample_keys);

  /* Break if pure or empty node. */
    size_t test_n_candidate = 0;
    for (const size_t & n : node_n_by_candidate) {
        test_n_candidate += n > 0;
        if (test_n_candidate == 2) break;
    }
    if (test_n_candidate != 2) return;

    const size_t n_sample_node = get_n_sample_node(node_key);
    const size_t n_candidate_value =
        data->get_n_unique_predictor_value(split_key);

  /* Use mid-point. See comments in best_decrease_by_value_smallq */
    auto update_best_value = [&](size_t j0){
        size_t j1 = j0 + 1;
      /* find the next candidate that is present in the node */
        while (j1 != n_candidate_value && node_n_by_candidate[j1] == 0) ++j1;
      /* should check that j1 != n_candidate_value ... but this should never
       * happen. */
        const double x0 = data->get_unique_predictor_value(split_key, j0);
        const double x1 = data->get_unique_predictor_value(split_key, j1);
        best_value = (x0 + x1) / 2;
        if (best_value == x1) best_value = x0;
    };

    tree_impl.best_decrease_by_real_value(
        split_key, n_sample_node, n_candidate_value,
        best_decrease, best_split_key, update_best_value
    );

    finalise_candidate_loop();

}


template <typename ImplT>
void Tree<ImplT>::best_decrease_by_value_unordered(
    const size_t split_key, const size_t node_key,
    const std::shared_ptr<const Data> data,
    const key_vector & sample_keys,
    double & best_decrease, size_t & best_split_key, double & best_value
) {

    ImplT & tree_impl = *static_cast<ImplT *>(this);

    const size_t n_sample_node = get_n_sample_node(node_key);

  /* Get the values of the predictor observed in this node. */
    dbl_vector candidate_values;
    data->get_all_values(
        candidate_values,
        sample_keys, split_key, start_pos[node_key], end_pos[node_key]
    );

    const size_t n_candidate_value = candidate_values.size();
  /* Break if pure or empty node. */
    if (n_candidate_value < 2) return;

    if (n_candidate_value >= std::numeric_limits<size_t>::digits)
        throw std::domain_error("Too many factor levels to enumerate all "
            "partitions.");

  /* Number of unique partitions is equal to the number of different ways
   * we can select factor levels modulo the equivalence of negation (i.e.
   * swapping selected and not-selected). */
    const size_t n_partition = 1ull << (n_candidate_value - 1);

  /* The input j is an enumeration of all partitions. Each bit in j corresponds
   * to an _offset_ into the candidate values. If the bit is 'on' then the
   * partition puts the corresponding candidate value to the right. The output
   * key is as per Tree::split_node (see comments for unordered). */
    auto to_partition_key = [&](size_t j){
        ull_bitenc key;
        for (size_t k = 0; k != n_candidate_value; ++k) {
            if (j & (1ull << k)) key.set(std::floor(candidate_values[k]) - 1);
        }
        return key;
    };

    tree_impl.best_decrease_by_partition(
        split_key, node_key, data, sample_keys,
        n_sample_node, n_partition, to_partition_key,
        best_decrease, best_split_key, best_value
    );

}


template <typename ImplT>
double Tree<ImplT>::best_statistic_by_value(
    const size_t split_key, const size_t node_key,
    const std::shared_ptr<const Data> data,
    const key_vector & sample_keys,
    double & best_statistic, size_t & best_split_key, double & best_value
) {

    ImplT & tree_impl = *static_cast<ImplT *>(this);

    const size_t n_sample_node = get_n_sample_node(node_key);

  /* All values of predictors that are present in this node are candidates for
   * splitting. */
    dbl_vector candidate_values;
    data->get_all_values(candidate_values, sample_keys, split_key,
                         start_pos[node_key], end_pos[node_key]);

    const size_t n_candidate_value = candidate_values.size();
  /* Break if pure or empty node. */
    if (n_candidate_value < 2) return -INFINITY;

    prepare_candidate_loop_via_value(split_key, node_key, data, sample_keys,
                                     candidate_values);

    /* best statistic and (split) value for this candidate */
    double this_statistic = -INFINITY, this_value = -INFINITY,
           this_p_value   = -INFINITY;
    auto update_this_value = [&](size_t j){
        const double x0 = candidate_values[j];
        const double x1 = candidate_values[j + 1];
        this_value = (x0 + x1) / 2;
        if (this_value == x1) this_value = x0;
    };
    tree_impl.best_statistic_by_real_value(
        n_sample_node, n_candidate_value, this_statistic, update_this_value,
        this_p_value
    );

    if (this_statistic > best_statistic) {
        best_statistic = this_statistic;
        best_value = this_value;
        best_split_key = split_key;
    }

    finalise_candidate_loop();

    return this_p_value;

}


template <typename ImplT>
void Tree<ImplT>::finalise_candidate_loop() {

    if (save_memory) {
      /* NOTE: release of memory may be implementation dependent */
        node_n_by_candidate.clear();
        node_n_by_candidate.shrink_to_fit();
    }

}



} /* namespace literanger */


#endif /* LITERANGER_TREE_DEFN_H */

