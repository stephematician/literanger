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
#ifndef LITERANGER_TREE_CLASSIFICATION_DEFN_H
#define LITERANGER_TREE_CLASSIFICATION_DEFN_H

/* class declaration */
#include "TreeClassification.decl.h"

/* standard library headers */
#include <algorithm>
#include <cmath>
#include <iterator>
#include <random>
#include <stdexcept>
#include <vector>

/* general literanger headers */
#include "utility.h"
/* required literanger class definitions */
#include "Data.defn.h"
#include "Tree.defn.h"


namespace literanger {

inline TreeClassification::TreeClassification(
    const dbl_vector_ptr response_weights,
    const TreeParameters & parameters,
    const bool save_memory
) :
    Tree(parameters, save_memory), response_weights(response_weights)
{

    switch (split_rule) {
    case HELLINGER: {
      /* Hellinger rule: applied only to binary classification  */
        if (n_response_value != 2)
            throw std::runtime_error("Cannot use Hellinger metric on "
                "non-binary data.");
    } break;
    case LOGRANK: case EXTRATREES: {
    } break;
    case MAXSTAT: case BETA: {
        throw std::invalid_argument("Unsupported split metric for "
            "classification.");
    } break;
    default: {
        throw std::invalid_argument("Invalid split metric.");
    } break; }

}


inline const std::unordered_map<size_t,key_vector> &
TreeClassification::get_leaf_keys() const { return leaf_keys; }


template <PredictionType prediction_type, typename result_type,
          enable_if_bagged<prediction_type>>
void TreeClassification::predict_from_inbag(
    const size_t node_key,
    result_type & result
) {

    using const_iterator = decltype(leaf_most_frequent)::const_iterator;
    const const_iterator most_frequent_it =
        leaf_most_frequent.find(node_key);
    const bool have_prediction = most_frequent_it !=
                                     leaf_most_frequent.cend();

    if (!have_prediction) {

        std::unordered_map<double,double> counts;
        counts.reserve(n_response_value);
        // TODO: check index here
        for (const size_t & response_key : leaf_keys.at(node_key))
            counts[response_key] += (*response_weights)[response_key];
        if (counts.empty()) return;
        leaf_most_frequent[node_key] = most_frequent_value(counts, gen);
        result = leaf_most_frequent[node_key];

    } else result = most_frequent_it->second;

}


template <PredictionType prediction_type, typename result_type,
          enable_if_inbag<prediction_type>>
void TreeClassification::predict_from_inbag(
    const size_t node_key,
    result_type & result
) {
    // TODO: check weighted - currently as per original ranger (ok?)
    std::uniform_int_distribution<> U_rng(0,
                                          leaf_keys.at(node_key).size() - 1);
    const size_t bag_key = U_rng(gen);
    result = leaf_keys.at(node_key)[bag_key];

}


template <PredictionType prediction_type, typename result_type,
          enable_if_nodes<prediction_type>>
void TreeClassification::predict_from_inbag(
    const size_t node_key,
    result_type & result
) {
    result = node_key;
}


inline void TreeClassification::resample_response_wise_impl(
    const std::shared_ptr<const Data> data,
    key_vector & sample_keys,
    count_vector & inbag_counts
) {

    using const_iterator = key_vector::const_iterator;
    const size_t n_sample = data->get_n_row();
    const std::vector<key_vector> & sample_keys_by_response =
        data->get_sample_keys_by_response();

    if (replace) {

        double start = 0;
        for (size_t j = 0; j != sample_fraction->size(); ++j) {

            const double end = start + (*sample_fraction)[j];
            const size_t n_inbag_j = n_sample * (round(end) - round(start));
            const size_t n_sample_j = sample_keys_by_response[j].size();

            std::uniform_int_distribution<size_t> U_rng(0, n_sample_j - 1);

            for (size_t k = 0; k != n_inbag_j; ++k) {
                const size_t draw = sample_keys_by_response[j][U_rng(gen)];
                sample_keys.emplace_back(draw);
                ++inbag_counts[draw];
            }

            start = end;

        }

    } else {

        double start = 0;
        key_vector sample_j;
        for (size_t j = 0; j != sample_fraction->size(); ++j) {

            const double end = start + (*sample_fraction)[j];
            const size_t n_inbag_j = n_sample * (round(end) - round(start));
            const size_t n_sample_j = sample_keys_by_response[j].size();
            sample_j.assign(n_sample_j, 0);

          /* Get shuffled values between 0 and number observed for current
           * response */
            std::iota(sample_j.begin(), sample_j.end(), 0);
            std::shuffle(sample_j.begin(), sample_j.end(), gen);

          /* Convert the draw from [0, 1, ...] to sample keys (row offsets) */
            for (size_t k = 0; k != n_sample_j; ++k)
                sample_j[k] = sample_keys_by_response[j][sample_j[k]];
          /* Copy results for the current response to the output */
            sample_keys.insert(std::end(sample_keys),
                               sample_j.cbegin(),
                               sample_j.cbegin() + n_inbag_j);
            const const_iterator start_it = sample_j.cbegin() + n_inbag_j;
            const const_iterator end_it = sample_j.cend();
            for (auto draw_it = start_it; draw_it != end_it; ++draw_it)
                ++inbag_counts[*draw_it];

            start = end;

        }

    }

}


inline void TreeClassification::new_growth(
    const std::shared_ptr<const Data> data
) {
    const size_t n_sample = data->get_n_row();
    leaf_keys.clear();
    leaf_most_frequent.clear();
  /* Guess for maximum number of leaf nodes */
    leaf_keys.reserve(std::ceil(n_sample / (double)min_split_n_sample));
    leaf_most_frequent.reserve(
        std::ceil(n_sample / (double)min_split_n_sample)
    );
}


inline void TreeClassification::add_terminal_node(
    const size_t node_key,
    const std::shared_ptr<const Data> data,
    const key_vector & sample_keys
) {
    const size_t start = start_pos[node_key];
    const size_t end = end_pos[node_key];
    leaf_keys[node_key].clear();
    leaf_keys[node_key].reserve(end - start);

    for (size_t j = start; j != end; ++j)
        leaf_keys[node_key].emplace_back(
            data->get_response_index()[sample_keys[j]]
        );
}


inline bool TreeClassification::compare_response(
    const std::shared_ptr<const Data> data,
    const size_t lhs_key,
    const size_t rhs_key
) const {
    return data->get_y(lhs_key, 0) == data->get_y(rhs_key, 0);
}


inline void TreeClassification::new_node_aggregates(
    const size_t node_key,
    const std::shared_ptr<const Data> data,
    const key_vector & sample_keys
) {
    const key_vector & response_keys = data->get_response_index();
    std::fill(node_n_by_response.begin(), node_n_by_response.end(), 0);

    for (size_t j = start_pos[node_key]; j != end_pos[node_key]; ++j) {
        const size_t sample_key = sample_keys[j];
        const size_t response_key = response_keys[sample_key];
        ++node_n_by_response[response_key];
    }
}


inline void TreeClassification::finalise_node_aggregates() { }


inline void TreeClassification::prepare_candidate_loop_via_value(
    const size_t split_key, const size_t node_key,
    const std::shared_ptr<const Data> data,
    const key_vector & sample_keys, const dbl_vector & candidate_values
) {

    const key_vector & response_keys = data->get_response_index();
    const size_t n_candidate_value = candidate_values.size();

    node_n_by_candidate_and_response.assign(
        n_candidate_value * n_response_value, 0
    );
    node_n_by_candidate.assign(n_candidate_value, 0);

    for (size_t j = start_pos[node_key]; j != end_pos[node_key]; ++j) {
        const size_t sample_key = sample_keys[j];
        const size_t response_key = response_keys[sample_key];
        const size_t offset = std::distance(
            candidate_values.cbegin(),
            std::lower_bound(candidate_values.cbegin(), candidate_values.cend(),
                             data->get_x(sample_key, split_key))
        );
        ++node_n_by_candidate[offset];
        ++node_n_by_candidate_and_response[offset * n_response_value +
                                               response_key];
    }

}


inline void TreeClassification::prepare_candidate_loop_via_index(
    const size_t split_key, const size_t node_key,
    const std::shared_ptr<const Data> data,
    const key_vector & sample_keys
) {

    const key_vector & response_keys = data->get_response_index();
    const size_t n_candidate_value =
        data->get_n_unique_predictor_value(split_key);

  /* Get counts by candidate (split) value, and by both candidate value and
   * response value. */
    node_n_by_candidate_and_response.assign(
        n_candidate_value * n_response_value, 0
    );
    node_n_by_candidate.assign(n_candidate_value, 0);

    for (size_t j = start_pos[node_key]; j != end_pos[node_key]; ++j) {
        const size_t sample_key = sample_keys[j];
        const size_t offset = data->get_index(sample_key, split_key);
        const size_t response_key = response_keys[sample_key];

        ++node_n_by_candidate[offset];
        ++node_n_by_candidate_and_response[offset * n_response_value +
                                               response_key];
    }

}


inline void TreeClassification::finalise_candidate_loop() {

    Tree::finalise_candidate_loop();

    if (save_memory) {
      /* NOTE: release of memory may be implementation dependent */
        node_n_by_candidate_and_response.clear();
        node_n_by_candidate_and_response.shrink_to_fit();
    }

}


template <typename UpdateT>
void TreeClassification::best_decrease_by_real_value(
    const size_t split_key,
    const size_t n_sample_node, const size_t n_candidate_value,
    double & best_decrease, size_t & best_split_key, UpdateT update_best_value
) const {

  /* NOTE: Pre-condition: n_candidate_value > 1 */
    size_t n_lhs = 0;
    count_vector node_n_by_response_lhs(n_response_value, 0);

    for (size_t j = 0; j != n_candidate_value - 1; ++j) {

        if (node_n_by_candidate[j] == 0) continue;

        n_lhs += node_n_by_candidate[j];
      /* Update the count, by response value, that lie to the left of the
       * current candidate (split) value. */
        for (size_t k = 0; k != n_response_value; ++k)
            node_n_by_response_lhs[k] +=
                node_n_by_candidate_and_response[j * n_response_value + k];

        if (n_lhs < min_leaf_n_sample) continue;

        const size_t n_rhs = n_sample_node - n_lhs;
        if (n_rhs < min_leaf_n_sample) break;

        const double decrease = evaluate_decrease(node_n_by_response_lhs,
                                                  n_lhs, n_rhs);

     /* If the decrease in node impurity has improved - then we update the best
      * split for the node. */
        if (decrease > best_decrease) {
            update_best_value(j);
            best_split_key = split_key;
            best_decrease = decrease;
        }

    }

}


template <typename CallableT>
void TreeClassification::best_decrease_by_partition(
    const size_t split_key, const size_t node_key,
    const std::shared_ptr<const Data> data,
    const key_vector & sample_keys,
    const size_t n_sample_node,
    const size_t n_partition, CallableT to_partition_key,
    double & best_decrease, size_t & best_split_key, double & best_value
) const {

    const key_vector & response_keys = data->get_response_index();

  /* Start from one (cannot have empty lhs) */
    for (size_t j = 1; j != n_partition; ++j) {

      /* Get the bit-encoded partition value */
        ull_bitenc partition_key = to_partition_key(j);

        count_vector node_n_by_response_lhs(n_response_value, 0);
        size_t n_lhs = 0;

        for (size_t k = start_pos[node_key]; k != end_pos[node_key]; ++k) {
            const size_t sample_key = sample_keys[k];
            const size_t response_key = response_keys[sample_key];
            const size_t level_bit = std::floor(
                data->get_x(sample_key, split_key) - 1
            );
            if (!partition_key.test(level_bit)) {
                ++n_lhs;
                ++node_n_by_response_lhs[response_key];
            }
        }
        if (n_lhs < min_leaf_n_sample) continue;

        const size_t n_rhs = n_sample_node - n_lhs;
        if (n_rhs < min_leaf_n_sample) continue;

        const double decrease = evaluate_decrease(node_n_by_response_lhs,
                                                  n_lhs, n_rhs);

        if (decrease > best_decrease) {
            (size_t &)best_value = partition_key.to_ullong();
            best_split_key = split_key;
            best_decrease = decrease;
        }

    }

}


template <typename UpdateT>
void TreeClassification::best_statistic_by_real_value(
    const size_t n_sample_node, const size_t n_candidate_value,
    double & this_decrease, UpdateT update_this_value, double & this_p_value
) { /* NOTE:: Pre-condition - split rule is valid */ }


inline double TreeClassification::evaluate_decrease(
    const count_vector & node_n_by_response_lhs,
    const size_t n_lhs, const size_t n_rhs
) const {

  /* NOTE:: Pre-condition - split rule is valid */

    switch (split_rule) {
    case HELLINGER: {
      /* TPR is the number of 1s on the right divided by the true number
       * of ones; FPR is the number of 0s on the right divided by the
       * true number of zeros */
        const double tpr = (node_n_by_response[1] - node_n_by_response_lhs[1]) /
            (double)node_n_by_response[1];
        const double fpr = (node_n_by_response[0] - node_n_by_response_lhs[0]) /
            (double)node_n_by_response[0];

      /* Decrease of impurity */
        const double a1 = sqrt(tpr) - sqrt(fpr);
        const double a2 = sqrt(1 - tpr) - sqrt(1 - fpr);
        return sqrt(a1 * a1 + a2 * a2);
    } break;
    case LOGRANK: case EXTRATREES: {
      /* Use (weighted) sum of square count to measure node impurity. */
        double sum_lhs_sq = 0;
        double sum_rhs_sq = 0;
        for (size_t k = 0; k != n_response_value; ++k) {
            const double node_n_by_response_rhs_k =
                node_n_by_response[k] - node_n_by_response_lhs[k];
            sum_lhs_sq += (*response_weights)[k] *
                node_n_by_response_lhs[k] * node_n_by_response_lhs[k];
            sum_rhs_sq += (*response_weights)[k] *
                node_n_by_response_rhs_k * node_n_by_response_rhs_k;
        }
        return sum_rhs_sq / n_rhs + sum_lhs_sq / n_lhs;
    } break;
    default: break; }

    return -INFINITY;

}


} /* namespace literanger */


#endif /* LITERANGER_TREE_CLASSIFICATION_DEFN_H */

