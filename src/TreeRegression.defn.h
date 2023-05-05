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
#ifndef LITERANGER_TREE_REGRESSION_DEFN_H
#define LITERANGER_TREE_REGRESSION_DEFN_H

/* class declaration */
#include "TreeRegression.decl.h"

/* standard library headers */
#include <algorithm>
#include <cassert>
#include <cmath>
#include <iterator>
#include <limits>
#include <random>
#include <stdexcept>

/* general literanger headers */
#include "utility_math.h"
/* required literanger class definitions */
#include "Data.defn.h"
#include "Tree.defn.h"


namespace literanger {

inline TreeRegression::TreeRegression(const double min_prop,
                                      const TreeParameters & parameters,
                                      const bool save_memory) :
    Tree(parameters, save_memory), min_prop(min_prop)
{ /* should check split rule here? */ }


template <PredictionType prediction_type, typename result_type,
          enable_if_bagged<prediction_type>>
void TreeRegression::predict_from_inbag(
    const size_t node_key,
    result_type & result
) {

    using const_iterator = decltype(leaf_mean)::const_iterator;
    const const_iterator mean_it = leaf_mean.find(node_key);
    const bool have_prediction = mean_it != leaf_mean.cend();

    if (!have_prediction) {

        double leaf_sum = 0;
        for (const double & response : leaf_values.at(node_key))
            leaf_sum += response;
        if (leaf_values.at(node_key).empty()) return;
        leaf_mean[node_key] = leaf_sum / leaf_values.at(node_key).size();
        result = leaf_mean[node_key];

    } else result = mean_it->second;

}


template <PredictionType prediction_type, typename result_type,
          enable_if_doove<prediction_type>>
void TreeRegression::predict_from_inbag(
    const size_t node_key,
    result_type & result
) {

    std::uniform_int_distribution<> U_rng(0,
                                          leaf_values.at(node_key).size() - 1);
    const size_t key = U_rng(gen);
    result = leaf_values.at(node_key)[key];

}


inline void TreeRegression::new_growth(
    const std::shared_ptr<const Data> data
) {
    const size_t n_sample = data->get_n_row();
    leaf_values.clear();
    leaf_mean.clear();
  /* Guess for maximum number of leaf nodes */
    leaf_values.reserve(std::ceil(n_sample / (double)min_split_n_sample));
    leaf_mean.reserve(
        std::ceil(n_sample / (double)min_split_n_sample)
    );
}


inline void TreeRegression::add_terminal_node(
    const size_t node_key,
    const std::shared_ptr<const Data> data,
    const key_vector & sample_keys
) {
    const size_t start = start_pos[node_key];
    const size_t end = end_pos[node_key];
    leaf_values[node_key].clear();
    leaf_values[node_key].reserve(end - start);

    for (size_t j = start; j != end; ++j)
        leaf_values[node_key].emplace_back(data->get_y(sample_keys[j], 0));
}


inline bool TreeRegression::compare_response(
    const std::shared_ptr<const Data> data,
    const size_t lhs_key,
    const size_t rhs_key
) const {
    return data->get_y(lhs_key, 0) == data->get_y(rhs_key, 0);
}


inline void TreeRegression::new_node_aggregates(
    const size_t node_key,
    const std::shared_ptr<const Data> data,
    const key_vector & sample_keys
) {
  /* Compute sum of responses and the sum of the response-squared in node, or in
   * the maximally selected rank statistic case, sum the scores and squared
   * scores. */
    node_sum = node_ssq = 0;
    if (split_rule != MAXSTAT) {
        for (size_t j = start_pos[node_key]; j != end_pos[node_key]; ++j) {
            const size_t sample_key = sample_keys[j];
            const double response = data->get_y(sample_key, 0);
            node_sum += response;
            if (split_rule == BETA) node_ssq += response * response;
        }
    } else {
        for (size_t j = start_pos[node_key]; j != end_pos[node_key]; ++j) {
            const size_t sample_key = sample_keys[j];
            response_scores.emplace_back(data->get_y(sample_key, 0));
        }
        response_scores = rank(response_scores);
        for (const double & score : response_scores) {
            node_sum += score;
            node_ssq += score * score;
        }
    }
}


inline void TreeRegression::finalise_node_aggregates() {
    response_scores.clear();
    if (save_memory) response_scores.shrink_to_fit();
}


inline void TreeRegression::prepare_candidate_loop_via_value(
    const size_t split_key, const size_t node_key,
    const std::shared_ptr<const Data> data,
    const key_vector & sample_keys, const dbl_vector & candidate_values
) {

    const size_t n_candidate_value = candidate_values.size();

    node_n_by_candidate.assign(n_candidate_value, 0);
    node_sum_by_candidate.assign(n_candidate_value, 0);
    if (split_rule == BETA) {
        node_ssq_by_candidate.assign(n_candidate_value, 0);
        response_by_candidate.assign(n_candidate_value, dbl_vector());
    }

    for (size_t j = start_pos[node_key]; j != end_pos[node_key]; ++j) {

        const size_t sample_key = sample_keys[j];
        const double response = split_rule != MAXSTAT ?
            data->get_y(sample_key, 0) : response_scores[j - start_pos[node_key]];
        const size_t offset = std::distance(
            candidate_values.cbegin(),
            std::lower_bound(candidate_values.cbegin(), candidate_values.cend(),
                             data->get_x(sample_key, split_key))
        );

        ++node_n_by_candidate[offset];
        node_sum_by_candidate[offset] += response;
        if (split_rule == BETA) {
            node_ssq_by_candidate[offset] += response * response;
            response_by_candidate[offset].push_back(response);
        }

    }

}


inline void TreeRegression::prepare_candidate_loop_via_index(
    const size_t split_key, const size_t node_key,
    const std::shared_ptr<const Data> data,
    const key_vector & sample_keys
) {

    const size_t n_candidate_value =
        data->get_n_unique_predictor_value(split_key);

    node_n_by_candidate.assign(n_candidate_value, 0);
    node_sum_by_candidate.assign(n_candidate_value, 0);
    if (split_rule == BETA) {
        node_ssq_by_candidate.assign(n_candidate_value, 0);
        response_by_candidate.assign(n_candidate_value, dbl_vector());
    }

    for (size_t j = start_pos[node_key]; j != end_pos[node_key]; ++j) {

        const size_t sample_key = sample_keys[j];
        const double response = data->get_y(sample_key, 0);
        const size_t offset = data->get_index(sample_key, split_key);

        ++node_n_by_candidate[offset];
        node_sum_by_candidate[offset] += response;
        if (split_rule == BETA) {
            node_ssq_by_candidate[offset] += response * response;
            response_by_candidate[offset].push_back(response);
        }

    }

}


inline void TreeRegression::finalise_candidate_loop() {

    Tree::finalise_candidate_loop();

    if (save_memory) {
      /* NOTE: release of memory may be implementation dependent */
        node_sum_by_candidate.clear();
        node_sum_by_candidate.shrink_to_fit();
        node_ssq_by_candidate.clear();
        node_ssq_by_candidate.shrink_to_fit();
        response_by_candidate.clear();
        response_by_candidate.shrink_to_fit();
    }

}


template <typename UpdateT>
void TreeRegression::best_decrease_by_real_value(
    const size_t split_key,
    const size_t n_sample_node, const size_t n_candidate_value,
    double & best_decrease, size_t & best_split_key, UpdateT update_best_value
) const {

    size_t n_lhs = 0;
    double sum_lhs = 0;
    double ssq_lhs = 0;
    assert(n_candidate_value > 1);

    for (size_t j = 0; j != n_candidate_value - 1; ++j) {

        if (node_n_by_candidate[j] == 0) continue;

        n_lhs += node_n_by_candidate[j];
        if (n_lhs < min_leaf_n_sample) continue;

        const size_t n_rhs = n_sample_node - n_lhs;
        if (n_rhs < min_leaf_n_sample) break;

      /* */
        sum_lhs += node_sum_by_candidate[j];
        if (split_rule == BETA) ssq_lhs += node_ssq_by_candidate[j];

        const double sum_rhs = node_sum - sum_lhs;
        const double ssq_rhs = node_ssq - ssq_lhs;
        const double decrease = evaluate_decrease(n_lhs, n_rhs,
                                                  sum_lhs, sum_rhs,
                                                  ssq_lhs, ssq_rhs);

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
void TreeRegression::best_decrease_by_partition(
    const size_t split_key, const size_t node_key,
    const std::shared_ptr<const Data> data,
    const key_vector & sample_keys,
    const size_t n_sample_node,
    const size_t n_partition, CallableT to_partition_key,
    double & best_decrease, size_t & best_split_key, double & best_value
) {

    if (split_rule == BETA) {
        node_n_by_candidate.assign(2, 0);
        response_by_candidate.assign(2, dbl_vector());
    }

  /* Start from one (cannot have empty lhs) */
    for (size_t j = 1; j != n_partition; ++j) {

      /* Get the bit-encoded partition value */
        ull_bitenc partition_key = to_partition_key(j);

        double ssq_lhs = 0;
        double sum_lhs = 0;
        size_t n_lhs = 0;

        for (size_t k = start_pos[node_key]; k != end_pos[node_key]; ++k) {
            const size_t sample_key = sample_keys[k];
            const size_t level_bit = std::floor(
                data->get_x(sample_key, split_key) - 1
            );
            if (!partition_key.test(level_bit)) {
                if (split_rule == BETA)
                    ssq_lhs += data->get_y(sample_key, 0) *
                                   data->get_y(sample_key, 0);
                sum_lhs += data->get_y(sample_key, 0);
                ++n_lhs;
            }
          /* Use the node_n_by_candidate and response_by_candidate containers to
           * store the left and right hand side values */
            if (split_rule == BETA) {
                const size_t j = (size_t)partition_key.test(level_bit);
                ++node_n_by_candidate[j];
                response_by_candidate[j].push_back(data->get_y(sample_key, 0));
            }
        }
        if (n_lhs < min_leaf_n_sample) continue;

        const size_t n_rhs = n_sample_node - n_lhs;
        if (n_rhs < min_leaf_n_sample) continue;

        const double sum_rhs = node_sum - sum_lhs;
        const double ssq_rhs = node_ssq - ssq_lhs;
        const double decrease = evaluate_decrease(n_lhs, n_rhs,
                                                  sum_lhs, sum_rhs,
                                                  ssq_lhs, ssq_rhs);

        if (decrease > best_decrease) {
            (size_t &)best_value = partition_key.to_ullong();
            best_split_key = split_key;
            best_decrease = decrease;
        }

    }

    if (save_memory) {
        node_n_by_candidate.clear();
        node_n_by_candidate.shrink_to_fit();
        response_by_candidate.clear();
        response_by_candidate.shrink_to_fit();
    }

}


template <typename UpdateT>
void TreeRegression::best_statistic_by_real_value(
    const size_t n_sample_node, const size_t n_candidate_value,
    double & this_decrease, UpdateT update_this_value, double & this_p_value
) {

    assert(n_sample_node > 1);

  /* smallest split to consider */
    const size_t min_split = std::max(0.0, n_sample_node * min_prop - 1);

    double sum_lhs = 0;
    size_t n_lhs = 0;
    assert(n_candidate_value > 1);

    for (size_t j = 0; j != n_candidate_value - 1; ++j) {

        if (node_n_by_candidate[j] == 0) continue;

        n_lhs += node_n_by_candidate[j];
        if (n_lhs < std::max(min_leaf_n_sample, min_split)) continue;

        const size_t n_rhs = n_sample_node - n_lhs;
        if (n_rhs < std::max(min_leaf_n_sample, min_split)) break;

        const double sum_rhs = node_sum - sum_lhs, ssq_lhs = 0, ssq_rhs = 0;

        const double decrease = evaluate_decrease(n_lhs, n_rhs,
                                                  sum_lhs, sum_rhs,
                                                  ssq_lhs, ssq_rhs);

        if (decrease > this_decrease) {
            update_this_value(j);
            this_decrease = decrease;
            const double p_value_Lausen92 = maxstat_p_value_Lausen92(
                this_decrease, min_prop
            );
            const double p_value_Lausen94 = maxstat_p_value_Lausen94(
                this_decrease, n_sample_node, node_n_by_candidate, j + 1
            );
            this_p_value = std::min(p_value_Lausen92, p_value_Lausen94);
        }

    }

}


inline double TreeRegression::evaluate_decrease(
    const size_t n_lhs, const size_t n_rhs,
    const double sum_lhs, const double sum_rhs,
    const double ssq_lhs, const double ssq_rhs
) const {

    switch (split_rule) {
    case BETA: {
      /* Need at least two observations per node to estimate parameters for
       * beta distribution. */
        if (n_lhs < 2 || n_rhs < 2) return -INFINITY;

      /* Using E[(x - E[x])^2] = E[X^2] - E[X]^2 and the unbiased estimate of
       * (sample) variance. */
        const double var_lhs = (ssq_lhs - (sum_lhs * sum_lhs) / n_lhs) /
                                   (n_lhs - 1.0);
        const double var_rhs = (ssq_rhs - (sum_rhs * sum_rhs) / n_rhs) /
                                   (n_rhs - 1.0);
        if (var_lhs < std::numeric_limits<double>::epsilon() ||
                var_rhs < std::numeric_limits<double>::epsilon())
            return -INFINITY;

      /* Mean and 'sample size' parameterisation of beta distribution - see
       * wikipedia mu/nu definitions */
        const double mu_lhs = sum_lhs / n_lhs;
        const double nu_lhs = mu_lhs * (1 - mu_lhs) / var_lhs - 1.0;
        const double mu_rhs = sum_rhs / n_rhs;
        const double nu_rhs = mu_rhs * (1 - mu_rhs) / var_rhs - 1.0;

      /* Evaluate the the log-likelihood on left- and right-hand side given the
       * current split candidate */
        double beta_lnL_lhs = 0, beta_lnL_rhs = 0;
        size_t count = 0; /* use this to track left vs right */
        for (size_t j = 0; j != node_n_by_candidate.size(); ++j) {
          /* Select the lhs or rhs parameters */
            const bool is_lhs = count < n_lhs;
            const double mu_j = is_lhs ? mu_lhs : mu_rhs;
            const double nu_j = is_lhs ? nu_lhs : nu_rhs;
            double & result = is_lhs ? beta_lnL_lhs : beta_lnL_rhs;
          /* Record the likelihood for each observation */
            for (const double & response : response_by_candidate[j])
                result += beta_log_likelihood(response, mu_j, nu_j);
            count += node_n_by_candidate[j];
        }

        const double result = beta_lnL_lhs + beta_lnL_rhs;

        return std::isnan(result) ? -INFINITY : result;

    } break;
    case EXTRATREES: case LOGRANK: {
        const double sum_lhs_sq = sum_lhs * sum_lhs;
        const double sum_rhs_sq = sum_rhs * sum_rhs;
        return sum_lhs_sq / n_lhs + sum_rhs_sq / n_rhs;
    } break;
    case MAXSTAT: {
        const double n = n_lhs + n_rhs;
        const double mu = node_sum / n;
        const double var = (node_ssq - node_sum * node_sum / n) / (n - 1);
        const double S = sum_lhs;
        const double E = n_lhs * mu;
        const double V = n_lhs * (double)n_rhs * var / n;
        return std::fabs((S - E) / std::sqrt(V));
    } break;
    case HELLINGER: {
        throw std::invalid_argument("Unsupported split metric for regression.");
    } break;
    default: { throw std::invalid_argument("Invalid split metric."); }
    }

    return -INFINITY;

}


} /* namespace literanger */


#endif /* LITERANGER_TREE_REGRESSION_DEFN_H */

