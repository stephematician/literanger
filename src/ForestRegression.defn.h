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
#ifndef LITERANGER_FOREST_REGRESSION_DEFN_H
#define LITERANGER_FOREST_REGRESSION_DEFN_H

/* class declaration */
#include "ForestRegression.decl.h"

/* standard library headers */
#include <algorithm>
#include <cmath>
#include <iterator>
#include <random>
#include <stdexcept>
#include <mutex>

/* requred literanger class definitions */
#include "Data.defn.h"
#include "Forest.defn.h"
#include "TreeRegression.defn.h"


namespace literanger {

inline ForestRegression::ForestRegression(
     const double min_prop,
    const std::vector<TreeParameters> tree_parameters, const bool save_memory
) :
    Forest(TREE_REGRESSION, tree_parameters, save_memory), min_prop(min_prop)
{ }


inline void ForestRegression::plant_tree(const std::shared_ptr<const Data> data,
                                         const TreeParameters & parameters) {
    trees.emplace_back(new TreeRegression(
      /* Regression-specific arguments */
        min_prop,
      /* Forwarded arguments */
        parameters, save_memory)
    );
}


inline void ForestRegression::new_growth(
    const std::shared_ptr<const Data> data
) {
    bool any_beta = false;
    for (const auto & parameters : tree_parameters) {
        any_beta |= parameters.split_rule == BETA;
    }
    if (any_beta) {
        for (size_t j = 0; j != data->get_n_row(); ++j) {
            if (data->get_y(j, 0) <= 0 || data->get_y(j, 0) >= 1)
                throw std::domain_error("Beta log-likelihood metric requires "
                    "regression data in the interval (0,1).");
        }
    }

    if (!save_memory) data->new_predictor_index();
}


inline void ForestRegression::finalise_growth(
    const std::shared_ptr<const Data> data
) { /* if (!save_memory) data->finalise_predictor_index(); */ }


inline void ForestRegression::new_oob_error(
    const std::shared_ptr<const Data> data, const size_t n_thread
) {
    oob_predictions.assign(data->get_n_row(), dbl_vector());
}


inline double ForestRegression::finalise_oob_error(
    const std::shared_ptr<const Data> data
) {

  /* for each observation; count the oob predictions by response */
    const size_t n_sample = data->get_n_row();
    std::vector<double> sums(n_sample);

    std::transform(
        oob_predictions.cbegin(), oob_predictions.cend(),
        sums.begin(),
        [](const decltype(oob_predictions)::const_iterator::value_type & oob_j){
            return std::accumulate(oob_j.cbegin(), oob_j.cend(), 0.0);
        }
    );

    double ssq = 0;
    size_t n_prediction = 0;

    for (size_t sample_key = 0; sample_key != n_sample; ++sample_key) {
        const size_t n_oob = oob_predictions[sample_key].size();
        if (n_oob > 0) {
            ++n_prediction;
            const double predicted = sums[sample_key] / n_oob;
            const double observed = data->get_y(sample_key, 0);
            ssq += (predicted - observed) * (predicted - observed);
        }
    }

    return ssq / n_prediction;

}


inline void ForestRegression::oob_one_tree(
    const size_t tree_key,
    const std::shared_ptr<const Data> data,
    const key_vector & oob_keys
) {

    TreeRegression & tree_impl =
        *static_cast<TreeRegression *>(trees[tree_key].get());

    const size_t n_oob = oob_keys.size();
    dbl_vector oob_values;
    oob_values.reserve(n_oob);

    for (auto key : oob_keys) {
        std::back_insert_iterator<dbl_vector> oob_inserter =
            std::back_inserter(oob_values);
        tree_impl.template predict<BAGGED>(data, key, oob_inserter);
    }

    {
        std::unique_lock<std::mutex> lock(mutex);
        for (size_t j = 0; j != n_oob; ++j)
            oob_predictions[oob_keys[j]].emplace_back(oob_values[j]);
    }

}


template <>
inline void ForestRegression::new_predictions<BAGGED>(
    const std::shared_ptr<const Data> data, const size_t n_thread
) {
    const size_t n_sample = data->get_n_row();
    predictions_to_bag.assign(n_sample, dbl_vector());
    for (auto & each_sample : predictions_to_bag) each_sample.reserve(n_tree);
    aggregate_predictions.assign(n_sample, 0);
}


template <PredictionType prediction_type, typename result_type,
          enable_if_bagged<prediction_type>>
void ForestRegression::finalise_predictions(
    result_type & result
) {
    result = aggregate_predictions;

    predictions_to_bag.clear();
    aggregate_predictions.clear();
    predictions_to_bag.shrink_to_fit();
    aggregate_predictions.shrink_to_fit();
}


template <>
inline void ForestRegression::predict_one_tree<BAGGED>(
    const size_t tree_key,
    const std::shared_ptr<const Data> data,
    const key_vector & sample_keys
) {

    TreeRegression & tree_impl =
        *static_cast<TreeRegression *>(trees[tree_key].get());

    const size_t n_predict = sample_keys.size();

    dbl_vector tree_predictions;
    tree_predictions.reserve(n_predict);

  /* Get the predictions for the tree */
    for (size_t key : sample_keys) {
        std::back_insert_iterator<dbl_vector> prediction_inserter =
            std::back_inserter(tree_predictions);
        tree_impl.predict<BAGGED>(data, key, prediction_inserter);
    }

  /* Copy the set of predictions for this tree to the container that will be
   * accessed in the aggregation step */
    {
        std::unique_lock<std::mutex> lock(mutex);
        for (size_t key : sample_keys)
            predictions_to_bag[key].emplace_back(tree_predictions[key]);
    }

}


template <>
inline void ForestRegression::aggregate_one_item<BAGGED>(
    const size_t item_key
) {
    const dbl_vector & responses = predictions_to_bag[item_key];
  /* Calculate mean */
    const double sum = std::accumulate(responses.cbegin(), responses.cend(),
                                       0.0);
    aggregate_predictions[item_key] = sum / responses.size();
}


template <>
inline void ForestRegression::new_predictions<INBAG>(
    const std::shared_ptr<const Data> data, const size_t n_thread
) {

    const size_t n_sample = data->get_n_row();
    prediction_keys_by_tree.assign(n_tree, key_vector());

  /* Randomly assign samples to trees */
    std::uniform_int_distribution<size_t> U_rng(0, n_tree  - 1);
    for (size_t sample_key = 0; sample_key != n_sample; ++sample_key) {
        const size_t tree_key = U_rng(gen);
        prediction_keys_by_tree[tree_key].push_back(sample_key);
    }

    aggregate_predictions.assign(n_sample, 0.0);

}


template <PredictionType prediction_type, typename result_type,
          enable_if_inbag<prediction_type>>
void ForestRegression::finalise_predictions(
    result_type & result
) {
    result = aggregate_predictions;

    prediction_keys_by_tree.clear();
    prediction_keys_by_tree.shrink_to_fit();
    aggregate_predictions.clear();
    aggregate_predictions.shrink_to_fit();
}


template <>
inline void ForestRegression::predict_one_tree<INBAG>(
    const size_t tree_key,
    const std::shared_ptr<const Data> data,
    const key_vector & sample_keys
) {
  /* Each work item is the set of predictions from one tree which were
   * randomly assigned in the new_prediction method */
    TreeRegression & tree_impl =
        *static_cast<TreeRegression *>(trees[tree_key].get());

    const size_t n_predict = prediction_keys_by_tree[tree_key].size();

    dbl_vector tree_predictions;
    tree_predictions.reserve(n_predict);

  /* Get the predictions for the tree */
    for (size_t key : prediction_keys_by_tree[tree_key]) {
        std::back_insert_iterator<dbl_vector> prediction_inserter =
            std::back_inserter(tree_predictions);
        tree_impl.predict<INBAG>(data, key, prediction_inserter);
    }

  /* Copy the predictions to aggregate-prediction container directly (the
   * aggregation step will do nothing) */
    {
        std::unique_lock<std::mutex> lock(mutex);
        for (size_t j = 0; j != n_predict; ++j) {
            const size_t sample_key = prediction_keys_by_tree[tree_key][j];
            aggregate_predictions[sample_key] = tree_predictions[j];
        }
    }

}


template <>
inline void ForestRegression::aggregate_one_item<INBAG>(
    const size_t item_key
) { }


template <>
inline void ForestRegression::new_predictions<NODES>(
    const std::shared_ptr<const Data> data, const size_t n_thread
) {
    const size_t n_sample = data->get_n_row();
    prediction_nodes.assign(n_sample, key_vector());
    for (auto & each_sample : prediction_nodes) each_sample.assign(n_tree, 0);
}


template <PredictionType prediction_type, typename result_type,
          enable_if_nodes<prediction_type>>
void ForestRegression::finalise_predictions(
    result_type & result
) {
    result = prediction_nodes;
    prediction_nodes.clear();
    prediction_nodes.shrink_to_fit();
}


template <>
inline void ForestRegression::predict_one_tree<NODES>(
    const size_t tree_key,
    const std::shared_ptr<const Data> data,
    const key_vector & sample_keys
) {

    TreeRegression & tree_impl =
        *static_cast<TreeRegression *>(trees[tree_key].get());

    const size_t n_predict = sample_keys.size();

    key_vector tree_predictions;
    tree_predictions.reserve(n_predict);

  /* Get the predictions for the tree */
    for (size_t key : sample_keys) {
        std::back_insert_iterator<key_vector> prediction_inserter =
            std::back_inserter(tree_predictions);
        tree_impl.predict<NODES>(data, key, prediction_inserter);
    }

  /* Copy the set of predictions for this tree to the container that will be
   * accessed in the aggregation step */
    {
        std::unique_lock<std::mutex> lock(mutex);
        for (size_t key : sample_keys)
            prediction_nodes[key][tree_key] = tree_predictions[key];
    }

}


template <>
inline void ForestRegression::aggregate_one_item<NODES>(
    const size_t item_key
) { }


} /* namespace literanger */


#endif /* LITERANGER_FOREST_REGRESSION_DEFN_H */

