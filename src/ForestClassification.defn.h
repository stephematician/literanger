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
#ifndef LITERANGER_FOREST_CLASSIFICATION_DEFN_H
#define LITERANGER_FOREST_CLASSIFICATION_DEFN_H

/* class declaration */
#include "ForestClassification.decl.h"

/* standard library headers */
#include <iterator>
#include <mutex>
#include <random>
#include <stdexcept>
#include <unordered_map>

/* general literanger headers */
#include "utility.h" // most_frequent_value
/* requred literanger class definitions */
#include "Data.defn.h"
#include "Forest.defn.h"
#include "TreeClassification.defn.h"


namespace literanger {

inline ForestClassification::ForestClassification(
    const dbl_vector_ptr response_values, const dbl_vector_ptr response_weights,
    const std::vector<TreeParameters> tree_parameters, const bool save_memory
) :
    Forest(TREE_CLASSIFICATION, tree_parameters, save_memory),
    response_values(response_values),
    response_weights(
        response_weights->size() != 0 ?
            response_weights :
            dbl_vector_ptr(new dbl_vector(response_values->size(), 1))
    )
{
    if (this->response_weights->size() != n_response_value)
        throw std::invalid_argument("Number of response weights does not match "
            "number of observed response values");

    bool any_hellinger = false;
    for (const auto & parameters : this->tree_parameters) {
        any_hellinger |= parameters.split_rule == HELLINGER;
    }
    if (any_hellinger && n_response_value != 2)
        throw std::invalid_argument("Hellinger metric only implemented for "
            "binary classification.");

}


inline void ForestClassification::plant_tree(
    const std::shared_ptr<const Data> data,
    const TreeParameters & parameters
) {
    trees.emplace_back(new TreeClassification(
      /* classification-specific arguments */
        response_weights,
      /* Forwarded arguments */
        parameters, save_memory)
    );
}


inline void ForestClassification::new_growth(
    const std::shared_ptr<const Data> data
) {
    bool any_by_response = false;
    for (const auto & parameters : tree_parameters) {
        any_by_response |= parameters.sample_fraction->size() > 1;
    }

    data->new_response_index(*response_values);
    if (any_by_response) data->new_sample_keys_by_response();
    if (!save_memory) data->new_predictor_index();
  /* TODO: could add in importance mode shuffling here? */
}


inline void ForestClassification::finalise_growth(
    const std::shared_ptr<const Data> data
) {
    data->finalise_sample_keys_by_response();
    data->finalise_response_index();
    // if (!save_memory) data->finalise_predictor_index();
}


inline void ForestClassification::new_oob_error(
    const std::shared_ptr<const Data> data, const size_t n_thread
) {
    oob_predictions.assign(data->get_n_row(), key_vector());
}


inline double ForestClassification::finalise_oob_error(
    const std::shared_ptr<const Data> data
) {
    // TODO: use response keys instead of values!
  /* for each observation; count the oob predictions by response */
    const size_t n_sample = data->get_n_row();
    std::vector<std::unordered_map<size_t,size_t>> counts { n_sample };

  /* loop over each sample key and record each prediction (by value) */
    for (size_t sample_key = 0; sample_key != n_sample; ++sample_key) {
        for (const size_t & response : oob_predictions[sample_key])
            ++counts[sample_key][response];
    }

    size_t n_misclassification = 0, n_prediction = 0;

    for (size_t sample_key = 0; sample_key != n_sample; ++sample_key) {
        const size_t observed = data->get_response_index()[sample_key];
        if (counts[sample_key].empty()) continue;
        const size_t predicted = most_frequent_value(counts[sample_key], gen);
        if (predicted != observed) ++n_misclassification;
        ++n_prediction;
    }

    oob_predictions.clear();
    oob_predictions.shrink_to_fit(); // only if save memory?

    return (double)n_misclassification / (double)n_prediction;

}


inline void ForestClassification::oob_one_tree(
    const size_t tree_key,
    const std::shared_ptr<const Data> data,
    const key_vector & oob_keys
) {

    TreeClassification & tree_impl =
        *static_cast<TreeClassification *>(trees[tree_key].get());

    const size_t n_oob = oob_keys.size();
    dbl_vector oob_values;
    oob_values.reserve(n_oob);

    for (const size_t & item_key : oob_keys) {
        std::back_insert_iterator<dbl_vector> oob_inserter =
            std::back_inserter(oob_values);
        tree_impl.predict<BAGGED>(data, item_key, oob_inserter);
    }

    {
        std::unique_lock<std::mutex> lock(mutex);
        for (size_t j = 0; j != n_oob; ++j)
            oob_predictions[oob_keys[j]].emplace_back(oob_values[j]);
    }

}


template <>
inline void ForestClassification::new_predictions<BAGGED>(
    const std::shared_ptr<const Data> data, const size_t n_thread
) {
    const size_t n_sample = data->get_n_row();
    predictions_to_bag.assign(n_sample, key_vector());
    for (auto & each_sample : predictions_to_bag) each_sample.reserve(n_tree);
    aggregate_predictions.assign(n_sample, 0);
}


template <PredictionType prediction_type, typename result_type,
          enable_if_bagged<prediction_type>>
void ForestClassification::finalise_predictions(
    result_type & result
) {
    result = aggregate_predictions;

    predictions_to_bag.clear();
    aggregate_predictions.clear();
    predictions_to_bag.shrink_to_fit();
    aggregate_predictions.shrink_to_fit();
}


template <>
inline void ForestClassification::predict_one_tree<BAGGED>(
    const size_t tree_key,
    const std::shared_ptr<const Data> data,
    const key_vector & sample_keys
) {

    TreeClassification & tree_impl =
        *static_cast<TreeClassification *>(trees[tree_key].get());

    const size_t n_predict = sample_keys.size();

    key_vector tree_predictions;
    tree_predictions.reserve(n_predict);

  /* Get the predictions for the tree */
    for (size_t key : sample_keys) {
        std::back_insert_iterator<key_vector> prediction_inserter =
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
inline void ForestClassification::aggregate_one_item<BAGGED>(
    const size_t item_key
) {
    std::unordered_map<size_t,size_t> counts;
    counts.reserve(n_response_value);

    for (const auto value : predictions_to_bag[item_key]) ++counts[value];
    aggregate_predictions[item_key] =
        (*response_values)[most_frequent_value(counts, gen)];
}


template <>
inline void ForestClassification::new_predictions<INBAG>(
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

    aggregate_predictions.assign(n_sample, 0);

}


template <PredictionType prediction_type, typename result_type,
          enable_if_inbag<prediction_type>>
void ForestClassification::finalise_predictions(
    result_type & result
) {
    result = aggregate_predictions;

    prediction_keys_by_tree.clear();
    prediction_keys_by_tree.shrink_to_fit();
    aggregate_predictions.clear();
    aggregate_predictions.shrink_to_fit();
}


template <>
inline void ForestClassification::predict_one_tree<INBAG>(
    const size_t tree_key,
    const std::shared_ptr<const Data> data,
    const key_vector & sample_keys
) {
  /* Each work item is the set of predictions from one tree which were
   * randomly assigned in the new_prediction method */
    TreeClassification & tree_impl =
        *static_cast<TreeClassification *>(trees[tree_key].get());

    const size_t n_predict = prediction_keys_by_tree[tree_key].size();

    key_vector tree_predictions;
    tree_predictions.reserve(n_predict);

  /* Get the predictions for the tree */
    for (size_t item_key : prediction_keys_by_tree[tree_key]) {
        std::back_insert_iterator<key_vector> prediction_inserter =
            std::back_inserter(tree_predictions);
        tree_impl.predict<INBAG>(data, item_key, prediction_inserter);
    }

  /* Copy the predictions to aggregate-prediction container directly (the
   * aggregation step will do nothing) */
    {
        std::unique_lock<std::mutex> lock(mutex);
        for (size_t j = 0; j != n_predict; ++j) {
            const size_t sample_key = prediction_keys_by_tree[tree_key][j];
            aggregate_predictions[sample_key] =
                (*response_values)[tree_predictions[j]];
        }
    }

}


template <>
inline void ForestClassification::aggregate_one_item<INBAG>(
    const size_t item_key
) { }


template <>
inline void ForestClassification::new_predictions<NODES>(
    const std::shared_ptr<const Data> data, const size_t n_thread
) {
    const size_t n_sample = data->get_n_row();
    prediction_nodes.assign(n_sample, key_vector());
    for (auto & each_sample : prediction_nodes) each_sample.assign(n_tree, 0);
}


template <PredictionType prediction_type, typename result_type,
          enable_if_nodes<prediction_type>>
void ForestClassification::finalise_predictions(
    result_type & result
) {
    result = prediction_nodes;
    prediction_nodes.clear();
    prediction_nodes.shrink_to_fit();
}


template <>
inline void ForestClassification::predict_one_tree<NODES>(
    const size_t tree_key,
    const std::shared_ptr<const Data> data,
    const key_vector & sample_keys
) {

    TreeClassification & tree_impl =
        *static_cast<TreeClassification *>(trees[tree_key].get());

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
inline void ForestClassification::aggregate_one_item<NODES>(
    const size_t item_key
) { }

} /* namespace literanger */


#endif /* LITERANGER_FOREST_CLASSIFICATION_DEFN_H */

