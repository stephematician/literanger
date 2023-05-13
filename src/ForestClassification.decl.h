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
#ifndef LITERANGER_FOREST_CLASSIFICATION_DECL_H
#define LITERANGER_FOREST_CLASSIFICATION_DECL_H

/* base class declaration */
#include "Forest.decl.h"

#include <cstddef>
#include <memory>
#include <vector>

/* general literanger headers */
#include "enum_types.h"
#include "globals.h"
/* required literanger class declarations */
#include "Data.decl.h"


namespace literanger {

/** A classification forest. */
struct ForestClassification : public Forest<ForestClassification> {

    friend class Forest<ForestClassification>;

    public:

        /** Construct a classification forest.
         * @param[in] response_values Unique values of the response in the order
         * that they appear in the data.
         * @param[in] response_weights An optional vector of weights for each
         * value of response in @p response_values.
         * @param[in] tree_parameters ... TODO:
         * @param[in] save_memory Indicator whether to aggressively release
         * memory and omit building an index (which takes up memory but speeds
         * up training). */
        ForestClassification(const dbl_vector_ptr response_values,
                             const dbl_vector_ptr response_weights,
                             const std::vector<TreeParameters> tree_parameters,
                             const bool save_memory);


    protected:

        /** Prepare a workspace for growth phase of the forest.
         *
         * This method prepares the data for growth, chiefly by calling methods
         * of the @p data object to construct helper containers that are
         * managed by @p data.
         *
         * The observed (training) values of the response are mapped to their
         * offset in the ForestClassification::response_values vector and
         * stored in a container managed by @p data.
         *
         * If the sampling (for each tree's in-bag cases) is to be performed
         * response-wise, then a container is constructed which stores the
         * sample-keys by their response, managed by @p data.
         *
         * If ForestClassification::save_memory is false then, similarly, the
         * values of the predictors are mapped to an index which is managed by
         * the @p data object.
         *
         * @param[in,out] data Data to train forest with, see literanger::Data
         * class for further details about format. Helper (mutable) containers
         * that are managed by this object are constructed, see details
         * above. */
        void new_growth(const std::shared_ptr<const Data> data);

        /** Finalise the workspace for growth of the forest.
         *
         * The helper containers that were used during growth are finalised via
         * calling methods of the @p data object.
         *
         * For details of the containers, see ForestClassification::new_growth.
         *
         * @param[in,out] data Data to train forest with, see literanger::Data
         * class for further details about format. Helper (mutable) containers
         * that are managed by this object are finalised. */
        void finalise_growth(const std::shared_ptr<const Data> data);

        /** Plant and grow (train) a single classification tree in the forest.
         * @param[in] data Data to train forest with, see literanger::Data class
         * for further details about format.
         * TODO: parameters
         * @param[in] draw_split_weights Weighting to apply to each predictor
         * when drawing the candidates for splitting nodes. */
        void plant_tree(const std::shared_ptr<const Data> data,
                        const TreeParameters & parameters);

        /** Prepare workspace for out-of-bag error estimation.
         *
         * Allocates space for the predictions of each tree for their out-of-bag
         * samples.
         *
         * @param[in] data Data to train forest with, see literanger::Data class
         * for further details about format.
         * @param[in] n_thread Number of workers (threads) to use. */
        void new_oob_error(const std::shared_ptr<const Data> data,
                           const size_t n_thread);

        /** Finalises out-of-bag error estimation.
         *
         * Calculates the out-of-bag error by bagging the predicted responses
         * from the trees that contributed to each case, then calculates the
         * overall misclassification rate.
         *
         * @param[in] data Data to train forest with, see literanger::Data class
         * for further details about format.
         * @returns The overall misclassification rate in out-of-bag samples. */
        double finalise_oob_error(const std::shared_ptr<const Data> data);

        /** Calculate the out-of-bag predictions for one tree in the forest.
         * @param[in] tree_key The index of the tree to elicit predictions from.
         * @param[in] data Data to train forest with, see literanger::Data class
         * for further details about format.
         * @param[in] oob_keys The sample keys for the out-of-bag samples
         * (cases) of the tree. */
        void oob_one_tree(const size_t tree_key,
                          const std::shared_ptr<const Data> data,
                          const key_vector & oob_keys);

        /** Prepare workspace for prediction using a trained forest.
         *
         * Constructs containers for the predictions as they are generated by
         * each tree.
         *
         * This method is templated by prediction type. For bagged prediction,
         * a container with the final aggregate predictions is constructed
         * (with zero entries). For predictions used in imputation, a container
         * which assigns samples to trees is generated.
         *
         * @param[in] data Data to train forest with, see literanger::Data class
         * for further details about format.
         * @param[in] n_thread Number of workers (threads) to use.
         * @tparam prediction_type The enumerated type of predictions to
         * calculate. */
        template <PredictionType prediction_type>
        void new_predictions(const std::shared_ptr<const Data> data,
                             const size_t n_thread);

        /** Finalise bagged predictions of the forest.
         * @param[out] result The bagged predictions for each case.
         * @tparam prediction_type The enumerated type of predictions to
         * calculate.
         * @tparam result_type The type for the returned data.
         * @tparam enable_if_bagged<prediction_type> Substitution success for
         * PredictionType::BAGGED - enables partial specialisation. */
        template <PredictionType prediction_type, typename result_type,
                  enable_if_bagged<prediction_type> = nullptr>
        void finalise_predictions(result_type & result);

        /** Finalise imputation predictions of the forest.
         * @param[out] result The drawn predictions for each case.
         * @tparam prediction_type The enumerated type of predictions to
         * calculate.
         * @tparam result_type The type for the returned data.
         * @tparam enable_if_inbag<prediction_type> Substitution success for
         * PredictionType::INBAG - enables partial specialisation. */
        template <PredictionType prediction_type, typename result_type,
                  enable_if_inbag<prediction_type> = nullptr>
        void finalise_predictions(result_type & result);

        /** Finalise imputation predictions of the forest.
         * @param[out] result The terminal nodes of each tree for each case.
         * @tparam prediction_type The enumerated type of predictions to
         * calculate.
         * @tparam result_type The type for the returned data.
         * @tparam enable_if_nodes<prediction_type> Substitution success for
         * PredictionType::NODES - enables partial specialisation. */
        template <PredictionType prediction_type, typename result_type,
                  enable_if_nodes<prediction_type> = nullptr>
        void finalise_predictions(result_type & result);

        /** Calculate the predictions from one tree in the forest.
         * @param[in] tree_key The index of the tree to elicit predictions from.
         * @param[in] data Data to train forest with, see literanger::Data class
         * for further details about format.
         * @param[in] sample_keys The keys for the cases to predict.
         * @tparam prediction_type The enumerated type of predictions to
         * calculate. */
        template <PredictionType prediction_type>
        void predict_one_tree(const size_t tree_key,
                              const std::shared_ptr<const Data> data,
                              const key_vector & sample_keys);

        /** Aggregate the predictions of one sample in the data set.
         * @param[in] item_key The key for which sample to aggregate. */
        template <PredictionType prediction_type>
        void aggregate_one_item(const size_t item_key);

        /** Values of the response in the order that they appear in the data. */
        const dbl_vector_ptr response_values;

        /** Weights of response values (classes) when computing node impurity.
         *
         * Same ordering as in ForestClassification::response_values. */
        const dbl_vector_ptr response_weights;

        /** The number of unique response values. */
        const size_t n_response_value = response_values->size();

        /** A (workspace) container of the predicted responses for each case
         * whenever that case was out-of-bag during training. */
        std::vector<key_vector> oob_predictions;

        /** A (workspace) container of the predicted responses by trees for
         * each case when prediction type is PredictionType::BAGGED. */
        std::vector<key_vector> predictions_to_bag;

        /** A (workspace) container of indices of cases that will be predicted
         * by each tree when prediction type is PredictionType::INBAG. */
        std::vector<key_vector> prediction_keys_by_tree;

        /** A (workspace) container of the predicted terminal nodes for each
         * case prediction type is PredictionType::NODES. */
        std::vector<key_vector> prediction_nodes;

        /** Container for the final bagged (or otherwise) predictions. */
        dbl_vector aggregate_predictions;


};


} /* namespace literanger */


#endif /* LITERANGER_FOREST_CLASSIFICATION_DECL_H */

