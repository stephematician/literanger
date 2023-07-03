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
#ifndef LITERANGER_TREE_CLASSIFICATION_DECL_H
#define LITERANGER_TREE_CLASSIFICATION_DECL_H

/* base class declaration */
#include "Tree.decl.h"

/* standard library headers */
#include <cstddef>
#include <memory>
#include <unordered_map>

/* general literanger headers */
#include "enum_types.h"
#include "globals.h"
/* required literanger class declarations */
#include "Data.decl.h"


namespace literanger {

struct TreeClassification : Tree<TreeClassification> {

    friend class Tree<TreeClassification>;

    public:

        /** Construct a classification tree.
         * @param[in] response_weights An vector of weights for each value of
         * the response (each class).
         * @param[in] parameters Parameters that describe the sampling,
         * drawing, and splitting for trees in a random forest.
         * @param[in] save_memory Indicator whether to aggressively release
         * memory and omit building an index (which takes up memory but speeds
         * up training). */
        TreeClassification(const dbl_vector_ptr response_weights,
                           const TreeParameters & parameters,
                           const bool save_memory);

        const std::unordered_map<size_t,key_vector> & get_leaf_keys() const;

        /** Predict response for a leaf (terminal) node.
         *
         * Uses the response that were in-bag during training of the tree to
         * predict or draw a response.
         *
         * @param[in] node_key Index (key) of node.
         * @param[out] result The predicted or otherwise-drawn value.
         * @tparam prediction_type The enumerated type of prediction to perform
         * e.g. bagged, imputation.
         * @tparam result_type The type of data to return; usually a single
         * value e.g. double.
         */
        template <PredictionType prediction_type, typename result_type,
                  enable_if_bagged<prediction_type> = nullptr>
        void predict_from_inbag(const size_t node_key,
                                result_type & result);

        template <PredictionType prediction_type, typename result_type,
                  enable_if_inbag<prediction_type> = nullptr>
        void predict_from_inbag(const size_t node_key,
                                result_type & result);

        template <PredictionType prediction_type, typename result_type,
                  enable_if_nodes<prediction_type> = nullptr>
        void predict_from_inbag(const size_t node_key,
                                result_type & result);

    protected:

        /** A container of the weight for each response value. */
        const dbl_vector_ptr response_weights;

        /** The number of unique response values. */
        const size_t n_response_value = response_weights->size();

        /** The count of each response value for the node being evaluated. */
        dbl_vector node_n_by_response = dbl_vector(n_response_value, 0);

        /** Count of the number of observations for each combination of
         * candidate split value and response. */
        mutable count_vector node_n_by_candidate_and_response;

        /** A map from node keys for each leaf node to the observation keys
         *  that was in-bag for the node during growth (training). */
        std::unordered_map<size_t,key_vector> leaf_keys;

        /** A map from node keys for each leaf node to the most frequent value
         *  that was in-bag for the node during growth (training). */
        mutable std::unordered_map<size_t,double> leaf_most_frequent;


    private:

        /** Implements the response-wise bootstrap sampling for a classification
         * tree.
         * @param[in] data Data to grow (or train) tree with. Contains
         * observations of predictors and the response, the former has
         * predictors across columns and observations by row, and the latter is
         * usually a column vector (or matrix).
         * @param[out] sample_keys A container of randomly drawn keys (i.e.
         * row-offsets) of observations.
         * @param[out] inbag_counts A container of counts of the number of
         * timees each observation appears in-bag.
         */
        void resample_response_wise_impl(const std::shared_ptr<const Data> data,
                                         key_vector & sample_keys,
                                         count_vector & inbag_counts);

        /** Prepare a classification tree for growth by reserving space for
         * terminal nodes.
         *
         * @param[in] data Data to grow (or train) tree with. Contains
         * observations of predictors and the response, the former has
         * predictors across columns and observations by row, and the latter is
         * usually a column vector (or matrix).
         */
        void new_growth(const std::shared_ptr<const Data> data);

        /** @copydoc TreeBase::add_terminal_node() */
        void add_terminal_node(const size_t node_key,
                               const std::shared_ptr<const Data> data,
                               const key_vector & sample_keys);

        /** @copydoc TreeBase::compare_response() */
        bool compare_response(const std::shared_ptr<const Data> data,
                              const size_t lhs_key, const size_t rhs_key) const;

        /** @copydoc Tree::new_node_aggregates() */
        void new_node_aggregates(
            const size_t node_key,
            const std::shared_ptr<const Data> data,
            const key_vector & sample_keys
        );

        /** @copydoc Tree::finalise_node_aggregates() */
        void finalise_node_aggregates();

        /** @copydoc Tree::prepare_candidate_loop_via_value() */
        void prepare_candidate_loop_via_value(
            const size_t split_key, const size_t node_key,
            const std::shared_ptr<const Data> data,
            const key_vector & sample_keys, const dbl_vector & candidate_values
        );

        /** @copydoc Tree::prepare_candidate_loop_via_index() */
        void prepare_candidate_loop_via_index(
            const size_t split_key, const size_t node_key,
            const std::shared_ptr<const Data> data,
            const key_vector & sample_keys
        );

        /** @copydoc Tree::finalise_loop_invariants() */
        void finalise_candidate_loop();

        /** Search the real-valued split candidates for the best decrease in
         * impurity and update the current best key, value, and decrease.
         *
         * @param[in] split_key Identifies which predictor to evaluate.
         * @param[in] n_sample_node The number of observations in the node.
         * @param[in] n_candidate_value The number of candidate values for
         * splitting.
         * @param[in,out] best_decrease The best decrease in node impurity
         * achieved by splitting.
         * @param[in,out] best_split_key The predictor which gave the best
         * decrease in node impurity.
         * @param[in] update_best_value A function that updates the best value
         * given an index into the candidate value vector.
         */
        template <typename UpdateT>
        void best_decrease_by_real_value(
            const size_t split_key,
            const size_t n_sample_node, const size_t n_candidate_value,
            double & best_decrease, size_t & best_split_key,
            UpdateT update_best_value
        ) const;

        /** Search the partition candidates for the best decrease in impurity
         * and update the current best key, value, and decrease.
         *
         * @param[in] split_key Identifies which predictor to evaluate.
         * @param[in] node_key Identifies the node to evaluate.
         * @param[in] data Data to train forest with. Contains observations of
         * predictors and the response, the former has predictors across
         * columns and observations by row, and the latter is usually a column
         * vector (or matrix).
         * @param[in] sample_keys The partially-sorted keys in the sample for
         * this tree.
         * @param[in] n_sample_node The number of observations in the node.
         * @param[in] n_partition The total number of partitions.
         * @param[in] to_partition_key A function that converts an integer index
         * to a partition bit-mask.
         * @param[in,out] best_decrease The best decrease in node impurity
         * achieved by splitting.
         * @param[in,out] best_split_key The predictor which gave the best
         * decrease in node impurity.
         * @param[in,out] best_value The value to split by that achieved the
         * best decrease in node impurity.
         */
        template <typename CallableT>
        void best_decrease_by_partition(
            const size_t split_key, const size_t node_key,
            const std::shared_ptr<const Data> data,
            const key_vector & sample_keys,
            const size_t n_sample_node,
            const size_t n_partition, CallableT to_partition_key,
            double & best_decrease, size_t & best_split_key, double & best_value
        ) const;

        template <typename UpdateT>
        void best_statistic_by_real_value(
            const size_t n_sample_node, const size_t n_candidate_value,
            double & this_decrease, UpdateT update_this_value,
            double & this_p_value
        );

        /** Evaluates the decrease in node impurity given the counts to the left
         * of the split.
         * @param[in] n_by_response_lhs The count by each response to the left
         * of the split.
         * @param[in] n_lhs The number of observations in the node to the left
         * of the split.
         * @param[in] n_rhs The number of observations in the node to the right
         * of the split.
         */
        double evaluate_decrease(
            const count_vector & n_by_response_lhs,
            const size_t n_lhs, const size_t n_rhs
        ) const;


};


} /* namespace literanger */


#endif /* LITERANGER_TREE_CLASSIFICATION_DECL_H */

