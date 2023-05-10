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
#ifndef LITERANGER_TREE_DECL_H
#define LITERANGER_TREE_DECL_H

/* base class declaration */
#include "TreeBase.decl.h"

/* standard library headers */
#include <cstddef>
#include <memory>

/* general literanger headers */
#include "enum_types.h"
#include "globals.h"
/* required literanger class declarations */
#include "Data.decl.h"


namespace literanger {

template <typename ImplT>
struct Tree : TreeBase {

    /* Allow access to base from implementation */
    friend ImplT;

    public:

        /** Predict or draw a response for a given case in the data.
         *
         * This method recurses through the tree to determine the leaf node that
         * a case belongs to. The method then calls the implementation-specific
         * method to predict or draw, depending on the prediction type, a value
         * given the inbag training data belonging to the leaf.
         *
         * @param[in] data Data to predict repsonse from. Contains the values of
         * predictors (by column) for each observation or new case (rows).
         * @param[in] sample_key The row (or case) to predict response for.
         * @param[out] result The predicted or otherwise-drawn value.
         * @tparam prediction_type The enumerated type of prediction to perform
         * (e.g. bagged).
         * @tparam result_type The type of data to return; usually a single
         * value e.g. double.
        */
        template <PredictionType prediction_type, typename result_type>
        void predict(const std::shared_ptr<const Data> data,
                     const size_t sample_key,
                     result_type & result);


    protected:

        /** Forwarding constructor.
         * @param args Arguments forwarded to TreeBase constructor.
         * @tparam ArgsT The arguments types of a TreeBase constuctor.
         */
        template <typename... ArgsT>
        Tree(ArgsT &&... args);


    private:

        /** @copydoc TreeBase::push_best_split() */
        bool push_best_split(
            const size_t node_key, const std::shared_ptr<const Data> data,
            const key_vector & sample_keys,
            const key_vector & split_candidate_keys
        );

        /**
         * Find the best-performing value to split a node on given an ordered
         * (numeric or factor) predictor and the extra-random trees rule.
         *
         * Searches a random sample (size `n_random_split`) of candidate values
         * taken from the interval spanning the observed values.
         *
         * @param[in] split_key The predictor to evaluate.
         * @param[in] node_key The node to evaluate.
         * @param[in] data Data to grow (or train) tree with. Contains
         * observations of predictors and the response, the former has
         * predictors across columns and observations by row, and the latter is
         * usually a column vector (or matrix).
         * @param[in] sample_keys The partially-sorted keys in the sample for
         * this tree.
         * @param[in,out] best_decrease The best decrease in node impurity
         * achieved by splitting.
         * @param[in,out] best_split_key The predictor which gave the best
         * decrease in node impurity.
         * @param[in,out] best_value The value to split by that achieved the
         * best decrease in node impurity.
         */
        void best_decrease_by_value_extratrees(
            const size_t split_key, const size_t node_key,
            const std::shared_ptr<const Data> data,
            const key_vector & sample_keys,
            double & best_decrease, size_t & best_split_key, double & best_value
        );

        /** Find the best-performing value to split a node on given an unordered
         * (factor) predictor and the extra-random trees rule.
         *
         * Searches a random sample (size `n_random_split`) of _all_ possible
         * partitions where at least one observed value is on the right-hand
         * side.
         *
         * @param[in] split_key The predictor to evaluate.
         * @param[in] node_key The node to evaluate.
         * @param[in] data Data to grow (or train) tree with. Contains
         * observations of predictors and the response, the former has
         * predictors across columns and observations by row, and the latter is
         * usually a column vector (or matrix).
         * @param[in] sample_keys The partially-sorted keys in the sample for
         * this tree.
         * @param[in,out] best_decrease The best decrease in node impurity
         * achieved by splitting.
         * @param[in,out] best_split_key The predictor which gave the best
         * decrease in node impurity.
         * @param[in,out] best_value The value to split by that achieved the
         * best decrease in node impurity.
         */
        void best_decrease_by_value_extratrees_unordered(
            const size_t split_key, const size_t node_key,
            const std::shared_ptr<const Data> data,
            const key_vector & sample_keys,
            double & best_decrease, size_t & best_split_key, double & best_value
        );

        /** Find the best-performing value to split a node on given an ordered
         * (factor or numeric) predictor for a SMALL sample-to-predictor ratio.
         *
         * Searches all mid-points between observed values.
         *
         * @param[in] split_key The predictor to evaluate.
         * @param[in] node_key The node to evaluate.
         * @param[in] data Data to grow (or train) tree with. Contains
         * observations of predictors and the response, the former has
         * predictors across columns and observations by row, and the latter is
         * usually a column vector (or matrix).
         * @param[in] sample_keys The partially-sorted keys in the sample for
         * this tree.
         * @param[in,out] best_decrease The best decrease in node impurity
         * achieved by splitting.
         * @param[in,out] best_split_key The predictor which gave the best
         * decrease in node impurity.
         * @param[in,out] best_value The value to split by that achieved the
         * best decrease in node impurity.
         */
        void best_decrease_by_value_smallq(
            const size_t split_key, const size_t node_key,
            const std::shared_ptr<const Data> data,
            const key_vector & sample_keys,
            double & best_decrease, size_t & best_split_key, double & best_value
        );

        /** Find the best-performing value to split a node on given an ordered
         * (factor or numeric) predictor for a LARGE sample-to-predictor ratio.
         *
         * Searches all mid-points between observed values.
         *
         * @param[in] split_key The predictor to evaluate.
         * @param[in] node_key The node to evaluate.
         * @param[in] data Data to grow (or train) tree with. Contains
         * observations of predictors and the response, the former has
         * predictors across columns and observations by row, and the latter is
         * usually a column vector (or matrix).
         * @param[in] sample_keys The partially-sorted keys in the sample for
         * this tree.
         * @param[in,out] best_decrease The best decrease in node impurity
         * achieved by splitting.
         * @param[in,out] best_split_key The predictor which gave the best
         * decrease in node impurity.
         * @param[in,out] best_value The value to split by that achieved the
         * best decrease in node impurity.
         */
        void best_decrease_by_value_largeq(
            const size_t split_key, const size_t node_key,
            const std::shared_ptr<const Data> data,
            const key_vector & sample_keys,
            double & best_decrease, size_t & best_split_key, double & best_value
        );

        /** Find the best-performing value to split a node on given an unordered
         * (factor) predictor.
         *
         * Searches all partitions of observed values.
         *
         * @param[in] split_key The predictor to evaluate.
         * @param[in] node_key The node to evaluate.
         * @param[in] data Data to grow (or train) tree with. Contains
         * observations of predictors and the response, the former has
         * predictors across columns and observations by row, and the latter is
         * usually a column vector (or matrix).
         * @param[in] sample_keys The partially-sorted keys in the sample for
         * this tree.
         * @param[in,out] best_decrease The best decrease in node impurity
         * achieved by splitting.
         * @param[in,out] best_split_key The predictor which gave the best
         * decrease in node impurity.
         * @param[in,out] best_value The value to split by that achieved the
         * best decrease in node impurity.
         */
        void best_decrease_by_value_unordered(
            const size_t split_key, const size_t node_key,
            const std::shared_ptr<const Data> data,
            const key_vector & sample_keys,
            double & best_decrease, size_t & best_split_key, double & best_value
        );

         /** Find the best-performing value to split a node on given an ordered
         * (factor or numeric) predictor via the max-stat rule.
         *
         * Searches mid-points between observed values.
         *
         * @param[in] split_key The predictor to evaluate.
         * @param[in] node_key The node to evaluate.
         * @param[in] data Data to grow (or train) tree with. Contains
         * observations of predictors and the response, the former has
         * predictors across columns and observations by row, and the latter is
         * usually a column vector (or matrix).
         * @param[in] sample_keys The partially-sorted keys in the sample for
         * this tree.
         * @param[in,out] best_statistic The best statistic for node impurity
         * achieved by splitting.
         * @param[in,out] best_split_key The predictor which gave the best
         * node impurity test statiatic.
         * @param[in,out] best_value The value to split by that achieved the
         * best node impurity test statistic.
         * @returns The p-value of the best statistic for the candidate
         * predictor. */
        double best_statistic_by_value(
            const size_t split_key, const size_t node_key,
            const std::shared_ptr<const Data> data,
            const key_vector & sample_keys,
            double & best_statistic, size_t & best_split_key,
            double & best_value
        );

        /** Prepare aggregate data used to search for best node split.
         * @param[in] node_key The node to evaluate.
         * @param[in] data Data to grow (or train) tree with. Contains
         * observations of predictors and the response, the former has
         * predictors across columns and observations by row, and the latter is
         * usually a column vector (or matrix).
         * @param[in] sample_keys The partially-sorted keys in the sample for
         * this tree.
         */
        virtual void new_node_aggregates(
            const size_t node_key,
            const std::shared_ptr<const Data> data,
            const key_vector & sample_keys
        ) = 0;

        /** Clean up any aggregate data for node */
        virtual void finalise_node_aggregates() = 0;

        /** Prepares the loop-invariants needed to evaluate the decrease for
         * each candidate value.
         *
         * The implementation should set the value of `node_n_by_candidate`.
         * Other invariants may also be evaluated. This is called before
         * best_decrease_by_real_value().
         *
         * @param[in] split_key The predictor to evaluate.
         * @param[in] node_key The node to evaluate.
         * @param[in] data Data to grow (or train) tree with.
         * @param[in] sample_keys The partially-sorted keys in the sample for
         * this tree.
         * @param[in] candidate_values The values that are candidates for
         * splitting.
         */
        virtual void prepare_candidate_loop_via_value(
            const size_t split_key, const size_t node_key,
            const std::shared_ptr<const Data> data,
            const key_vector & sample_keys, const dbl_vector & candidate_values
        ) = 0;

        /** Prepares the loop-invariants needed to evaluate the decrease for
         * each candidate value.
         *
         * The implementation should set the value of `node_n_by_candidate`.
         * Other invariants may also be evaluated. This is called before
         * `best_decrease_by_real_value`. Candidate values are not passed in
         * this interface, instead the candidate values should be acquired from
         * the data (usually via `data->get_index`).
         *
         * @param[in] split_key The predictor to evaluate.
         * @param[in] node_key The node to evaluate.
         * @param[in] data Data to grow (or train) tree with.
         * @param[in] sample_keys The partially-sorted keys in the sample for
         * this tree.
         */
        virtual void prepare_candidate_loop_via_index(
            const size_t split_key, const size_t node_key,
            const std::shared_ptr<const Data> data,
            const key_vector & sample_keys
        ) = 0;

        virtual void finalise_candidate_loop();


};


} /* namespace literanger */


#endif /* LITERANGER_TREE_DECL_H */

