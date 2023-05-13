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
#ifndef LITERANGER_FOREST_DECL_H
#define LITERANGER_FOREST_DECL_H

/* base class declaration */
#include "ForestBase.decl.h"

/* standard library headers */
#include <cstddef>
#include <memory>

/* general literanger headers */
#include "enum_types.h"
#include "globals.h"
#include "utility.h" // toggle_print
/* required literanger class declarations */
#include "Data.decl.h"


namespace literanger {

#define BEGIN_WORKER \
    try {
#define END_WORKER                                    \
    } catch (std::exception & e) {                    \
        {                                             \
            std::unique_lock<std::mutex> lock(mutex); \
            interrupted |= true;                      \
            condition_variable.notify_one();          \
        }                                             \
        throw;                                        \
    }

/** Random forest (interface).
 * @tparam ImplT A type that implements Forest interface.  */
template <typename ImplT>
struct Forest : public ForestBase {

    public:

        /** @copydoc ForestBase::plant()
         *
         * This method will construct each tree in the forest in the main thread
         * and then split the forest into intervals of neighbouring trees,
         * growing each tree using one worker (thread) per interval.
         *
         * If @p compute_oob_error is true then the out-of-bag error is
         * estimated by computing a bagged prediction for each observation
         * using the predictions from trees for which that observation was
         * out-of-bag during training.
         *
         * @see Forest::grow_interval */
        void plant(const std::shared_ptr<const Data> data,
                   const dbl_vector_ptr case_weights,
                   const size_t seed, const size_t n_thread,
                   const bool compute_oob_error,
                   const interruptor & user_interrupt,
                   double & oob_error,
                   toggle_print & print_out);

        /** Predict responses using random forest.
         *
         * This method will split the forest into intervals of neighbouring
         * trees and obtain the requested predictions from each tree using one
         * worker (thread) per interval. Similarly, aggregation of the
         * predictions for each interval is performed using one worker per
         * interval. The predictions are then finalised on the main thread.
         *
         * This function is templated by the enumerated prediction type
         * (see literanger::PredictionType). The type is either the traditional
         * bagged prediction or, as needed for imputation algorithms based on
         * random forest, a random draw from the in-bag values of a random
         * (with-replacement) tree in the forest (e.g. Doove et al 2014).
         *
         * Doove, L. L., Van Buuren, S., & Dusseldorp, E. (2014). Recursive
         * partitioning for missing data imputation in the presence of
         * interaction effects. Computational statistics & data analysis, 72,
         * 92-104.
         *
         * @see Forest::predict_interval
         * @see Forest::aggregate_interval
         *
         * @param[in] data The new cases, i.e. observations of predictors, for
         * which predicted responses will be drawn, see literanger::Data class
         * for further details about format.
         * @param[in] seed The seed for the pseudo-random number generator
         * engine.
         * @param[in] n_thread Number of threads to use when growing and
         * predicting.
         * @param[in] user_interrupt An operator that checks for user interrupt.
         * @param[out] result The predictions for each case.
         * @param[out] print_out A toggle-able printer for outputting progress
         * when training or predicting.
         * @tparam prediction_type The enumerated type of prediction to perform
         * (e.g. bagged).
         * @tparam result_type The type of data returned, usually a
         * container. */
        template <PredictionType prediction_type, typename result_type>
        void predict(const std::shared_ptr<const Data> data,
                     const size_t seed, const size_t n_thread,
                     const interruptor & user_interrupt,
                     result_type & result, toggle_print & print_out);


    protected:

        /** Forwarding constructor.
         * @param args Arguments forwarded to ForestBase constructor.
         * @tparam ArgsT The arguments types of a ForestBase constuctor. */
        template <typename... ArgsT>
        Forest(ArgsT &&... args);


    private:

        /** Grows trees in a given interval of the forest.
         *
         * Intervals of neighbouring trees in the forest are grown (trained)
         * by independent workers (threads). This method is run by the worker to
         * iterate over the trees and call the tree's implementation-specific
         * TreeBase::grow() (train) method. If requested (when
         * @p compute_oob_error is true) this also calls the forest
         * type-specific oob_one_tree() for each tree.
         *
         * @see ForestClassification::oob_one_tree
         * @see ForestRegression::oob_one_tree
         *
         * @param[in] work_index The index for this worker into a vector which
         * specifies all the intervals of trees (or work).
         * @param[in] data Data to train forest with, see literanger::Data class
         * for further details about format.
         * @param[in] case_weights The weight for each case (row) in training.
         * @param[in] compute_oob_error Indicator for whether to calculate
         * (otherwise skip calculating) the out-of-bag error. */
        void grow_interval(const size_t work_index,
                           const std::shared_ptr<const Data> data,
                           const dbl_vector_ptr case_weights,
                           const bool compute_oob_error);

        /** Predict responses in a given interval of the forest.
         *
         * The predictions from intervals of neighbouring trees are performed
         * by independent workers (threads). This method is run by the worker to
         * iterate over the trees and call the prediction- and tree-type
         * specific predict_one_tree() method. The predict_one_tree() method
         * will take care of all the predictions required by a tree and store
         * them in a suitable container in the Forest object. The final
         * predictions are the result of a separate aggregation step, see
         * Forest::aggregate_interval.
         *
         * This function is templated by prediction type, much like the method
         * that calls it, see Forest::predict.
         *
         * @param[in] work_index The index for this worker into a vector which
         * specifies all the intervals of trees (or work).
         * @param[in] data The new cases, i.e. observations of predictors, for
         * which predicted responses will be drawn, see literanger::Data class
         * for further details about format.
         * @tparam prediction_type The enumerated type of predictions to
         * calculate. */
        template <PredictionType prediction_type>
        void predict_interval(const size_t work_index,
                              const std::shared_ptr<const Data> data);


};


} /* namespace literanger */


#endif /* LITERANGER_FOREST_DECL_H */

