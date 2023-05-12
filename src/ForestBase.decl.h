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
#ifndef LITERANGER_FOREST_BASE_DECL_H
#define LITERANGER_FOREST_BASE_DECL_H

/* standard library headers */
#include <condition_variable>
#include <cstddef>
#include <memory>
#include <mutex>
#include <random>
#include <string>
#include <vector>

/* general literanger headers */
#include "enum_types.h"
#include "globals.h"
#include "utility.h" // toggle_print
#include "utility_interrupt.h" // interruptor
/* required literanger class declarations */
#include "Data.decl.h"
#include "TreeBase.decl.h"
#include "TreeParameters.h"


namespace literanger {

/** Abstract base of a random forest interface. */
struct ForestBase {

    public:

        using dbl_vector_ptr = std::shared_ptr<dbl_vector>;
        using name_vector = std::vector<std::string>;
        using key_vector_ptr = std::shared_ptr<key_vector>;

        /** Non-copyable. @param[in] rhs right-hand side of copy. */
        ForestBase(const ForestBase & x) = delete;
        /** Non-assignable. @param[in] rhs right-hand side of assignment. */
        ForestBase & operator=(const ForestBase & x) = delete;
        /** Virtual destructor for pure-abstract class. */
        virtual ~ForestBase() = default;

        /** Seed the pseudo-random number generator engine.
         * @param[in] seed Value to seed ForestBase::gen with. */
        void seed_gen(const size_t seed);

        /** Plant and grow (train) trees in a random forest using supplied data.
         * @param[in] data Data to train forest with, see literanger::Data class
         * for further details about format.
         * @param[in] case_weights The weight for each case (row) in training.
         * @param[in] seed The seed for the pseudo-random number generator
         * engine.
         * @param[in] n_thread Number of threads to use when growing and
         * predicting.
         * @param[in] compute_oob_error Indicator of whether to estimate the
         * out-of-bag error or not.
         * @param[in] user_interrupt An operator that checks for user interrupt.
         * @param[out] oob_error The value of the out-of-bag error if requested.
         * @param[out] print_out A toggle-able printer for outputting progress
         * when training or predicting. */
        virtual void plant(const std::shared_ptr<const Data> data,
                           const dbl_vector_ptr case_weights,
                           const size_t seed,
                           const size_t n_thread,
                           const bool compute_oob_error,
                           const interruptor & user_interrupt,
                           double & oob_error,
                           toggle_print & print_out) = 0;

        const std::vector<TreeParameters> & get_tree_parameters() const;


    protected:

        /** Construct a random forest object.
         * @param[in] tree_type The type of tree (classification or regression)
         * to grow in the random forest.
         * @param[in] tree_parameters ...
         * @param[in] save_memory Indicator whether to aggressively release
         * memory and omit building an index (which takes up memory but speeds
         * up training). */
        ForestBase(const TreeType tree_type,
                   const std::vector<TreeParameters> tree_parameters,
                   const bool save_memory);

        /** Show the proportion of completed events in a phase.
         * @param[in] operation A suffix string that describes the current
         * process (e.g. "Growing trees").
         * @param[in] max_events The total number of events in the process.
         * @param[in] n_thread Number of threads to use when growing and
         * predicting.
         * @param[in] user_interrupt An operator that checks for user interrupt.
         * @param[out] print_out A toggle-able printer for outputting progress
         * when training or predicting.
         */
        void show_progress(std::string operation, const size_t max_events,
                           const size_t n_thread,
                           const interruptor & user_interrupt,
                           toggle_print & print_out);

        /** The type of tree grown in the forest.
         *
         * FIXME: This is used in place of RTTI to dynamically case for
         * prediction. */
        const TreeType tree_type;

        /** Number of trees in forest. */
        const size_t n_tree;

        // TODO: allow non-const???
        /** The (generic) parameters for each tree in the forest. */
        const std::vector<TreeParameters> tree_parameters;

        /** Aggressively release resources and do not construct predictor
         * indices. */
        const bool save_memory;

        /** Pseudo-random number generator for bootstrapping and also for
         * seeding each tree's pseudo-rng during the growth phase. */
        std::mt19937_64 gen;

        /** Count of the completed events in a 'queue', e.g. the number of
         * trees currently grown. */
        size_t event_count;

        /** Indicator of whether a 'queue' has been interrupted. */
        bool interrupted;

        /** Mutex for updating event_count or interrupted members */
        std::mutex mutex;

        /** Condition variable for the progress report loop. */
        std::condition_variable condition_variable;

        /** Intervals of work to perform in each thread. */
        count_vector work_intervals;

        /** A container for the trees in the forest. */
        std::vector<std::unique_ptr<TreeBase>> trees;


};


/** Make a unique ForestBase resource by forwarding arguments.
 * @param[in] args Arguments forwarded to a random forest constructor
 * @returns A unique pointer to the constructed random forest.
 * @tparam T The derived random forest type to construct.
 * @tparam ArgsT The argument types of a constructor for the derived type.
 */
template <typename T, typename... ArgsT>
std::unique_ptr<ForestBase> make_forest(ArgsT &&... args);


} /* namespace literanger */


#endif /* LITERANGER_FOREST_BASE_DECL_H */