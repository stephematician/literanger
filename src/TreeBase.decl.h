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
#ifndef LITERANGER_TREE_BASE_DECL_H
#define LITERANGER_TREE_BASE_DECL_H

/* standard library headers */
#include <cstddef>
#include <memory>
#include <random>
#include <utility>
#include <vector>

/* general literanger headers */
#include "enum_types.h"
#include "globals.h"
/* required literanger class declarations */
#include "Data.decl.h"
#include "TreeParameters.h"


namespace literanger {

/** Abstract base of a tree-interface
 *
 * TODO: better documentation
 */
struct TreeBase {

    public:

        using key_vector_ptr = std::shared_ptr<key_vector>;
        using dbl_vector_ptr = std::shared_ptr<dbl_vector>;

        /** Non-copyable. @param[in] rhs right-hand side of copy. */
        TreeBase(const TreeBase & x) = delete;
        /** Non-assignable. @param[in] rhs right-hand side of assignment. */
        TreeBase & operator=(const TreeBase & x) = delete;
        /** Virtual destructor for pure-abstract class. */
        virtual ~TreeBase() = default;

        /** Seed the pseudo-random number generator engine.
         * @param[in] seed Value to seed TreeBase::gen with. */
        void seed_gen(const size_t seed);

        /** @name Access node data. */
        /**@{*/
        const key_vector & get_split_keys() const;
        const dbl_vector & get_split_values() const;
        const key_vector & get_left_children() const;
        const key_vector & get_right_children() const;
        /**@}*/

        /** Grow (train) the tree using supplied data.
         * @param[in] data Data to train forest with, see literanger::Data class
         * for further details about format.
         * @param[in] case_weights The weight for each case (row) in training.
         * @param[in] compute_oob_error Indicator of whether to return the
         * out-of-bag keys or not.
         * @returns A vector of out-of-bags keys (empty if not requested). */
        key_vector grow(const std::shared_ptr<const Data> data,
                        const dbl_vector_ptr case_weights,
                        const bool compute_oob_error);

        /** Get the number of samples contained in a node.
         * @param[in] node_key The node to query.
         * @returns The number of samples in the node. */
        size_t get_n_sample_node(const size_t node_key) const;


    protected:

        /** Construct a tree object.
         * @param[in] parameters Parameters that describe the sampling,
         * drawing, and splitting for trees in a random forest.
         * @param[in] save_memory Indicator whether to aggressively release
         * memory and omit building an index (which takes up memory but speeds
         * up training). */
        TreeBase(const TreeParameters & parameters, const bool save_memory);

        /** @name Generic tree parameters.
         * Parameters that describe the sampling, drawing, and splitting for
         * trees in a random forest.
         * @see literanger::TreeParameters */
        /*@{*/
        /** @copydoc TreeParameters::n_predictor */
        const size_t n_predictor;
        /** @copydoc TreeParameters::is_ordered */
        const std::shared_ptr<const std::vector<bool>> is_ordered;
        /** @copydoc TreeParameters::replace */
        const bool replace;
        /** @copydoc TreeParameters::sample_fraction */
        const dbl_vector_ptr sample_fraction;
        /** @copydoc TreeParameters::n_try */
        const size_t n_try;
        /** @copydoc TreeParameters::draw_always_predictor_keys */
        const key_vector_ptr draw_always_predictor_keys;
        /** @copydoc TreeParameters::draw_predictor_weights */
        const dbl_vector_ptr draw_predictor_weights;
        /** @copydoc TreeParameters::split_rule */
        const SplitRule split_rule;
        /** @copydoc TreeParameters::min_metric_decrease */
        const double min_metric_decrease;
        /** @copydoc TreeParameters::max_depth */
        const size_t max_depth;
        /** @copydoc TreeParameters::min_split_n_sample */
        const size_t min_split_n_sample;
        /** @copydoc TreeParameters::min_leaf_n_sample */
        const size_t min_leaf_n_sample;
        /** @copydoc TreeParameters::n_random_split */
        const size_t n_random_split;
        /**@}*/

        /** Aggressively release resources and use a unique value mapping. */
        const bool save_memory;

        /** Pseudo-random number generator for all sampling. */
        std::mt19937_64 gen;

        /** The predictor key for each node that identifies the variable to
         * split by. */
        key_vector split_keys;

        /** The value for each node that determines whether a data point belongs
         * in the left or right child (given the predictor). */
        dbl_vector split_values;

        /** A pair of containers for left and right child-node keys */
        std::pair<key_vector,key_vector> child_node_keys;

        /** Reference to the left child-node keys */
        key_vector & left_children = child_node_keys.first;

        /** Reference to the left child-node keys */
        key_vector & right_children = child_node_keys.second;

        /** The starting offset of the observations within a container of
         * partially-sorted observation keys for each node. */
        count_vector start_pos;

        /** The past-the-end offset of the observations within a container of
         * partially-sorted observation keys for each node. */
        count_vector end_pos;

        /** Count of the number of observations for each candidate split
         * value. */
        count_vector node_n_by_candidate;


    private:

        /** */
        void push_back_empty_node();

        /**
         * Bootstrap/draw a sample from the set of keys `[0, 1, 2, ..., N-1]`
         * and optionally return the values _not_ drawn.
         *
         * @param[in] n_sample Both the number of keys to (randomly) draw and
         * the size of the set of keys.
         * @param[in] get_oob_keys Indicator for returning the out-of-bag
         * keys.
         * @param[out] sample_keys The randomly-drawn keys from the set.
         * @param[out] oob_keys The 'out-of-bag' keys - i.e. the keys that
         * aren't in `sample_keys`.
         */
        void resample_unweighted(const size_t n_sample, const bool get_oob_keys,
                                 key_vector & sample_keys,
                                 key_vector & oob_keys);

        /**
         * Boostrap/draw a sample from a set of keys `[0, 1, 2, ..., N-1]` where
         * each key has a user-provided probability of selection, and optionally
         * return the values _not_ drawn.
         *
         * @param[in] n_sample Both the number of keys to (randomly) draw and
         * the size of the set of keys.
         * @param[in] weights The weights or probabilities for each key.
         * @param[in] get_oob_keys Indicator for returning the out-of-bag
         * keys.
         * @param[out] sample_keys The randomly-drawn keys from the set.
         * @param[out] oob_keys The 'out-of-bag' keys - i.e. the keys that
         * aren't in `sample_keys`.
         */
        void resample_weighted(
            const size_t n_sample, const dbl_vector_ptr weights,
            const bool get_oob_keys,
            key_vector & sample_keys, key_vector & oob_keys
        );

        /**
         * Bootsrap/draw a sample from the set of keys `[0, 1, 2, ..., N-1]`
         * with a user-specified fraction for each response value, and
         * optionally return the values _not_ drawn.
         *
         * @param[in] data Data to use for growth (training); the number of rows
         * is used to identify the size of the set to draw from.
         * @param[in] get_oob_keys Indicator for returning the out-of-bag
         * keys.
         * @param[out] sample_keys The randomly-drawn keys from the set.
         * @param[out] oob_keys The 'out-of-bag' keys - i.e. the keys that
         * aren't in `sample_keys`.
         */
        void resample_response_wise(
            const std::shared_ptr<const Data> data,
            const bool get_oob_keys,
            key_vector & sample_keys, key_vector & oob_keys
        );

        /**
         * Base-class implementation of response-wise resampling does nothing.
         *
         * @param[in] data Data used for growth (training).
         * @param[out] sample_keys A container of randomly drawn keys (i.e.
         * row-offsets) of observations.
         * @param[out] inbag_counts A container of counts of the number of
         * timees each observation appears in-bag.
         */
        virtual void resample_response_wise_impl(
            const std::shared_ptr<const Data> data,
            key_vector & sample_keys,
            count_vector & inbag_counts
        );

        /**
         * Split a node using rules for selecting candidate predictors,
         * evaluating decrease in impurity, and selecting candidate values to
         * split by.
         *
         * @param[in] node_key The key (offset of split_vars vector etc) of the
         * current node.
         * @param[in] depth The current depth of the tree.
         * @param[in] last_left_node_key The most-recently generated left node
         * at the current depth.
         * @param[in] data Data to used for growth (training).
         * @param[out] sample_keys The partially-ordered keys where any key to
         * the right of another key is placed later in the container.
         * @returns Indicator for whether a split was performed: this is the
         * opposite of original ranger.
         */
        bool split_node(const size_t node_key,
                        const size_t depth,
                        const size_t last_left_node_key,
                        const std::shared_ptr<const Data> data,
                        key_vector & sample_keys);

        /** Draw candidate predictors for splitting.
         * @returns A vector of predictor keys (column offsets) that are
         * candidates for splitting.
         */
        key_vector draw_candidates();

        /** Prepare a tree for growth by reserving space for terminal nodes.
         *
         * @param[in] data Data to grow (or train) tree with. Contains
         * observations of predictors and the response, the former has
         * predictors across columns and observations by row, and the latter is
         * usually a column vector (or matrix).
         */
        virtual void new_growth(const std::shared_ptr<const Data> data) = 0;

        /** Default finalisation (do nothing) for growth phase.
         *
         * Implementation may use this to do any post-processing for terminal
         * nodes.
         */
        virtual void finalise_growth();

        /** */
        virtual void push_back_empty_node_impl();

        /** Store the observed values in the leaf (terminal) node container.
         *
         * @param[in] node_key The key for a new leaf (terminal) node.
         * @param[in] data Data to train forest with. Contains observations of
         * predictors and the response, the former has predictors across
         * columns and observations by row, and the latter is usually a column
         * vector (or matrix).
         * @param[in] sample_keys Container of partially ordered observation
         * keys (row offsets) used to grow the tree; any node that is left of
         * another is found later in the container.
         */
        virtual void add_terminal_node(const size_t node_key,
                                       const std::shared_ptr<const Data> data,
                                       const key_vector & sample_keys) = 0;

        /** Compare two responses for equality.
         * @param[in] data Data for a random forest (prediction or training).
         * Contains observations of predictors and the response, the latter is
         * usually a column vector (or matrix) with one observation (or case)
         * per row.
         * @param[in] lhs_key The row-offset of the left-hand-side of the
         * comparison.
         * @param[in] rhs_key THe row-offset of the right-hand-side of the
         * comparison.
         * @returns True if the response values are numerically equal.
        */
        virtual bool compare_response(
            const std::shared_ptr<const Data> data,
            const size_t lhs_key, const size_t rhs_key
        ) const = 0;

        /** Add the best-performing split for a specified node; if no split
         * decreases impurity then do nothing.
         * @param[in] node_key The node to evaluate.
         * @param[in] data Data to use for growth (training).
         * @param[in] sample_keys The partially-sorted keys in the sample for
         * this tree.
         * @param[in] split_candidate_keys Identifies the predictors that are
         * candidates for splitting.
         * @returns Whether a split was added.
         */
        virtual bool push_best_split(
            const size_t node_key, const std::shared_ptr<const Data> data,
            const key_vector & sample_keys,
            const key_vector & split_candidate_keys
        ) = 0;


};


} /* namespace literanger */


#endif /* LITERANGER_TREE_BASE_DECL_H */

