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
#ifndef LITERANGER_TREE_PARAMETERS_H
#define LITERANGER_TREE_PARAMETERS_H

#include <cstddef>
#include <memory>
#include <string>
#include <vector>

/* general literanger headers */
#include "enum_types.h"
#include "globals.h"


namespace literanger {

/** Generic parameters for a tree in a random forest.
 *
 * Parameters that describe the sampling, drawing, and splitting of a tree from
 * a random forest. A vector of these parameters is passed to the forest
 * constructor which dictates how many trees and what the values of the
 * parameters for each tree are. */
struct TreeParameters {

    using name_vector = std::vector<std::string>;
    using key_vector_ptr = std::shared_ptr<key_vector>;
    using dbl_vector_ptr = std::shared_ptr<dbl_vector>;

    /** Generic tree parameter constructor
     * @param[in] n_predictor Number of predictors in the tree model.
     * @param[in] is_unordered Container of indicators for ordered predictors.
     * @param[in] replace Whether to sample with replacement when training.
     * @param[in] sample_fraction The fraction of samples to use to train
     * each tree.
     * @param[in] n_try The number of candidate predictors for each split.
     * @param[in] draw_always_predictor_keys The key of each predictor that will
     * always be a candidate for splitting.
     * @param[in] draw_predictor_weights Weights for each predictor when drawing
     * candidates.
     * @param[in] split_rule The rule for identifying the best split.
     * @param[in] min_metric_decrease The minimum change in the metric for an
     * acceptable split.
     * @param[in] max_depth The maximum depth of any tree.
     * @param[in] min_split_n_sample The minimum number of in-bag samples in a
     * node that may be split in the growth phase.
     * @param[in] min_leaf_n_sample The minimum number of in-bag samples in a
     * leaf node in the growth phase.
     * @param[in] n_random_split The number of values to draw when splitting
     * via the extratrees rule. */
    TreeParameters(
        const size_t n_predictor,
        const std::shared_ptr<std::vector<bool>> is_ordered,
        const bool replace, const dbl_vector_ptr sample_fraction,
        const size_t n_try, const key_vector_ptr draw_always_predictor_keys,
        const dbl_vector_ptr draw_predictor_weights,
        const SplitRule split_rule, const double min_metric_decrease,
        const size_t max_depth, const size_t min_split_n_sample,
        const size_t min_leaf_n_sample, const size_t n_random_split
    );

    /** @name Specification of predictors in tree or forest. */
    /**@{*/
    /** Number of predictors in the random forest model. */
    size_t n_predictor;
    /** Indicator for ordered predictors */
    std::shared_ptr<std::vector<bool>> is_ordered;
    /**@}*/

    /** @name Resampling training data for growing (training) a tree. */
    /**@{*/
    /** Indicator for sampling with replacement when when fitting. */
    bool replace;
    /** The fraction of samples to use when fitting each tree (scalar) or, when
     * when a vector is supplied, the response-specific fractions. */
    dbl_vector_ptr sample_fraction;
    /**@}*/

    /** @name Drawing candidate predictors for node splitting. */
    /**@{*/
    /** Number of randomly-drawn predictors amongst the candidates at each
      * node split */
    size_t n_try;
    /** Predictors that are always candidates for splitting. */
    key_vector_ptr draw_always_predictor_keys;
    /** Weights for each predictor that determine probability of selection as a
     * candidate for splitting (see std::discrete_distribution). */
    dbl_vector_ptr draw_predictor_weights;
    /**@}*/

    /** @name Node-splitting rules. */
    /**@{*/
    /** Rule for selecting the predictor and value to split on. */
    SplitRule split_rule;
    /** Minimum decrease in metric that will be accepted when splitting. */
    double min_metric_decrease;
    /** Maximum depth of the trees in the forest. */
    size_t max_depth;
    /** Minimum number of in-bag samples a node must have to consider for
     * splitting. */
    size_t min_split_n_sample;
    /** Minimum number of in-bag samples in a leaf node. */
    size_t min_leaf_n_sample;
    /** Number of random splits to draw when using extra-random trees
      * algorithm. */
    size_t n_random_split;
    /**@}*/


};


inline TreeParameters::TreeParameters(
    const size_t n_predictor,
    const std::shared_ptr<std::vector<bool>> is_ordered,
    const bool replace, const dbl_vector_ptr sample_fraction,
    const size_t n_try, const key_vector_ptr draw_always_predictor_keys,
    const dbl_vector_ptr draw_predictor_weights,
    const SplitRule split_rule, const double min_metric_decrease,
    const size_t max_depth, const size_t min_split_n_sample,
    const size_t min_leaf_n_sample, const size_t n_random_split
) :
    n_predictor(n_predictor), is_ordered(is_ordered),
    replace(replace), sample_fraction(sample_fraction),
    n_try(n_try), draw_always_predictor_keys(draw_always_predictor_keys),
    draw_predictor_weights(draw_predictor_weights),
    split_rule(split_rule), min_metric_decrease(min_metric_decrease),
    max_depth(max_depth), min_split_n_sample(min_split_n_sample),
    min_leaf_n_sample(min_leaf_n_sample), n_random_split(n_random_split)
{
    if (this->n_try == 0) throw std::domain_error("'n_try' must be positive.");
    if (this->split_rule == EXTRATREES && this->n_random_split == 0)
        throw std::domain_error("'n_random_split' must be positive.");
    if (this->n_try > this->n_predictor)
        throw std::domain_error("'n_try' can not be larger than number of "
            "predictors (columns).");
}


} /* namespace literanger */


#endif /* LITERANGER_TREE_PARAMETERS_H */

