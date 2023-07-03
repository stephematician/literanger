/*-------------------------------------------------------------------------------
 * This file is part of 'literanger'. literanger was adapted from the 'ranger'
 * package for R Statistical Software <https://www.r-project.org>. ranger was
 * authored by Marvin N Wright with the GNU General Public License version 3.
 * The adaptation was performed by Stephen Wade in 2023. literanger carries the
 * same license, terms, and permissions as ranger.
 *
 * literanger is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * literanger is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with literanger. If not, see <http://www.gnu.org/licenses/>.
 *
 * Written by:
 *
 *   Stephen Wade
 *   Cancer Council New South Wales
 *   Woolloomooloo NSW 2011
 *   Australia
 *-------------------------------------------------------------------------------
 */
#ifndef LITERANGER_CPP11_TRAIN_DECL_H
#define LITERANGER_CPP11_TRAIN_DECL_H

/* standard library headers */
#include <cstddef>
#include <memory>
#include <string>
#include <vector>

/* cpp11 and R headers */
#include "cpp11.hpp"

/* general literanger headers */
#include "enum_types.h"

/** Train a random forest.
 *
 * See file R/train.R in the R package for further details.
 *
 * @param[in] x The predictor data where each column is a predictor and each row
 * is an observation (or case).
 * @param[in] y The response data, single column, where each row is an
 * observation, may be expanded to include survival data (with censoring) in
 * future.
 * @param[in] sparse_x Optional, NULL to ignore; the predictor
 * data with the same layout as in \p x, but more efficient (memory-wise) for
 * matrices with many zero elements.
 * @param[in] case_weights Optional. The weights for sampling training,
 * observations with larger weights will be selected with higher probability
 * during bootstrap sampling for each tree.
 * @param[in] tree_type The type of tree in the forest - currently either
 * "classification" or "regression".
 * @param[in] n_tree The number of trees to train in the forest.
 * @param[in] predictor_names The names of predictor (independent) variables,
 * must have the same length as the number of columns in @p x or @p sparse_x.
 * @param[in] names_of_unordered The names of the predictors that are unorderded
 * and hence values for splitting are taken from all possible partititons.
 * @param[in] replace Indicator for whether to draw with (or without)
 * replacement when generating the sample for training each tree.
 * @param[in] sample_fraction If scalar-valued; the fraction of data to sample
 * from when training each tree. If vector-valued, for classification trees
 * only, the fraction to sample by the response category (class).
 * @param[in] n_try The number of predictors to randomly draw as candidates when
 * splitting a node.
 * @param[in] draw_predictor_weights For predictor-drawing weights shared by all
 * trees; a vector of _non-negative_ weights for each predictor. For
 * tree-specific predictor-drawing weights; a list of size @p n_tree containing
 * (non-negative) vectors with length equal to the number of predictors.
 * @param[in] names_of_always_draw Names of predictor (variables) to be selected
 * _in addition_ to the @p n_try predictors drawn as candidates to split by.
 * @param[in] split_rule The splitting rule used to select and evaluate
 * candidate values (of a predictor) when splitting, either "gini", "variance",
 * "extratrees", "hellinger", "maxstat", or "beta" (valid values depend upon
 * tree type).
 * @param[in] max_depth Maximum depth of a tree, i.e. one implies that one split
 * is allowed, zero implies unlimited depth.
 * @param[in] min_split_n_sample Smallest number of in-bag samples a node must
 * have to be split.
 * @param[in] min_leaf_n_sample Minimum number of in-bag samples in a leaf node.
 * @param[in] response_weights "lassification" @p tree_type only; weights for
 * the response classes (in order they appear _in data_) when evaluating node
 * metrics, and used by each tree when tallying votes.
 * @param[in] n_random_split "extratrees" @p split_rule only;the number of
 * random candidate-values to draw when splitting a node.
 * @param[in] alpha "maxstat" @p split_rule only: significance threshold to
 * allow splitting.
 * @param[in] min_prop "maxstat" @p split_rule only: lower quantile of covariate
 * distribution to be considered for splitting.
 * @param[in] seed Seed used for (seeding) pseudo-random number generators for
 * each tree.
 * @param[in] save_memory Indicator for aggressively releasing memory and for
 * skipping building an index for the predictors (used to speed up training).
 * @param[in] n_thread The number of threads to use for training, if zero then
 * the number of threads will be the value returned by
 * `std::thread::hardware_concurrency`.
 * @param[in] verbose Indicator for additional printed output while growing and
 * predicting.
 *
 * @returns A list (in R) with:
 * -  `predictor_names`: the names of the predictor variables in the model.
 * -  `names_of_unordered`: the names of predictors that are unordered.
 * -  `tree_type`: the type of tree in the forest.
 * -  `n_tree`: the number of trees that were trained.
 * -  `n_try`: the number of predictors drawn as candidates for each split.
 * -  `split_rule`: the name of the split metric used.
 * -  `max_depth`: the maximum allowed depth of a tree in the forest.
 * -  `min_metric_decrease`: the minimum decrease in the metric for an
 *     acceptable split, equal to negative `alpha` for maximally selected rank
 *     statistic split metric, else zero.
 * -   `min_split_n_sample`: the minimum number of in-bag samples in a node
 *     prior to splitting.
 * -   `min_leaf_n_sample`: the minumum number of in-bag samples in a leaf node.
 * -   `seed`: the seed supplied to the C++ library.
 * -   `oob_error`: the misclassification rate or the mean square error using
 *     out-of-bag samples.}
 * -   `n_random_split`: "extratrees" split rule only; the number of candidate
 *     _values_ drawn when splitting.
 * -   `response_values`: classification only; the values of the response in
 *     the order they appear in the data.
 * -   `cpp11_ptr`: an external pointer to the trained forest.
 */
cpp11::list cpp11_train(
    cpp11::doubles_matrix<> x, cpp11::doubles_matrix<> y, cpp11::sexp sparse_x,
    cpp11::doubles case_weights,
    std::string tree_type, const size_t n_tree,
    cpp11::strings predictor_names, cpp11::strings names_of_unordered,
    const bool replace, cpp11::doubles sample_fraction,
    size_t n_try,
    cpp11::list draw_predictor_weights, cpp11::strings names_of_always_draw,
    std::string split_rule, const size_t max_depth, size_t min_split_n_sample,
    size_t min_leaf_n_sample,
  /* tree-type specific */
    cpp11::doubles response_weights,
  /* split-rule specific */
    const size_t n_random_split, const double alpha, const double min_prop,
    const size_t seed, const bool save_memory, const size_t n_thread,
    const bool verbose
);

/* Helpers to set default values. */
void set_n_try(size_t & n_try, cpp11::strings predictor_names);
void set_min_split_n_sample(size_t & min_split_n_sample,
                            const literanger::TreeType tree_type);
void set_min_leaf_n_sample(size_t & min_leaf_n_sample,
                           const literanger::TreeType tree_type);
void set_min_metric_decrease(double & min_metric_decrease,
                             const literanger::SplitRule split_rule,
                             const double alpha);

/** Update the (shared resource) container of the weights used when drawing
 * predictors for splitting.
 * @param[in,out] draw_predictor_weights Weights for each predictor when drawing
 * candidates.
 * @param[in] n_predictor The number of predictors in the random forest
 * model.
 * @param[in] n_try The number of candidate predictors for each split.
 * @param[in] draw_always_predictor_keys A container of keys for predictors
 * that are always candidates for splitting. */
void set_draw_predictor_weights(
    std::shared_ptr<std::vector<double>> draw_predictor_weights,
    const size_t n_predictor, const size_t n_try,
    const std::vector<size_t> & draw_always_predictor_keys
);


#endif /* LITERANGER_CPP11_TRAIN_DECL_H */

