#ifndef LITERANGER_CPP11_TRAIN_DECL_H
#define LITERANGER_CPP11_TRAIN_DECL_H

/* standard library headers */
#include <cstddef>
#include <string>
#include <vector>

/* cpp11 and R headers */
#include "cpp11.hpp"

/* generic literanger headers */
#include "enum_types.h"


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

