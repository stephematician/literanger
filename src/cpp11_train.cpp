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

/* call declaration */
#include "cpp11_train.decl.h"

/* standard library headers */
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <thread>
#include <unordered_map>

/* general literanger headers */
#include "enum_types.h"
#include "utility.h"
/* literanger R package headers */
#include "cpp11_utility.h"

/* required literanger class headers */
#include "DataR.h"
#include "DataSparse.h"
#include "ForestClassification.h"
#include "ForestRegression.h"


[[cpp11::register]]
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
    cpp11::doubles response_weights,
    const size_t n_random_split, const double alpha, const double min_prop,
    const size_t seed, const bool save_memory, const size_t n_thread,
    const bool verbose
) {

    using namespace literanger;
    using namespace cpp11::literals;

    cpp11::writable::list result;

    std::unique_ptr<ForestBase> forest{ };
    std::shared_ptr<Data> data { };

    toggle_print print_out { verbose, Rprintf };
    R_user_interruptor user_interrupt { };

    const TreeType enum_tree_type(as_tree_type(tree_type));

  /* Convert the parameters for the forest to standard library types and set
   * default values. */
    const auto std_predictor_names = as_vector<std::string>(predictor_names);
    const auto std_names_of_unordered = as_vector<std::string>(
        names_of_unordered
    );
    const auto n_predictor = std_predictor_names.size();

    const auto std_sample_fraction = as_vector_ptr<double>(
        sample_fraction
    );

    set_n_try(n_try, predictor_names);
    const auto std_names_of_always_draw = as_vector<std::string>(
        names_of_always_draw
    );
    const auto std_draw_predictor_weights =
        as_nested_ptr<double,cpp11::doubles>(draw_predictor_weights);

    const SplitRule enum_split_rule(as_split_rule(split_rule));
    double min_metric_decrease;
    set_min_metric_decrease(min_metric_decrease, enum_split_rule, alpha);

    set_min_split_n_sample(min_split_n_sample, enum_tree_type);
    set_min_leaf_n_sample(min_leaf_n_sample, enum_tree_type);

  /* Construct the container for the parameters of each tree in the forest. */
    std::vector<TreeParameters> tree_parameters;
    const auto is_ordered = make_is_ordered(std_predictor_names,
                                            std_names_of_unordered);
    const auto draw_always_predictor_keys = make_draw_always_predictor_keys(
        std_predictor_names, std_names_of_always_draw, n_try
    );

    { const auto empty = std::shared_ptr<dbl_vector>(new dbl_vector());

        for (size_t j = 0; j != n_tree; ++j) {
            std::shared_ptr<dbl_vector> draw_predictor_weights_j;
            switch (std_draw_predictor_weights.size()) {
            case 0: { draw_predictor_weights_j = empty;
            } break;
            case 1: { draw_predictor_weights_j = std_draw_predictor_weights[0];
            } break;
            default: { draw_predictor_weights_j = std_draw_predictor_weights[j];
            }
            }
            set_draw_predictor_weights(
                draw_predictor_weights_j, n_predictor, n_try,
                *draw_always_predictor_keys
            );
            tree_parameters.emplace_back(
                n_predictor,
                is_ordered,
                replace, std_sample_fraction,
                n_try, draw_always_predictor_keys, draw_predictor_weights_j,
                enum_split_rule, min_metric_decrease, max_depth,
                min_split_n_sample, min_leaf_n_sample, n_random_split
            );
        }
    }


  /* Construct the data used for training */
    const bool use_sparse = sparse_x != R_NilValue;

    if (use_sparse) {

        cpp11::integers sp_Dim  = { sparse_x.attr("Dim") };

        if ((size_t)sp_Dim[1] != (size_t)predictor_names.size())
            throw std::domain_error("Mismatch between length of "
                "'predictor_names' and 'x'.");
    } else {
        if ((size_t)x.ncol() != (size_t)predictor_names.size())
            throw std::domain_error("Mismatch between length of "
                "'predictor_names' and 'x'.");
    }

    if (use_sparse) {
        data = std::shared_ptr<Data>(
            new DataSparse(cpp11::as_integers({ sparse_x.attr("Dim")}),
                           cpp11::as_integers({ sparse_x.attr("i")}),
                           cpp11::as_integers({ sparse_x.attr("p")}),
                           cpp11::as_doubles({ sparse_x.attr("x")}),
                           y)
        );
    } else {
        data = std::shared_ptr<Data>(new DataR(x, y));
    }


  /* Create the random forest object. */
    using dbl_vector_ptr = std::shared_ptr<dbl_vector>;
    dbl_vector_ptr response_values { };

    switch (enum_tree_type) {
    case TREE_CLASSIFICATION: {
        response_values = dbl_vector_ptr(
            new dbl_vector(data->get_response_values())
        );
        forest = make_forest<ForestClassification>(
            response_values,
            as_vector_ptr<double>(response_weights),
            tree_parameters, save_memory
        );
    } break;
    case TREE_REGRESSION: {
        forest = make_forest<ForestRegression>(
            min_prop, tree_parameters, save_memory
        );
    } break;
    default: throw std::invalid_argument("Unsupported tree type.");
    }


  /* Now train the forest */
    const size_t plant_n_thread = n_thread == DEFAULT_N_THREAD ?
        std::thread::hardware_concurrency() : n_thread;
    if (plant_n_thread == 0)
        throw std::domain_error("'n_thread' must be positive.");

    double oob_error;
    forest->plant(
        data, as_vector_ptr<double>(case_weights), seed, plant_n_thread,
        true, user_interrupt, oob_error, print_out
    );
    // TODO: per-tree case weights?


  /* Store the results (selected arguments) */
    result.push_back({ "predictor_names"_nm = predictor_names });
    result.push_back({ "names_of_unordered"_nm = names_of_unordered });
    result.push_back({ "tree_type"_nm = tree_type });
    result.push_back({ "n_tree"_nm = n_tree });
    result.push_back({ "n_try"_nm = n_try });
    result.push_back({ "split_rule"_nm = split_rule });
    result.push_back({ "max_depth"_nm = max_depth });
    result.push_back({ "min_metric_decrease"_nm = min_metric_decrease });
    result.push_back({ "min_split_n_sample"_nm = min_split_n_sample });
    result.push_back({ "min_leaf_n_sample"_nm = min_leaf_n_sample });
    // TODO:  min_prop ?
    result.push_back({ "seed"_nm = seed });
    result.push_back({ "oob_error"_nm = oob_error });
    // TODO: ??? confusion matrix

    if (enum_split_rule == EXTRATREES)
        result.push_back({ "n_random_split"_nm = n_random_split });

    if (enum_tree_type == TREE_CLASSIFICATION) {
        result.push_back({ "response_values"_nm = *response_values });
    }

    result.push_back({
        "cpp11_ptr"_nm = cpp11::external_pointer<ForestBase>(forest.release())
    });

    return result;

}


/* Helpers to set default values. */

void set_n_try(size_t & n_try, cpp11::strings predictor_names) {
    if (n_try != 0) return;
    n_try = (size_t)std::max(1., std::sqrt((double)(predictor_names.size())));
}


void set_min_split_n_sample(size_t & min_split_n_sample,
                            const literanger::TreeType tree_type) {
    using namespace literanger;
    #if !defined(__GNUC__) || __GNUC__ >= 5
      using umap_key_t = TreeType;
    #else
      using umap_key_t = size_t;
    #endif

    if (min_split_n_sample != 0) return;

    static std::unordered_map<umap_key_t,size_t> table = {
        { TreeType::TREE_CLASSIFICATION,
          DEFAULT_MIN_SPLIT_N_SAMPLE_CLASSIFICATION },
        { TreeType::TREE_REGRESSION, DEFAULT_MIN_SPLIT_N_SAMPLE_REGRESSION }
    };

    min_split_n_sample = table[tree_type];
}


void set_min_leaf_n_sample(size_t & min_leaf_n_sample,
                           const literanger::TreeType tree_type) {
    using namespace literanger;
    #if !defined(__GNUC__) || __GNUC__ >= 5
      using umap_key_t = TreeType;
    #else
      using umap_key_t = size_t;
    #endif

    if (min_leaf_n_sample != 0) return;

    static std::unordered_map<umap_key_t,size_t> table = {
        { TreeType::TREE_CLASSIFICATION,
          DEFAULT_MIN_LEAF_N_SAMPLE_CLASSIFICATION },
        { TreeType::TREE_REGRESSION, DEFAULT_MIN_LEAF_N_SAMPLE_REGRESSION }
    };

    min_leaf_n_sample = table[tree_type];
}


void set_min_metric_decrease(double & min_metric_decrease,
                             const literanger::SplitRule split_rule,
                             const double alpha) {
    using namespace literanger;

    switch (split_rule) {
    case EXTRATREES: case LOGRANK: case HELLINGER: {
        min_metric_decrease = 0;
    } break;
    case BETA: {
        min_metric_decrease = -std::numeric_limits<double>::max();
    } break;
    case MAXSTAT: {
        min_metric_decrease = -alpha;
    } break;
    default: throw std::runtime_error("Unexpected value of split metric.");
    }

}

void set_draw_predictor_weights(
    std::shared_ptr<std::vector<double>> draw_predictor_weights,
    const size_t n_predictor, const size_t n_try,
    const std::vector<size_t> & draw_always_predictor_keys
) {

    if (draw_predictor_weights->empty()) return;

    if (draw_predictor_weights->size() != n_predictor)
        throw std::invalid_argument("Number of draw-predictor weights not "
            "equal to number of predictors.");

  /* indicator variable for belonging to always-draw predictor */
    std::vector<bool> is_always(n_predictor, false);
    for (auto key : draw_always_predictor_keys) is_always[key] = true;

    size_t n_zero_weight = 0;

    for (size_t j = 0; j != n_predictor; ++j) {
        double & w = (*draw_predictor_weights)[j];
        if (w < 0)
            throw std::domain_error("One or more draw-predictor weights not "
                "in range [0,Inf).");
        w = w != 0 && !is_always[j] ? w : (++n_zero_weight, 0);
    }

    if (n_predictor - n_zero_weight < n_try)
        throw std::invalid_argument("Too many zeros in draw-predictor weights. "
            "Need at least n_try variables to split at.");

}

