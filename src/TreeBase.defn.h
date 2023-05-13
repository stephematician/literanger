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
#ifndef LITERANGER_TREE_BASE_DEFN_H
#define LITERANGER_TREE_BASE_DEFN_H

/* class declaration */
#include "TreeBase.decl.h"

/* standard library headers */
#include <algorithm>
#include <iterator>
#include <stdexcept>

/* general literanger headers */
#include "utility_draw.h"
/* required literanger class definitions */
#include "Data.defn.h"


namespace literanger {

inline TreeBase::TreeBase(const TreeParameters & parameters,
                          const bool save_memory) :
    n_predictor(parameters.n_predictor), is_ordered(parameters.is_ordered),
    replace(parameters.replace), sample_fraction(parameters.sample_fraction),
    n_try(parameters.n_try),
    draw_always_predictor_keys(parameters.draw_always_predictor_keys),
    draw_predictor_weights(parameters.draw_predictor_weights),
    split_rule(parameters.split_rule),
    min_metric_decrease(parameters.min_metric_decrease),
    max_depth(parameters.max_depth),
    min_split_n_sample(parameters.min_split_n_sample),
    min_leaf_n_sample(parameters.min_leaf_n_sample),
    n_random_split(parameters.n_random_split),
    save_memory(save_memory)
{ }


inline void TreeBase::seed_gen(const size_t seed) { gen.seed(seed); }


inline const key_vector & TreeBase::get_split_keys() const {
    return split_keys;
}


inline const dbl_vector & TreeBase::get_split_values() const {
    return split_values;
}


inline const key_vector & TreeBase::get_left_children() const {
    return left_children;
}


inline const key_vector & TreeBase::get_right_children() const {
    return right_children;
}


inline key_vector TreeBase::grow(
    const std::shared_ptr<const Data> data,
    const dbl_vector_ptr case_weights,
    const bool compute_oob_error
) {

    const size_t n_sample = data->get_n_row();
    key_vector sample_keys { };
    key_vector oob_keys { };

  /* Should be starting with empty tree */
    if (split_keys.size() != 0)
        throw std::runtime_error("Expected to start with empty tree.");

  /* Implementation-specific initialisation - usually for any data used by
   * the split_node implementation */
    new_growth(data);

  /* Construct first node - will be modified by first call to split_node */
    push_back_empty_node();

    const bool response_wise = sample_fraction->size() > 1;
    const bool weighted = !case_weights->empty();
    if (weighted && response_wise)
        throw std::invalid_argument("Cannot have both weighted and "
            "response-wise (class-wise) weighting.");

  /* Depending on the sampling strategy, get the keys for the observed data
   * to be used for growing this tree; optionally also get the keys for the
   * out-of-bag data if we're computing OOB error */
    if (weighted) {
        resample_weighted(n_sample, case_weights, compute_oob_error,
                          sample_keys, oob_keys);
    } else if (response_wise) {
        resample_response_wise(data, compute_oob_error, sample_keys, oob_keys);
    } else {
        resample_unweighted(n_sample, compute_oob_error, sample_keys, oob_keys);
    }

  /* Now we iteratively split nodes */
    size_t depth = 0, last_left_node_key = 0;
    start_pos[0] = 0;
    end_pos[0] = sample_keys.size();

    for (size_t n_open_node = 1, node_key = 0; n_open_node != 0; ++node_key) {
        const bool did_split = split_node(
            node_key, depth, last_left_node_key, data, sample_keys
        );
        if (!did_split) { --n_open_node; } else {
            ++n_open_node;
            if (node_key >= last_left_node_key) {
                last_left_node_key = split_keys.size() - 2;
                ++depth;
            }
        }
    }

  /* Implementation specific finalisation */ /* NOTE: might be redundant? */
    finalise_growth();

    return oob_keys;

}


inline size_t TreeBase::get_n_sample_node(const size_t node_key) const {
    return end_pos[node_key] - start_pos[node_key];
}


inline void TreeBase::push_back_empty_node() {

    split_keys.emplace_back(0);
    split_values.emplace_back(0);
    left_children.emplace_back(0);
    right_children.emplace_back(0);
    start_pos.emplace_back(0);
    end_pos.emplace_back(0);

    push_back_empty_node_impl();

}


inline void TreeBase::resample_unweighted(const size_t n_sample,
                                          const bool get_oob_keys,
                                          key_vector & sample_keys,
                                          key_vector & oob_keys) {

    const size_t n_sample_inbag = (size_t)(n_sample * (*sample_fraction)[0]);

    sample_keys.clear();
    if (get_oob_keys) oob_keys.clear();

    if (replace) {

        count_vector inbag_counts = count_vector(n_sample, 0);
        draw_replace(n_sample_inbag, n_sample, gen, sample_keys, inbag_counts);

        if (get_oob_keys) {
          /* Reserves more than expected number of out-of-bag samples */
          /* NOTE: M.N. Wright's ranger used exp(-frac) + 0.1: our adaptation
           * uses exp(-frac + 0.15) so that no more than n_sample is reserved */
            const double fraction = (double)n_sample_inbag / (double)n_sample;
            oob_keys.reserve(n_sample * std::exp(-fraction + 0.15));
            for (size_t j = 0; j != n_sample; ++j) {
                if (inbag_counts[j] == 0) oob_keys.push_back(j);
            }
        }

    } else {

      /* FIXME: should test if this is faster than using draw_no_replace? */
        sample_keys.assign(n_sample, 0);
        std::iota(sample_keys.begin(), sample_keys.end(), 0);
        std::shuffle(sample_keys.begin(), sample_keys.end(), gen);
        if (get_oob_keys) {
            oob_keys.reserve(n_sample - n_sample_inbag);
            std::copy(sample_keys.cbegin() + n_sample_inbag,
                      sample_keys.cend(),
                      std::back_inserter(oob_keys));
        }
        sample_keys.resize(n_sample_inbag);

    }

}


inline void TreeBase::resample_weighted(const size_t n_sample,
                                        const dbl_vector_ptr weights,
                                        const bool get_oob_keys,
                                        key_vector & sample_keys,
                                        key_vector & oob_keys) {

    if (weights->size() != n_sample)
        throw std::invalid_argument("Case weights must have the same length "
            "as number of rows in data.");
    const size_t n_sample_inbag = (size_t)(n_sample * (*sample_fraction)[0]);
    count_vector inbag_counts = count_vector(n_sample, 0);

    sample_keys.clear();
    if (get_oob_keys) oob_keys.clear();

    if (replace) {
        draw_replace_weighted(n_sample_inbag, *weights, gen,
                              sample_keys, inbag_counts);

    } else
        draw_no_replace_weighted(n_sample_inbag, *weights, gen,
                                 sample_keys, inbag_counts);

    if (get_oob_keys) {

        const double fraction = (double)n_sample_inbag / (double)n_sample;
        oob_keys.reserve( /* See NOTE in resample_unweighted */
            replace ? (size_t)(n_sample * std::exp(-fraction + 0.15)) :
                      n_sample - n_sample_inbag
        );
        for (size_t j = 0; j != n_sample; ++j) {
            if (inbag_counts[j] == 0) oob_keys.emplace_back(j);
        }

    }

}


inline void TreeBase::resample_response_wise(
    const std::shared_ptr<const Data> data,
    const bool get_oob_keys,
    key_vector & sample_keys,
    key_vector & oob_keys
) {

    const size_t n_sample = data->get_n_row();
    count_vector inbag_counts = count_vector(n_sample, 0);

    sample_keys.clear();
    if (get_oob_keys) oob_keys.clear();

  /* Implementation-specific response-wise bootstrap/draw */
    resample_response_wise_impl(data, sample_keys, inbag_counts);

    const size_t n_sample_inbag = sample_keys.size();

    if (get_oob_keys) {

        const double fraction = (double)n_sample_inbag / (double)n_sample;
        oob_keys.reserve( /* See NOTE in resample_unweighted */
            replace ? (size_t)(n_sample * std::exp(-fraction + 0.15)) :
                      n_sample - n_sample_inbag
        );
        for (size_t j = 0; j != n_sample; ++j) {
            if (inbag_counts[j] == 0) oob_keys.emplace_back(j);
        }

    }

}


inline void TreeBase::resample_response_wise_impl(
    const std::shared_ptr<const Data> data,
    key_vector & sample_keys,
    count_vector & inbag_counts
) {
    throw std::invalid_argument("Response-wise sampling not supported for this "
        "tree type.");
}


inline bool TreeBase::split_node(const size_t node_key,
                                 const size_t last_left_node_key,
                                 const size_t depth,
                                 const std::shared_ptr<const Data> data,
                                 key_vector & sample_keys) {

    const size_t n_sample_node = get_n_sample_node(node_key);

    if (max_depth && depth > max_depth)
        throw std::runtime_error("Cannot split a node that is already at "
            "maximum depth of tree.");

    { /* Test if we have reached a terminal node */
        const bool too_deep = node_key >= last_left_node_key && max_depth &&
            depth == max_depth;

        if (n_sample_node < min_split_n_sample || too_deep) {
            add_terminal_node(node_key, data, sample_keys);
            return false;
        }
    }

    { /* Test if the node is 'pure' - i.e. all values equal */
        bool pure = true;
        const size_t start_key = sample_keys[start_pos[node_key]];
        for (size_t j = start_pos[node_key]; j != end_pos[node_key]; ++j) {
            const size_t test_key = sample_keys[j];
            if (!compare_response(data, start_key, test_key)) {
                pure = false;
                break;
            }
        }
        if (pure) {
            add_terminal_node(node_key, data, sample_keys);
            return false;
        }
    }

    { /* Draw a random subset of variables to possibly split at - then find best
       * split (implementation-specific) */
        key_vector split_candidate_keys = draw_candidates();
        const bool split_found = push_best_split(
            node_key, data, sample_keys, split_candidate_keys
        );
        /* node_key, split_candidate_keys, n_node_sample?? */
        if (!split_found) {
            add_terminal_node(node_key, data, sample_keys);
            return false;
        }
    }

  /* If we reach here, then find_best_split will have appended to split_keys
   * and split_values */

    const size_t split_key = split_keys[node_key];
    const double split_value = split_values[node_key];
  /* NOTE: No longer unpermute split_key */

  /* Initialise left and right children */
    const size_t left_key = split_keys.size();
    left_children[node_key] = left_key;
    push_back_empty_node();
    start_pos[left_key] = start_pos[node_key];

    const size_t right_key = split_keys.size();
    right_children[node_key] = right_key;
    push_back_empty_node();
    start_pos[right_key] = end_pos[node_key];

  /* Partially sorts the keys within the current node; any (sample) key whose
   * data is 'left' of the split will be placed earlier in the sequence than any
   * key 'right' of the split */
    if ((*is_ordered)[split_key]) {
      /* ordered factor: 'left' of the split if value <= `split_value`. */
        for (size_t j = start_pos[node_key]; j < start_pos[right_key]; ) {
            const size_t key = sample_keys[j];
            if (data->get_x(key, split_key) <= split_value) {
              /* */
                ++j;
            } else {
                // If going to right, move to right end
                --start_pos[right_key];
                std::swap(sample_keys[j], sample_keys[start_pos[right_key]]);
            }
        }
    } else {
      /* Unordered factor (partitioning): `partition_key` is a bit-wise encoding
       * of which factor levels are left (or not) - i.e. if an observed factor
       * has numeric value of 5, and bit 5 of the key is zero - the observation
       * is to the left of the split.
       *
       * NOTE: Casting of double to ull - unsafe? */
        size_t j = start_pos[node_key];
        const ull_bitenc partition_key = *((unsigned long long *)(&split_value));
        while (j < start_pos[right_key]) {
            const size_t key = sample_keys[j];
            const size_t obs_bit = std::floor(data->get_x(key, split_key) - 1);
            if (!partition_key.test(obs_bit)) {
                ++j;
            } else {
                --start_pos[right_key];
                std::swap(sample_keys[j], sample_keys[start_pos[right_key]]);
            }
        }
    }

  /* End position of the left child is the start position of the right */
    end_pos[left_key] = start_pos[right_key];
    end_pos[right_key] = end_pos[node_key];
  /* Indicate that we succeeded in splitting */
    return true;

}


inline key_vector TreeBase::draw_candidates() {

    key_vector result;
    count_vector inbag_counts = count_vector(n_predictor, 0);

    if (draw_predictor_weights->empty()) {
        draw_no_replace(n_try, n_predictor, *draw_always_predictor_keys,
                        gen, result,
                        inbag_counts);
    } else {
        draw_no_replace_weighted(n_try, *draw_predictor_weights, gen,
                                 result, inbag_counts);
    }

    result.reserve(result.size() + draw_always_predictor_keys->size());
    std::copy(draw_always_predictor_keys->cbegin(),
              draw_always_predictor_keys->cend(),
              std::back_inserter(result));

    return result;

}


inline void TreeBase::finalise_growth() { /* Default does nothing */ }


inline void TreeBase::push_back_empty_node_impl() { /* Default does nothing */ }


} /* namespace literanger */


#endif /* LITERANGER_TREE_BASE_DEFN_H */

