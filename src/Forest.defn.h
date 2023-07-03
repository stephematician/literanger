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
#ifndef LITERANGER_FOREST_DEFN_H
#define LITERANGER_FOREST_DEFN_H

/* class declaration */
#include "Forest.decl.h"

/* standard library headers */
#include <algorithm>
#include <future>
#include <mutex>
#include <random>
#include <stdexcept>
#include <utility>
#include <vector>

/* general literanger headers */
#include "utility.h" // equal_split
/* required literanger class definitions */
#include "Data.defn.h"
#include "ForestBase.defn.h"
#include "TreeBase.defn.h"


namespace literanger {

/* construction call definition */

template <typename ImplT>
template <typename... ArgsT>
Forest<ImplT>::Forest(ArgsT &&... args) :
    ForestBase(std::forward<ArgsT>(args)...)
{ }


/* Interface definitions */

template <typename ImplT>
void Forest<ImplT>::plant(const std::shared_ptr<const Data> data,
/* TODO: maybe vector of case-weights? */
                          const dbl_vector_ptr case_weights,
                          const size_t seed,
                          const size_t n_thread,
                          const bool compute_oob_error,
                          const interruptor & user_interrupt,
                          double & oob_error,
                          toggle_print & print_out) {

    ImplT & forest_impl = *(static_cast<ImplT*>(this));

    const size_t n_sample = data->get_n_row();

    /* should check all? */
   // if ((double)n_sample * (*sample_fraction)[0] < 1)
   //     throw std::domain_error("sample_fraction too small (results in zero "
   //         "samples).");

    if (!case_weights->empty() && case_weights->size() != n_sample)
        throw std::invalid_argument("Number of case weights not equal to "
            "number of samples.");

    print_out("Growing trees...\n");
    seed_gen(seed);

    for (size_t j = 0; j != n_tree; ++j) {
        const dbl_vector_ptr sample_fraction =
            tree_parameters[j].sample_fraction;
        for (double p : *sample_fraction) {
            if ((double)n_sample * p < 1)
                throw std::domain_error("'sample_fraction' too small (results "
                    "in zero samples).");
        }
        forest_impl.plant_tree(data, tree_parameters[j]);
    }

    {
        std::uniform_int_distribution<size_t> U_rng { };
        for (size_t j = 0; j != n_tree; ++j) {
            const size_t seed_j = seed == 0 ? U_rng(gen) : (j + 1) * seed;
            trees[j]->seed_gen(seed_j);
        }
    }

    const size_t n_grow_thread = std::min(n_tree, n_thread);
  /* Set the ranges for the threads */
    equal_split(work_intervals, 0, n_tree - 1, n_grow_thread);

    interrupted = false;
    event_count = 0;

    std::vector<std::future<void>> work_result;
    work_result.reserve(n_grow_thread);

    forest_impl.new_growth(data);

    if (compute_oob_error) forest_impl.new_oob_error(data, n_grow_thread);

  /* Start growing trees in threads */
    for (size_t work_index = 0; work_index != n_grow_thread; ++work_index)
        work_result.push_back(std::async(
            std::launch::async,
            &Forest<ImplT>::grow_interval,
            this, work_index, data, case_weights,
            compute_oob_error)
        );

  /* Block until all tree-growth threads have finished */
    show_progress("Growing trees...", n_tree, n_grow_thread,
                  user_interrupt, print_out);
    for (auto & result : work_result) { result.wait(); result.get(); }

    if (interrupted) throw std::runtime_error("User interrupt.");

    if (compute_oob_error) oob_error = forest_impl.finalise_oob_error(data);

    forest_impl.finalise_growth(data);

}


template <typename ImplT>
template <PredictionType prediction_type, typename result_type>
void Forest<ImplT>::predict(const std::shared_ptr<const Data> data,
                            const size_t seed,
                            const size_t n_thread,
                            const interruptor & user_interrupt,
                            result_type & result,
                            toggle_print & print_out) {

    ImplT & forest_impl = *(static_cast<ImplT*>(this));

    print_out("Predicting...\n");
    seed_gen(seed);

    {
        std::uniform_int_distribution<size_t> U_rng { };
        for (size_t j = 0; j != n_tree; ++j) {
            const size_t seed_j = seed == 0 ? U_rng(gen) : (j + 1) * seed;
            trees[j]->seed_gen(seed_j); //
        }
    }

  /* Set the ranges for the threads */
    const size_t n_predict_call = n_tree;
    const size_t n_predict_thread = std::min(n_predict_call, n_thread);
    equal_split(work_intervals, 0, n_predict_call - 1, n_predict_thread);

    interrupted = false;
    event_count = 0;

    std::vector<std::future<void>> work_result;
    work_result.reserve(n_predict_thread);

  /* Initialise a workspace for predictions */
    forest_impl.template new_predictions<prediction_type>(data,
                                                          n_predict_thread);

  /* Generate and store predictions for each tree */
    for (size_t work_index = 0; work_index != n_predict_thread; ++work_index)
        work_result.push_back(std::async(
            std::launch::async,
            &Forest<ImplT>::predict_interval<prediction_type>,
            this, work_index, data)
        );

    show_progress("Predicting...", n_predict_call, n_predict_thread,
                  user_interrupt, print_out);
    for (auto & result : work_result) { result.wait(); result.get(); }

    if (interrupted) throw std::runtime_error("User interrupt.");

  /* Aggregate over trees for each observation/sample */
    for (size_t item_key = 0; item_key != data->get_n_row(); ++item_key)
        forest_impl.template aggregate_one_item<prediction_type>(item_key);

  /* Clean up the workspace and return the prediction. */
    forest_impl.template finalise_predictions<prediction_type>(result);

}


template <typename ImplT>
void Forest<ImplT>::grow_interval(
    const size_t work_index,
    const std::shared_ptr<const Data> data,
    const dbl_vector_ptr case_weights,
    const bool compute_oob_error
) {
    BEGIN_WORKER

    ImplT & forest_impl = *(static_cast<ImplT*>(this));

    if (work_index >= work_intervals.size() - 1) return;

    const size_t start = work_intervals[work_index];
    const size_t end = work_intervals[work_index + 1];

    for (size_t tree_key = start; tree_key != end; ++tree_key) {
        key_vector oob_keys = trees[tree_key]->grow(data, case_weights,
                                                    compute_oob_error);

        if (compute_oob_error)
            forest_impl.oob_one_tree(tree_key, data, oob_keys);

        std::unique_lock<std::mutex> lock(mutex);
        if (interrupted) { condition_variable.notify_one(); return; }
        ++event_count;
        condition_variable.notify_one();

    }

    END_WORKER
}


template <typename ImplT>
template <PredictionType prediction_type>
void Forest<ImplT>::predict_interval(
    const size_t work_index,
    const std::shared_ptr<const Data> data
) {
    BEGIN_WORKER

    ImplT & forest_impl = *(static_cast<ImplT*>(this));

    if (work_index >= work_intervals.size() - 1) return;

    const size_t start = work_intervals[work_index];
    const size_t end = work_intervals[work_index + 1];

    key_vector sample_keys(data->get_n_row(), 0);
    std::iota(sample_keys.begin(), sample_keys.end(), 0);

    for (size_t tree_key = start; tree_key != end; ++tree_key) {

      /* this needs to store the prediction in the workspace */
        forest_impl.template predict_one_tree<prediction_type>(
            tree_key, data, sample_keys
        );

        std::unique_lock<std::mutex> lock(mutex);
        if (interrupted) { condition_variable.notify_one(); return; }
        ++event_count;
        condition_variable.notify_one();

    }

    END_WORKER
}


} /* namespace literanger */


#endif /* LITERANGER_FOREST_DEFN_H */

