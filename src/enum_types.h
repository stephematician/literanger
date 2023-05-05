/* This file was adapted from the C++ core of the "ranger" package for R
 * Statistical Software.
 *
 * Adaptation was authored by Stephen Wade. The same license terms as the
 * original c++ core of the ranger package apply to the adaptation.
 *
 * License statement for c++ core of ranger:
 *
 * Copyright (c) [2014-2018] [Marvin N. Wright]
 *
 * This software may be modified and distributed under the terms of the MIT
 * license.
 *
 * Please note that the C++ core of ranger is distributed under MIT license and
 * the R package "ranger" under GPL3 license.
 */
#ifndef LITERANGER_ENUM_TYPES_H
#define LITERANGER_ENUM_TYPES_H

/* standard library headers */
#include <cstddef>
#include <iterator>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <unordered_map>


namespace literanger {

/** Enumerated tree types. */
enum TreeType { TREE_CLASSIFICATION, TREE_REGRESSION };


/** Enumerated rules for selecting a predictor to split on. */
enum SplitRule { LOGRANK, MAXSTAT, EXTRATREES, BETA, HELLINGER };


/** Enumerated types of prediction. */
enum PredictionType {
    BAGGED, /**< Each predicted value is bootstrap-aggregated over all trees */
    DOOVE   /**< Each predicted value comes from one randomly-sampled tree */
};


template <PredictionType prediction_type>
using enable_if_bagged =
    typename std::enable_if<prediction_type == BAGGED, std::nullptr_t>::type;

template <PredictionType prediction_type>
using enable_if_doove =
    typename std::enable_if<prediction_type == DOOVE, std::nullptr_t>::type;


/* Declarations */

/** Convert a string to enumerated tree type.
 * @param x "classification" or "regression" only supported. */
TreeType as_tree_type(std::string x);
/** Convert a string to enumerated splitting rule.
 * @param x e.g. "gini", "variance", etc. */
SplitRule as_split_rule(std::string x);
/** Convert a string to enumerated prediction type.
 * @param x "bagged" or "doove" only. */
PredictionType as_prediction_type(std::string x);


/* Definitions */

inline TreeType as_tree_type(std::string x) {

    static std::unordered_map<std::string,TreeType> table = {
        { "classification", TreeType::TREE_CLASSIFICATION },
        { "regression", TreeType::TREE_REGRESSION }
    };

    const auto it = table.find(x);

    if (it == std::end(table))
        throw std::invalid_argument("Invalid tree type.");

    return it->second;

}


inline SplitRule as_split_rule(std::string x) {

    static std::unordered_map<std::string,SplitRule> table = {
        { "gini", SplitRule::LOGRANK },
        { "variance", SplitRule::LOGRANK },
        { "maxstat", SplitRule::MAXSTAT },
        { "extratrees", SplitRule::EXTRATREES },
        { "beta", SplitRule::BETA },
        { "hellinger", SplitRule::HELLINGER }
    };

    auto it = table.find(x);

    if (it == std::end(table))
        throw std::invalid_argument("Invalid split metric.");

    return it->second;

}


inline PredictionType as_prediction_type(std::string x) {

    static std::unordered_map<std::string,PredictionType> table = {
        {"bagged", PredictionType::BAGGED},
        {"doove", PredictionType::DOOVE}
    };

    auto it = table.find(x);

    if (it == std::end(table))
        throw std::invalid_argument("Invalid prediction type.");

    return it->second;

}


} /* namespace literanger */


#endif /* LITERANGER_ENUM_TYPES_H */

