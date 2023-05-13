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
#ifndef LITERANGER_CPP11_UTILITY_H
#define LITERANGER_CPP11_UTILITY_H

/* standard library headers */
#include <algorithm>
#include <memory>
#include <vector>

/* R and cpp11 headers */
#include "cpp11.hpp"
#include "cpp11/R.hpp"
#include "R.h"

/* required literanger class definition */
#include "utility_interrupt.h"


namespace literanger {

/** Test for user interruption. @param x Ignored. */
static void chk_int_fn(void * x) { R_CheckUserInterrupt(); }


/* Declarations */

/** User interrupt from main R thread */
struct R_user_interruptor : public interruptor {
    /** Call user-interrupt check at top level of execution. */
    bool operator()() const;
};


/** Convert 1-dimensional cpp11 object into a vector.
 *
 * @tparam FromT the type of the cpp11 to convert
 * @tparam T the value type of the returned std::vector */
template <typename T, typename FromT>
std::vector<T> as_vector(const FromT x);


/** Convert 1-dimensional cpp11 object into a ptr<vector>.
 *
 * @tparam FromT the type of the cpp11 to convert
 * @tparam T the value type of the returned std::vector */
template <typename T, typename FromT,
          template <typename...> class PtrT = std::shared_ptr>
PtrT<std::vector<T>> as_vector_ptr(const FromT x);


/** Convert a 'nested' container of cpp11 objects into vector<vector>.
 *
 */
template <typename T, typename NestedT, typename FromT>
std::vector<std::vector<T>> as_nested(const FromT x);


/** Convert a 'nested' container of cpp11 objects into vector<ptr<vector>>.
 *
 */
template <typename T, typename NestedT, typename FromT,
          template <typename...> class PtrT = std::shared_ptr>
std::vector<PtrT<std::vector<T>>> as_nested_ptr(const FromT x);


/* Definitions */

inline bool R_user_interruptor::operator()() const {
    return (R_ToplevelExec(chk_int_fn, NULL) == FALSE);
}


template <typename T, typename FromT>
std::vector<T> as_vector(const FromT x) {
    return cpp11::as_cpp<std::vector<T>>(x);
}


template <typename T, typename FromT,
          template <typename...> class PtrT>
PtrT<std::vector<T>> as_vector_ptr(const FromT x) {
    return PtrT<std::vector<T>>( new std::vector<T>(as_vector<T>(x)) );
}


template <typename T, typename NestedT, typename FromT>
std::vector<std::vector<T>> as_nested(const FromT x) {
    std::vector<std::vector<T>> value(x.size());
    std::transform(
        x.cbegin(), x.cend(), value.begin(),
        [](const typename FromT::const_iterator::value_type & item){
            return as_vector<T>(NestedT(item));
        }
    );
    return value;
}


template <typename T, typename NestedT, typename FromT,
          template <typename...> class PtrT>
std::vector<PtrT<std::vector<T>>> as_nested_ptr(const FromT x) {
    std::vector<PtrT<std::vector<T>>> value(x.size());
    std::transform(
        x.cbegin(), x.cend(), value.begin(),
        [](const typename FromT::const_iterator::value_type & item){
            return as_vector_ptr<T>(NestedT(item));
        }
    );
    return value;
}


} /* namespace literanger */


#endif /* LITERANGER_CPP11_UTILITY_H */

