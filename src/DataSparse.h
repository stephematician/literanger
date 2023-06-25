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
#ifndef LITERANGER_DATA_SPARSE_H
#define LITERANGER_DATA_SPARSE_H

/* base class definition */
#include "Data.h"

/* standard library headers */
#include <algorithm>
#include <cstddef>
#include <stdexcept>

/* R and cpp11 headers */
#include "cpp11.hpp"


namespace literanger {
/* Convert a double (valued) compressed column-major matrix (dgCmatrix) */

/** Data for random forests using sparse-matrix predictor. */
struct DataSparse : public Data {

    public:

        /** Construct data from compressed column sparse matrix and R matrix.
         * @param[in] dim The dimensions of the matrix of predictor values,
         * one column per predictor and one observation (or case) per row
         * @param[in] i Row (indices) of all non-zero values from a compressed
         * column sparse representation of the predictor matrix.
         * @param[in] p Consecutive elements of this vector define the span of
         * each column within \p i.
         * @param[in] v The value of the matrix with row equal to the
         * corresponding element of \p i, and column determined by which
         * consecutive elements of \p v the element number is bound below and
         * above by.
         * @param[in] y Response data with one observation (or case) per row and
         * one (or more) response values per column. */
        template <typename CountT, typename ValueT>
        DataSparse(const CountT dim,
                   const CountT i, const CountT p, const ValueT v,
                   const cpp11::doubles_matrix<> y);

        /** @copydoc Data::~Data */
        virtual ~DataSparse() override = default;

        /** @copydoc Data::get_x */
        double get_x(const size_t sample_key,
                     const size_t predictor_key,
                     const bool permute = false) const override;

        /** @copydoc Data::get_y */
        double get_y(const size_t sample_key,
                     const size_t column) const override;


    private:

        /** Reference to the sparse matrix of predictor values managed by R. */
        const cpp11::integers x_i;
        const cpp11::integers x_p;
        const cpp11::doubles x_v;
        /** Reference to the matrix of response values managed by R. */
        const cpp11::doubles_matrix<> y;


};


/* Member definitions */

template <typename CountT, typename ValueT>
DataSparse::DataSparse(
    const CountT dim, const CountT i, const CountT p, const ValueT v,
    const cpp11::doubles_matrix<> y
) : Data(dim[0], dim[1]), x_i(i), x_p(p), x_v(v), y(y) {

    if (y.nrow() != dim[0])
        throw std::invalid_argument("Mismatch between number of observations "
            "in 'x' and 'y'");

}


inline double DataSparse::get_x(const size_t sample_key,
                                const size_t predictor_key,
                                const bool permute) const {
    /* test this TODO: */
    using int_t = cpp11::integers::value_type;
    const int_t j_start = x_p[predictor_key];
    const int_t j_end = x_p[predictor_key + 1l];
    if (j_start == j_end) return 0.0;

    const int_t row_offset = as_row_offset(sample_key, permute);

    const cpp11::integers::const_iterator i_ptr = std::lower_bound(
        x_i.cbegin() + j_start, x_i.cbegin() + j_end, row_offset
    );
    return i_ptr == (x_i.cbegin() + j_end) || *i_ptr != row_offset ?
        0.0 : x_v[i_ptr - x_i.cbegin()];
}


inline double DataSparse::get_y(const size_t sample_key,
                                const size_t column) const {
    return y(sample_key, column);
}


} /* namespace literanger */


#endif /* LITERANGER_DATA_SPARSE_H */

