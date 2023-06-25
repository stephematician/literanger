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
#ifndef LITERANGER_DATA_R_H
#define LITERANGER_DATA_R_H

/* base class definition */
#include "Data.h"

/* standard library headers */
#include <cstddef>
#include <stdexcept>

/* R and cpp11 headers */
#include "cpp11.hpp"


namespace literanger {

/** Data for random forests using matrix predictor and response. */
struct DataR : public Data {

    public:

        /** Construct data from R matrices.
         * @param[in] x Predictor data with one predictor per column and
         * one observation (or case) per row.
         * @param[in] y Response data with one observation (or case) per row and
         * one (or more) response values per column. */
        DataR(const cpp11::doubles_matrix<> x,
              const cpp11::doubles_matrix<> y);

        /** @copydoc Data::~Data */
        virtual ~DataR() override = default;

        /** @copydoc Data::get_x */
        double get_x(const size_t sample_key,
                     const size_t predictor_key,
                     const bool permute = false) const override;

        /** @copydoc Data::get_y */
        double get_y(const size_t sample_key,
                     const size_t column) const override;


    private:

        /** Reference to the matrix of predictor values managed by R. */
        const cpp11::doubles_matrix<> x;
        /** Reference to the matrix of response values managed by R. */
        const cpp11::doubles_matrix<> y;


};


/* Member definitions */

inline DataR::DataR(const cpp11::doubles_matrix<> x,
                    const cpp11::doubles_matrix<> y) :
    Data(y.nrow(), x.ncol()), x(x), y(y) {

    if (y.nrow() != x.nrow())
        throw std::invalid_argument("Mismatch between number of observations "
            "in 'x' and 'y'");

}


inline double DataR::get_x(const size_t sample_key,
                           const size_t predictor_key,
                           const bool permute) const {
    return x(as_row_offset(sample_key, permute), predictor_key);
}


inline double DataR::get_y(const size_t sample_key, const size_t column) const {
    return y(sample_key, column);
}


} /* namespace literanger */


#endif /* LITERANGER_DATA_R_H */

