/*-------------------------------------------------------------------------------
 * This file is part of 'literanger'. literanger was adapted from the 'ranger'
 * package for R Statistical Software <https://www.r-project.org>. ranger was
 * authored by Marvin N Wright with the GNU General Public License version 3.
 * The adaptation was performed by Stephen Wade in 2023. literanger carries the
 * same license, terms, and
 * permissions as ranger.
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
#include "cpp11_predict.decl.h"

/* standard library headers */
#include <memory>
#include <stdexcept>
#include <string>
#include <thread>

/* cpp11 and R headers */
#include "cpp11.hpp"

/* eigen (3) sparse data headers */
#include "eigen3/Eigen/Sparse"

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
cpp11::list cpp11_predict(
    cpp11::list object,
    cpp11::doubles_matrix<> x, cpp11::sexp sparse_x,
    std::string prediction_type, const size_t seed,
    const size_t n_thread, const bool verbose
) {

    using namespace literanger;
    using namespace cpp11::literals;
    cpp11::writable::list result;

    const std::string tree_type = cpp11::as_cpp<std::string>(object["tree_type"]);
    cpp11::external_pointer<ForestBase> forest = object["cpp11_ptr"];

    std::shared_ptr<Data> data { };

    toggle_print print_out { verbose, Rprintf };
    R_user_interruptor user_interrupt { };

    Eigen::SparseMatrix<double> eigen_x;
    const bool use_sparse = sparse_x != R_NilValue;

    if (use_sparse) {
      /* Convert a double (valued) compressed column-major matrix (dgCmatrix) to
       * Eigen::sparseMatrix<double>. */
        cpp11::integers sp_i    = { sparse_x.attr("i") };
        cpp11::integers sp_p    = { sparse_x.attr("p") };
        cpp11::doubles sp_x    = { sparse_x.attr("x") };
        cpp11::integers sp_Dim = { sparse_x.attr("Dim") };

        if (!sp_Dim[1])
            throw std::invalid_argument("Invalid dimension for sparse matrix.");

        eigen_x.resize(sp_Dim[0], sp_Dim[1]);
        eigen_x.reserve(sp_i.size());
        for (size_t j_out = 0; j_out != (size_t)sp_Dim[1]; ++j_out) {
            for (size_t j = sp_p[j_out]; j != (size_t)sp_p[j_out+1]; ++j)
                eigen_x.insert(sp_i[j], j_out) = sp_x[j];
        }
    }

    const cpp11::writable::doubles_matrix<> y(
        use_sparse ? eigen_x.rows() : x.nrow(), 1
    );

    if (use_sparse) {
        data = std::shared_ptr<Data>(new DataSparse(eigen_x, y));
    } else {
        data = std::shared_ptr<Data>(new DataR(x, y));
    }

  /* Perform predictions */
    const size_t predict_n_thread = n_thread == DEFAULT_N_THREAD ?
        std::thread::hardware_concurrency() : n_thread;
    if (predict_n_thread == 0)
        throw std::domain_error("'n_thread' must be positive.");

    switch(as_tree_type(tree_type)) {
    /* NOTE: return type will eventually depend upon tree and prediction type */
    case TREE_CLASSIFICATION: {
        auto & forest_impl = dynamic_cast<ForestClassification &>(*forest);
        switch(as_prediction_type(prediction_type)) {
        case BAGGED: {
            dbl_vector predictions;
            forest_impl.template predict<BAGGED>(data, seed, predict_n_thread,
                                                 user_interrupt,
                                                 predictions, print_out);
            result.push_back({"values"_nm=predictions});
        } break;
        case INBAG: {
            dbl_vector predictions;
            forest_impl.template predict<INBAG>(data, seed, predict_n_thread,
                                                user_interrupt,
                                                predictions, print_out);
            result.push_back({"values"_nm=predictions});
        } break;
       // case NODES: {
       //     std::vector<key_vector> predictions;
       //     forest_impl.template predict<NODES>(data, seed, predict_n_thread,
       //                                         user_interrupt,
       //                                         predictions, print_out);
       //     result.push_back({"values"_nm=predictions});
       // } break; }
        default: throw std::invalid_argument("Unsupported prediction type."); }

    } break;
    case TREE_REGRESSION: {
        auto & forest_impl = dynamic_cast<ForestRegression &>(*forest);
        switch(as_prediction_type(prediction_type)) {
        case BAGGED: {
            dbl_vector predictions;
            forest_impl.template predict<BAGGED>(data, seed, predict_n_thread,
                                                 user_interrupt,
                                                 predictions, print_out);
            result.push_back({"values"_nm=predictions});
        } break;
        case INBAG: {
            dbl_vector predictions;
            forest_impl.template predict<INBAG>(data, seed, predict_n_thread,
                                                user_interrupt,
                                                predictions, print_out);
            result.push_back({"values"_nm=predictions});
        } break;
       // case NODES: {
       //     std::vector<key_vector> predictions;
       //     forest_impl.template predict<NODES>(data, seed, predict_n_thread,
       //                                         user_interrupt,
       //                                         predictions, print_out);
       //     result.push_back({"values"_nm=predictions});
       // } break; }
        default: throw std::invalid_argument("Unsupported prediction type."); }

    } break;
    default: throw std::invalid_argument("Unsupported tree type.");
    }

    return result;

}

