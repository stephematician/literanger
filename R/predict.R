# -------------------------------------------------------------------------------
# This has been adapted from Ranger.
#
# Ranger is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Ranger is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Ranger. If not, see <http://www.gnu.org/licenses/>.
#
# Written by:
#
#   Marvin N. Wright
# Institut fuer Medizinische Biometrie und Statistik
# Universitaet zu Luebeck
# Ratzeburger Allee 160
# 23562 Luebeck
# Germany
#
# http://www.imbs-luebeck.de
# -------------------------------------------------------------------------------
#
# The adaptation was performed by Stephen Wade. The changes include eliminating
# features that were not required for the multiple imputation procedure. These
# include (for now):
#
# 1.  `type` argument; the predictions are always drawn from the original
#     (in-bag) data for the node.
# 2.  Probability and survival forests are not supported.
# 3.  GWAS data (gwaa.data) not supported.
# 4.  No standard error estimate provided.

#' Ranger prediction
#'
#' Predict with new data and a ranger object from Ranger.
#'
#' The predicted class (classification) or numeric value (regression) is
#' returned by drawing a tree and then a sample from the tree for each row of
#' the data.
#'
#' Factor levels are returned as numerics. To retrieve the corresponding factor
#' levels, use `rf$forest$levels`, if `rf` is the ranger object.
#'
#' @param object A trained random forest `literanger` object.
#' @param newdata New test data of class `data.frame` or `dgCmatrix``.
#' @param seed Random seed. Default is `NULL`, which generates the seed from
#'   `R`. Set to `0` to ignore the `R` seed.
#' @param num.threads Number of threads. Default is number of CPUs available.
#' @param verbose Verbose output on or off.
#'
#' @return Object of class `literanger.prediction` with elements:
#' \describe{
#'   \item{`predictions`}{Predicted (drawn) classes/value for classification and
#'     regression.}
#'   \item{`num.trees`}{Number of trees.}
#'   \item{`num.independent.variables`}{Number of independent variables.}
#'   \item{`treetype`}{Type of forest/tree. Classification, regression or
#'     survival.}
#'   \item{`num.samples`}{Number of samples.}
#' }
#'
#' @examples
#' ## Classification forest
#' train.idx <- sample(nrow(iris), 2/3 * nrow(iris))
#' iris.train <- iris[train.idx, ]
#' iris.test <- iris[-train.idx, ]
#' rg.iris <- train_lite(data=iris.train, response_name="Species")
#' pred.iris <- predict(rg.iris, newdata=iris.test)
#' table(iris.test$Species, pred.iris$values)
#'
#' @references
#' -   Wright, M. N. & Ziegler, A. (2017). ranger: A Fast Implementation of
#'     Random Forests for High Dimensional Data in C++ and R. J Stat Softw
#'     77:1-17. \doi{10.18637/jss.v077.i01}.
#'
#' @seealso \code{\link{train}}
#' @author Marvin N. Wright
#' @export
#' @md
predict.literanger <- function(
    object, newdata=NULL, prediction_type=c("bagged", "doove"),
    seed=sample.int(n=.Machine$integer.max, size=1),
    n_thread=0, verbose=F, ...
) {

  # if (is.null(forest$num.trees) ||
  #       is.null(forest$child.nodeIDs) || is.null(forest$split.varIDs) ||
  #       is.null(forest$split.values) || is.null(forest$independent.variable.names) ||
  #       is.null(forest$treetype)) {
  #   stop("Error: Invalid forest object.")
  # }

  # Non-quantile prediction ?
    if (is.null(newdata))
        stop("Argument 'newdata' is required for prediction.")
    x <- newdata

  # Sparse matrix data
    if (inherits(x, "Matrix") && !inherits(x, "dgCMatrix"))
        stop("Currently only sparse data of class 'dgCMatrix' ",
             "supported.")

    input_is_matrix <- is.matrix(x) || inherits(x, "Matrix")
    var_names_x <- if (input_is_matrix) colnames(x) else names(x)

    predictor_names <- object$predictor_names
    if (input_is_matrix && !all(predictor_names %in% var_names_x))
        stop("One or more independent variables not found in data.")

    if (!identical(var_names_x, predictor_names)) {
        if (input_is_matrix) {
            x <- x[, var_names_x %in% predictor_names, drop=F]
        } else {
            x <- x[predictor_names]
        }
    }

  # Recode characters
    if (!input_is_matrix) {
        char_columns <- sapply(x, is.character)
        x[char_columns] <- lapply(x[char_columns], factor)
    }

  # Recode factors if forest grown in 'ordered' mode
    predictor_levels <- object$predictor_levels
    factor_ind <- !sapply(predictor_levels, is.null)
    if (length(factor_ind) > 0) {

        if (input_is_matrix) stop("Expected newdata to be a data.frame",
            "factors")

        x[factor_ind] <- mapply(
            function(x_j, levels_j) {
                if (!identical(setdiff(levels(x_j), levels_j), character()))
                    warning("Predictor levels found that were",
                        "not present during training")
                factor(x_j, levels=union(levels_j, levels(x_j)), exclude=NULL)
            },
            x[factor_ind], object$predictor_levels[factor_ind], SIMPLIFY=F
        )

    }

  # Check missing values
    if (any(is.na(x))) {
        offending_columns <- predictor_names[colSums(is.na(x)) > 0]
        stop("Missing values in the predictors: ",
             paste0(offending_columns, collapse = ", "), ".", call.=F)
    }
 
    prediction_type <- match.arg(prediction_type)

  # Num threads; default 0 => detect from system in C++.
    if (!is.numeric(n_thread) || n_thread < 0)
        stop("Invalid value for n_thread")

  # Use sparse matrix
    if (inherits(x, "dgCMatrix")) {
        sparse_x <- x
        x <- matrix(numeric())
    } else {
        sparse_x <- NULL
        if (is.data.frame(x)) x <- data.matrix(x)
        if (storage.mode(x) %in% "integer") storage.mode(x) <- "double"
    }

    result <- cpp11_predict(
        object, x, sparse_x,
        prediction_type, seed, n_thread,
        verbose
    )

    if (length(result) == 0) stop("User interrupt or internal error.")

    result$tree_type <- object$tree_type

    if (is.numeric(result$values) && !is.null(object$response_levels))
        result$values <- factor(result$values,
                                levels=seq_along(object$response_levels),
                                labels=object$response_levels)
    
    class(result) <- "literanger.prediction"

    invisible(result)

}

