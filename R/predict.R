# ------------------------------------------------------------------------------
# This file is part of 'literanger'. literanger was adapted from the 'ranger'
# package for R statistical software. ranger was authored by Marvin N Wright
# with the GNU General Public License version 3. The adaptation was performed by
# Stephen Wade in 2023. literanger carries the same license, terms, and
# permissions as ranger.
#
# literanger is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# literanger is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with literanger. If not, see <http://www.gnu.org/licenses/>.
#
# Written by:
#
#   Stephen Wade
#   Cancer Council New South Wales
#   Woolloomooloo NSW 2011
#   Australia
# ------------------------------------------------------------------------------

#' Literanger prediction
#'
#' 'literanger' provides different types of prediction that may be used in
#' multiple imputation algorithms with random forests. The usual prediction is
#' the 'bagged' prediction, the most frequent value (or the mean) of the in-bag
#' samples in a terminal node. Doove et al (2014) propose a prediction that
#' better matches the predictive distribution as needed for multiple imputation;
#' take a random draw from the observations in the terminal node from a randomly
#' drawn tree in the forest for each predicted value needed. Alternatively, the
#' usual most-frequent-value or mean of the in-bag responses can be used as in
#' missForest (Stekhoven et al, 2014) or miceRanger
#' <https://cran.r-project.org/package=miceRanger> and missRanger
#' <https://cran.r-project.org/package=missRanger>.
#'
#' Forests trained by literanger retain information about the in-bag responses
#' in each terminal node, thus facilitating efficient predictions within a
#' variation on multiple imputation proposed by Doove et al (2014). This type of
#' prediction can be selected by setting `prediction_type="inbag"`, or the usual
#' prediction for classification and regression forests, the most-frequent-value
#' and mean of in bag samples respectively, is given by setting
#' `prediction_type="bagged"`.
#'
#' A list is returned. The `values` item contains the predicted classes or
#' values (classification and regression forests, respectively). Factor levels
#' are returned as factors with the levels as per the original training data.
#'
#' Compared to the original package ranger, literanger excludes certain
#' features:
#'
#' -   Probability, survival, and quantile regression forests.
#' -   Support for class gwaa.data.
#' -   Standard error estimation.
#'
#' @param object A trained random forest `literanger` object.
#' @param newdata Data of class `data.frame`, `matrix`, or `dgCMatrix`
#' (Matrix), for the latter two; must have column names; all predictors named in
#' `object$predictor_names` must be present.
#' @param prediction_type Name of the prediction algorithm; "bagged" is the
#' most-frequent value among in-bag samples for classification, or the mean of
#' in-bag responses for regression; "inbag" predicts by drawing one in-bag
#' response from a random tree for each row; "nodes" (currently unsupported)
#' returns the node keys (ids) of the terminal node from every tree for each
#' row.
#' @param seed Random seed, an integer between 1 and `.Machine$integer.max`.
#' Default generates the seed from `R`, set to `0` to ignore the `R` seed and
#' use a C++ `std::random_device`.
#' @param n_thread Number of threads. Default is determined by system, typically
#' the number of cores.
#' @param verbose Show computation status and estimated runtime.
#' @param ... Ignored.
#'
#' @return Object of class `literanger_prediction` with elements:
#' \describe{
#'   \item{`values`}{Predicted (drawn) classes/value for classification and
#'     regression.}
#'   \item{`tree_type`}{Number of trees.}
#'   \item{`seed`}{The seed supplied to the C++ library.}
#' }
#'
#' @examples
#' ## Classification forest
#' train_idx <- sample(nrow(iris), 2/3 * nrow(iris))
#' iris_train <- iris[ train_idx, ]
#' iris_test  <- iris[-train_idx, ]
#' rf_iris <- train(data=iris_train, response_name="Species")
#' pred_iris_bagged <- predict(rf_iris, newdata=iris_test,
#'                             prediction_type="bagged")
#' pred_iris_inbag  <- predict(rf_iris, newdata=iris_test,
#'                             prediction_type="inbag")
#' # compare bagged vs actual test values
#' table(iris_test$Species, pred_iris_bagged$values)
#' # compare bagged prediction vs in-bag draw
#' table(pred_iris_bagged$values, pred_iris_inbag$values)
#'
#' @author Stephen Wade <stephematician@gmail.com>, Marvin N Wright (original
#' ranger package)
#'
#' @references
#'
#' -   Doove, L. L., Van Buuren, S., & Dusseldorp, E. (2014). Recursive
#'     partitioning for missing data imputation in the presence of interaction
#'     effects. _Computational Statistics & Data Analysis_, 72, 92-104.
#'     \doi{10.1016/j.csda.2013.10.025}.
#' -   Stekhoven, D.J. and Buehlmann, P. (2012). MissForest--non-parametric
#'     missing value imputation for mixed-type data. _Bioinformatics_, 28(1),
#'     112-118. \doi{10.1093/bioinformatics/btr597}.
#' -   Wright, M. N., & Ziegler, A. (2017a). ranger: A Fast Implementation of
#'     Random Forests for High Dimensional Data in C++ and R. _Journal of
#'     Statistical Software_, 77, 1-17. \doi{10.18637/jss.v077.i01}.
#'
#' @seealso \code{\link{train}}
#'
#' @export
#' @md
predict.literanger <- function(
    object, newdata=NULL, prediction_type=c("bagged", "inbag", "nodes"),
    seed=1L + sample.int(n=.Machine$integer.max - 1L, size=1),
    n_thread=0, verbose=FALSE, ...
) {

    if (is.null(newdata))
        stop("Argument 'newdata' is required for prediction.")
    x <- newdata

    if (inherits(x, "Matrix") && !inherits(x, "dgCMatrix"))
        stop("Currently only sparse data of class 'dgCMatrix' ",
             "supported.")

    is_df <- is.data.frame(x)
    x_names <- if(is_df) names(x) else colnames(x)

    predictor_names <- object$predictor_names
    if (!all(predictor_names %in% x_names))
        stop("One or more predictors not found in data.")

    if (!identical(x_names, predictor_names)) {
        if (!is_df) {
            x <- x[, predictor_names, drop=FALSE]
        } else {
            x <- x[predictor_names]
        }
    }

  # Recode characters
    if (is_df) {
        char_columns <- sapply(x, is.character)
        x[char_columns] <- lapply(x[char_columns], factor)
    }

  # Recode factors if forest grown in 'ordered' mode
    predictor_levels <- (object$predictor_levels)[predictor_names]
    factor_ind <- !sapply(predictor_levels, is.null)
    if (length(factor_ind) > 0) {

        if (!is_df) stop("'newdata' must be data.frame with factor class for: ",
            paste(names(predictor_levels)[factor_ind], collapse=", "), ".")

        x[factor_ind] <- mapply(
            function(x_j, levels_j) {
                if (!identical(setdiff(levels(x_j), levels_j), character()))
                    warning("Predictor levels found that were not present ",
                        "during training")
                factor(x_j, levels=union(levels_j, levels(x_j)), exclude=NULL)
            },
            x[factor_ind], predictor_levels[factor_ind], SIMPLIFY=FALSE
        )

    }

  # Check missing values
    if (any(is.na(x))) {
        offending_columns <- predictor_names[colSums(is.na(x)) > 0]
        stop("Missing values in the predictors: ",
             paste0(offending_columns, collapse = ", "), ".", call.=FALSE)
    }

    prediction_type <- match.arg(prediction_type)

    if (round(seed) != seed) warning("Rounding 'seed' to nearest integer.")
    seed <- as.integer(round(seed))

  # Num threads; default 0 => detect from system in C++.
    if (!is.numeric(n_thread) || n_thread < 0)
        stop("Invalid value for 'n_thread'")

  # Use sparse matrix
    if (inherits(x, "dgCMatrix")) {
        sparse_x <- x
        x <- matrix(numeric())
    } else {
        sparse_x <- NULL
        if (is.data.frame(x)) x <- data.matrix(x)
        if (storage.mode(x) %in% "integer") storage.mode(x) <- "double"
    }

  # Call the library's prediction function using cpp11
    result <- cpp11_predict(
        object, x, sparse_x,
        prediction_type, seed, n_thread,
        verbose
    )

    if (length(result) == 0) stop("User interrupt or internal error.")

    result$tree_type <- object$tree_type
    result$seed <- seed

    if (is.numeric(result$values) && !is.null(object$response_levels))
        result$values <- factor(result$values,
                                levels=seq_along(object$response_levels),
                                labels=object$response_levels)

    class(result) <- "literanger_prediction"

    invisible(result)

}

