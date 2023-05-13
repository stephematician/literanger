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


#' Train forest using ranger for multiple imputation algorithms.
#'
#' 'literanger' trains random forests for use in multiple imputation problems
#' via an adaptation of the 'ranger' R package. ranger is a fast
#' implementation of random forests (Breiman, 2001) or recursive partitioning,
#' particularly suited for high dimensional data (Wright et al, 2017a).
#' literanger supports prediction used in algorithms such as "Multiple
#' Imputation via Chained Equations" (Van Buuren, 2007).
#'
#' literanger trains classification and regression forests using the original
#' Random Forest (Breiman, 2001) or extremely randomized trees (Geurts et
#' al, 2006) algorithms. The trained forest retains information about the in-bag
#' responses in each terminal node, thus facilitating a variation on the
#' algorithm for multiple imputation with random forests proposed by Doove et
#' al (2014). This algorithm should match the predictive distribution more
#' closely than using predictive mean matching.
#'
#' The default split metric for classification trees is the Gini impurity, which
#' can be extended to use the extra-randomized trees rule (Geurts et al, 2006).
#' For binary responses, the Hellinger distance metric may be used instead
#' (Cieslak et al, 2012).
#'
#' The default split metric for regression trees is the estimated variance,
#' which can be extended to include the extra-randomized trees rule, too.
#' Alternatively, the beta log-likelihood (Wright et al, 2017b) or maximally
#' selected rank statistics (Wright et al, 2019) are available.
#'
#' When the `data` and `response_name` arguments are supplied the response
#' variable is identified by its corresponding column name. The type of response
#' may be used to determine the type of tree. If the response is a factor then
#' classification trees are used. If the response is numeric then regression
#' trees are used. The `classification` argument can be used to override the
#' default tree type when the response is numeric. Alternatively, use `x` and
#' `y` arguments to specify response and predictor; this can avoid conversions
#' and save memory. If memory usage issues persist, consider setting
#' `save_memory=TRUE` but be aware that this option slows down the tree growing.
#'
#' The `min_split_n_sample` rule can be used to control the minimum number of
#' in-bag samples required to split a node; thus, as in the original algorithm,
#' nodes with fewer samples than `min_split_n_sample` are possible. To put a
#' floor under the number of samples per node, the `min_leaf_n_sample`
#' argument is used.
#'
#' When drawing candidate predictors for splitting a node on, the predictors
#' identified by `names_of_always_draw` are included _in addition_ to the
#' `n_try` predictors that are randomly drawn. Another way to modify the way
#' predictors are selected is via the `draw_predictor_weights` argument, which
#' are normalised and interpreted as probabilities when drawing candidates. The
#' weights are assigned _in the order they appear in the data_. Weights assigned
#' by `draw_predictor_weights` to variables in `names_of_always_draw` are
#' ignored. The usage of `draw_predictor_weights` can increase the computation
#' times for large forests.
#'
#' Unordered-factor predictors can be handled in 3 different ways by using
#' `unordered_predictors`:
#'
#' -   For "ignore" all factors are regarded ordered;
#' -   For "partition" all possible 2-partitions are considered for splitting.
#' -   For "order" and 2-class classification the factor levels are ordered by
#'     their proportion falling in the second class, for regression by their
#'     mean response, as described in Hastie et al. (2009), chapter 9.2.4. For
#'     multi-class classification the factor levels are ordered by the first
#'     principal component of the weighted covariance matrix of the contingency
#'     table (Coppersmith et al, 1999).
#'
#' The use of "order" is recommended, as it computationally fast and
#' can handle an unlimited number of factor levels. Note that the factors are
#' only reordered once and not again in each split.
#'
#' Compared to the original package ranger, literanger excludes certain
#' features:
#'
#' -   Formula interface.
#' -   Probability, survival, and quantile regression forests.
#' -   Support for class gwaa.data.
#' -   Measures of variable importance.
#' -   Regularisation of importance.
#' -   Access to in-bag data via R.
#' -   Support for user-specified hold-out data.
#'
#' @param data Training data of class `data.frame`, `matrix`, or `dgCMatrix`
#' (Matrix), for the latter two; must have column names.
#' @param response_name Name of response (dependent) variable if `data` was
#' provided.
#' @param predictor_names Names of predictor (independent) variables if `data`
#' was provided; default is all variables that are not the response.
#' @param x Predictor data (independent variables), alternative interface to
#' `data` and `response_name`.
#' @param y Response vector (dependent variable), alternative interface to
#' `data` and `response_name`.
#' @param case_weights Weights for sampling of training observations.
#' Observations with larger weights will be selected with higher probability
#' in the bootstrap (or sub-sampled) samples for each tree.
#' @param classification Set to `TRUE` to grow a classification forest if the
#' response is numeric (including if data is a matrix), else, a regression
#' forest is grown.
#' @param n_tree Number of trees (default 10).
#' @param replace Sample with replacement to train each tree.
#' @param sample_fraction Fraction of observations to sample to train each tree.
#' Default is 1 for sampling with replacement and 0.632 for sampling without
#' replacement. For classification, this can be a vector of class-specific
#' values.
#' @param n_try Number of variables (predictors) to draw that are candidates for
#' splitting each node by. Default is the (rounded down) square root of the
#' number of predictors. Alternatively, a single argument function returning an
#' integer, given the number of predictors.
#' @param draw_predictor_weights For predictor-drawing weights shared by all
#' trees; a numeric vector of _non-negative_ weights for each predictor. For
#' tree-specific predictor-drawing weights; a list of size `n_tree` containing
#' (non-negative) vectors with length equal to the number of predictors.
#' @param names_of_always_draw Character vector with predictor (variable) names
#' to be selected _in addition_ to the `n_try` predictors drawn as candidates to
#' split by.
#' @param split_rule Splitting rule. For classification estimation "gini",
#' "extratrees" or "hellinger" with default "gini". For regression "variance",
#' "extratrees", "maxstat" or "beta" with default "variance".
#' @param max_depth Maximal tree depth. A value of NULL or 0 (the default)
#' corresponds to unlimited depth, 1 to tree stumps (1 split per tree).
#' @param min_split_n_sample Minimal number of in-bag samples a node must have
#' in order to be split. Default 1 for classification and 5 for regression.
#' @param min_leaf_n_sample Minimum number of in-bag samples in a leaf node.
#' @param unordered_predictors Handling of unordered factor predictors. One of
#' "ignore", "order" and "partition". For the "extratrees" splitting rule the
#' default is "partition" for all other splitting rules "ignore".
#' @param response_weights Classification only: Weights for the response classes
#' (in order of the factor levels) in the splitting rule e.g. cost-sensitive
#' learning. Weights are also used by each tree to determine majority vote.
#' @param n_random_split "extratrees" split metric only: Number of random splits
#' to consider for each candidate splitting variable, default is 1.
#' @param alpha "maxstat" splitting rule only: Significance threshold to allow
#' splitting, default is 0.5, must be in the interval `(0,1]`.
#' @param min_prop "maxstat" splitting rule only: Lower quantile of covariate
#' distribution to be considered for splitting, default is 0.1, must be in the
#' interval `[0,0.5]`.
#' @param seed Random seed, an integer between 1 and `.Machine$integer.max`.
#' Default generates the seed from `R`, set to `0` to ignore the `R` seed and
#' use a C++ `std::random_device`.
#' @param save_memory Use memory saving (but slower) splitting mode. Warning:
#' This option slows down the tree growing, use only if you encounter memory
#' problems.
#' @param n_thread Number of threads. Default is determined by system, typically
#' the number of cores.
#' @param verbose Show computation status and estimated runtime.
#'
#' @return Object of class `literanger` with elements:
#' \describe{
#'   \item{`predictor_names`}{The names of the predictor variables in the
#'     model.}
#'   \item{`names_of_unordered`}{The names of predictors that are unordered.}
#'   \item{`tree_type`}{The type of tree in the forest.}
#'   \item{`n_tree`}{The number of trees that were trained.}
#'   \item{`n_try`}{The number of predictors drawn as candidates for each
#'     split.}
#'   \item{`split_rule`}{The name of the split metric used.}
#'   \item{`max_depth`}{The maximum allowed depth of a tree in the forest.}
#'   \item{`min_metric_decrease`}{The minimum decrease in the metric for an
#'     acceptable split (equal to negative @p alpha for maximally selected
#'     rank statistics, else zero).}
#'   \item{`min_split_n_sample`}{The minimum number of in-bag samples in a node
#'     prior to splitting.}
#'   \item{`min_leaf_n_sample`}{The minimum number of in-bag samples in a leaf
#'     node.}
#'   \item{`seed`}{The seed supplied to the C++ library.}
#'   \item{`oob_error`}{The misclassification rate or the mean square error
#'     using out-of-bag samples.}
#'   \item{`cpp11_ptr`}{An external pointer to the trained forest. DO NOT
#'     MODIFY.}
#'   \item{`response_values`}{Classification only: the values of the response in
#'     the order they appear in the data.}
#'   \item{`response_levels`}{Classification only: the labels for the response
#'     in the order they appear in the data.}
#' }
#'
#' @examples
#' ## Classification forest with default settings
#' train(data=iris, response_name="Species")
#'
#' ## Prediction
#' train_idx <- sample(nrow(iris), 2/3 * nrow(iris))
#' iris_train <- iris[train_idx, ]
#' iris_test <- iris[-train_idx, ]
#' rg_iris <- train(data=iris_train, response_name="Species")
#' pred_iris <- predict(rg_iris, newdata=iris_test)
#' table(iris_test$Species, pred_iris$values)
#'
#' @author Stephen Wade <stephematician@gmail.com>, Marvin N Wright (original
#' ranger package)
#'
#' @references
#'
#' -   Breiman, L. (2001). Random forests. _Machine Learning_, 45, 5-32.
#'     \doi{10.1023/A:1010933404324}.
#' -   Cieslak, D. A., Hoens, T. R., Chawla, N. V., & Kegelmeyer, W. P. (2012).
#'     Hellinger distance decision trees are robust and skew-insensitive. _Data
#'     Mining and Knowledge Discovery_, 24, 136-158.
#'     \doi{10.1007/s10618-011-0222-1}.
#' -   Coppersmith, D., Hong, S. J., & Hosking, J. R. (1999). Partitioning
#'     nominal attributes in decision trees. _Data Mining and Knowledge
#'     Discovery_, 3, 197-217. \doi{10.1023/A:1009869804967}.
#' -   Doove, L. L., Van Buuren, S., & Dusseldorp, E. (2014). Recursive
#'     partitioning for missing data imputation in the presence of interaction
#'     effects. _Computational Statistics & Data Analysis_, 72, 92-104.
#'     \doi{10.1016/j.csda.2013.10.025}.
#' -   Geurts, P., Ernst, D., & Wehenkel, L. (2006). Extremely randomized trees.
#'     _Machine Learning_, 63, 3-42. \doi{10.1007/s10994-006-6226-1}.
#' -   Hastie, T., Tibshirani, R., Friedman, J. H., & Friedman, J. H. (2009).
#'     The elements of statistical learning: data mining, inference, and
#'     prediction (Vol. 2). New York: Springer.
#'     \doi{10.1007/978-0-387-21606-5}.
#' -   Van Buuren, S. (2007). Multiple imputation of discrete and continuous
#'     data by fully conditional specification. _Statistical Methods in Medical
#'     Research_, 16(3), 219-242. \doi{10.1177/0962280206074463}.
#' -   Weinhold, L., Schmid, M., Wright, M. N., & Berger, M. (2019). A random
#'     forest approach for modeling bounded outcomes. _arXiv preprint_,
#'     arXiv:1901.06211. \doi{10.48550/arXiv.1901.06211}.
#' -   Wright, M. N., & Ziegler, A. (2017a). ranger: A Fast Implementation of
#'     Random Forests for High Dimensional Data in C++ and R. _Journal of
#'     Statistical Software_, 77, 1-17. \doi{10.18637/jss.v077.i01}.
#' -   Wright, M. N., Dankowski, T., & Ziegler, A. (2017b). Unbiased split
#'     variable selection for random survival forests using maximally selected
#'     rank statistics. _Statistics in medicine_, 36(8), 1272-1284.
#'     \doi{10.1002/sim.7212}.
#'
#' @seealso \code{\link{predict.literanger}}
#'
#' @export
#' @md
train <- function(
    data=NULL, response_name=character(), predictor_names=character(),
    x=NULL, y=NULL,
    case_weights=numeric(),
    classification=NULL, n_tree=10,
    replace=TRUE, sample_fraction=ifelse(replace, 1, 0.632),
    n_try=NULL,
    draw_predictor_weights=numeric(), names_of_always_draw=character(),
    split_rule=NULL, max_depth=0, min_split_n_sample=0, min_leaf_n_sample=0,
    unordered_predictors=NULL,
    response_weights=numeric(),
    n_random_split=1, alpha=0.5, min_prop=0.1,
    seed=1L + sample.int(n=.Machine$integer.max - 1L, size=1),
    save_memory=FALSE, n_thread=0, verbose=FALSE
) {

    if (is.null(data) && !is.null(x) && !is.null(y)) {
        predictor_names <- colnames(x)
    } else if (!is.null(data) && is.null(x) && is.null(y)) {
        if (!is.character(response_name))
            stop("'response_name' must be a character vector.")
        if (!is.character(predictor_names))
            stop("'predictor_names' must be a character vector (empty for ",
                "default).")
        if (!any(response_name %in% colnames(data)) &&
                !all(response_name %in% colnames(data)))
            stop("'response_name' variable not found in data.")
        if (!all(predictor_names %in% colnames(data)))
            stop("'predictor_names' contained variables not found in data.")
        if (identical(predictor_names, character()))
            predictor_names <- setdiff(colnames(data), response_name)
        y <- data[, response_name, drop=TRUE]
        x <- data[, predictor_names, drop=FALSE]
    } else
        stop("Only one of 'data' or 'x'/'y' may be supplied.")

  # Sparse matrix data
    if (inherits(x, "Matrix") && !inherits(x, "dgCMatrix"))
        stop("Currently only sparse data of class 'dgCMatrix' supported.")

  # Check missing values
    if (any(is.na(x))) {
        offending_columns <- predictor_names[colSums(is.na(x)) > 0]
        stop("Missing values in the predictors: ",
             paste0(offending_columns, collapse=", "), ".", call.=FALSE)
    }
    if (any(is.na(y)))
        stop("Missing values in the response.", call.=FALSE)

  # Check response levels
    if (is.factor(y))  {
        if (nlevels(y) != nlevels(droplevels(y))) {
            dropped_levels <- setdiff(levels(y), levels(droplevels(y)))
            warning("Dropped unused factor level(s) in response variable: ",
                    paste0(dropped_levels, collapse=", "), ".", call.=FALSE)
        }
    }

  # Treetype
    if (is.factor(y) || is.logical(y)) {
        tree_type <- "classification"
    } else if (is.numeric(y) && (is.null(ncol(y)) || ncol(y) == 1)) {
        if (!is.null(classification) && classification) {
            tree_type <- "classification"
        } else {
            tree_type <- "regression"
        }
    } else stop("Unsupported type of response variable.")

    n_predictor <- length(predictor_names)
    if (n_predictor < 1)
        stop("No predictors found - missing column names?")

  # Case weights: NULL for no weights or all weights equal
    if (length(unique(case_weights)) == 1) case_weights <- numeric()

  # Number of trees
    if (!is.numeric(n_tree) || n_tree < 1) stop("Invalid value for 'n_tree'.")

  # Sample fraction
    if (!is.numeric(sample_fraction))
        stop("Invalid value for 'sample_fraction'. Please give a value in ",
             "(0,1] or a vector of values in [0,1].")

  # n_try as a function
    if (is.function(n_try)) {

        if (length(formals(n_try)) > 1)
            stop("'n_try' function must have exactly one argument.")

      # Evaluate function
        n_try <- try(n_try(n_predictor), silent=TRUE)

        if (inherits(n_try, "try-error")) {
            message("The 'n_try' function produced the error: ", n_try)
            stop("'n_try' function evaluation resulted in an error.")
        }

      # Check for a single numeric
        if (!is.numeric(n_try) || length(n_try) != 1) {
            stop("'n_try' function must return a single integer or numeric.")
        } else
            n_try <- as.integer(n_try)

      # Check for limits
        if (n_try < 1 || n_try > n_predictor)
            stop("'n_try' function must evaluate to a value not less ",
                 "than 1 and not greater than the number of predictors ( = ",
                 n_predictor, " ).")
    }

  # Check (numeric or otherwise) n_try
    if (is.null(n_try)) {
        n_try <- 0
    } else if (!is.numeric(n_try) || n_try < 0)
        stop("Invalid value for 'n_try'.")

    if (length(sample_fraction) > 1) {

        if (!(tree_type %in% "classification"))
            stop("Invalid value for 'sample_fraction'. Vector values only ",
                 "valid for classification forests.")
        if (any(sample_fraction < 0) || any(sample_fraction > 1))
            stop("Invalid value for 'sample_fraction'. Please give a value in ",
                 "(0,1] or a vector of values in [0,1].")
        if (length(sample_fraction) != nlevels(y))
            stop("Invalid value for 'sample_fraction'. Expecting ", nlevels(y),
                 " values, provided ", length(sample_fraction), ".")
        if (!replace & any(sample_fraction * length(y) > table(y))) {
            j <- which(sample_fraction * length(y) > table(y))[1]
            stop("Not enough samples in class ", names(j),
                 "; available: ", table(y)[j],
                 ", requested: ", (sample_fraction * length(y))[j], ".")
        }
        if (!identical(case_weights, numeric()))
            stop("Combination of 'case_weights' argument and class-wise ",
                 "sampling not supported.")
        # Fix order (C++ needs in order as classes appear in data)
        sample_fraction <- sample_fraction[as.numeric(unique(y))]

    } else if (sample_fraction <= 0 || sample_fraction > 1) {
        stop("Invalid value for 'sample_fraction'. Please give a value in ",
             "(0,1] or a vector of values in [0,1].")
    } else if(!replace && length(case_weights) > 0 &&
                  sum(case_weights > 0) < sample_fraction * nrow(x))
        stop("Fewer non-zero case weights than observations to sample.")

  # Split select weights: NULL for no weights
    if (identical(draw_predictor_weights, numeric())) {
        draw_predictor_weights <- list()
    } else if (is.numeric(draw_predictor_weights)) {
        if (length(draw_predictor_weights) != n_predictor)
        stop("Size of 'draw_predictor_weights' (numeric) not equal to number ",
             "of predictors.")
        draw_predictor_weights <- list(draw_predictor_weights)
    } else if (is.list(draw_predictor_weights)) {
        if (length(draw_predictor_weights) != n_tree)
            stop("Size of 'draw_predictor_weights' (list) not equal to number ",
                 "of trees.")
        if (!all(sapply(draw_predictor_weights, length) == n_predictor))
            stop("One or more vectors in 'draw_predictor_weights' (list) has ",
                 "size not equal to number of predictors.")
    } else stop("Invalid weighting for 'draw_predictor_weights'.")

    missing_always_draw <- setdiff(names_of_always_draw, predictor_names)
    if (length(missing_always_draw) > 0)
        stop("Always drawn variable(s) not found in data: ",
             paste(missing_always_draw, collapse=", "), ".")

  # Split rule
    if (is.null(split_rule))
        split_rule <- c("classification"="gini",
                        "regression"="variance")[tree_type]

    if (tree_type %in% "classification") {
        split_rule <- match.arg(split_rule,
                                c("gini", "extratrees", "hellinger"))
    } else if (tree_type %in% "regression") {
        split_rule <- match.arg(split_rule,
                                c("variance", "extratrees", "maxstat", "beta"))
    }

    if (split_rule == "beta") {
        ## Check for 0..1 outcome
        if (min(y) <= 0 || max(y) >= 1)
            stop("Beta log-likelihood metric applicable to regression data ",
                 "in the interval (0,1).")
    } else if (split_rule == "hellinger" && ((is.factor(y) && nlevels(y) > 2) ||
                   (length(unique(y)) > 2)))
        stop(paste("Hellinger split metric only implemented for binary",
                   "classification."))

  # Importance mode
  #  importance_mode <- match.arg(importance_mode)

  # Minimum node size
    if (!is.numeric(min_split_n_sample) || min_split_n_sample < 0)
        stop("Invalid value for 'min_split_n_sample'.")

  # Tree depth
    if (!is.numeric(max_depth) || max_depth < 0)
        stop("Invalid value for 'max_depth'.")

  # Minimum node size
    if (!is.numeric(min_leaf_n_sample) || min_leaf_n_sample < 0)
        stop("Invalid value for 'min_leaf_n_sample'.")

  # Respect unordered factors
    if (split_rule %in% "extratrees") {
        unordered_predictors <- match.arg(
            unordered_predictors, c("partition", "ignore", "order")
        )
    } else {
        unordered_predictors <- match.arg(
            unordered_predictors, c("ignore", "order", "partition")
        )
    }

    if (split_rule == "extratrees" && save_memory &&
            unordered_predictors == "partition")
        stop("'save_memory' option not possible in extra trees with ",
             "unordered predictors.")

    if (unordered_predictors == "order" && tree_type == "regression" &&
            split_rule == "maxstat")
      warning("The 'order' mode for unordered factor handling with the ",
              "'maxstat' splitrule is experimental.")

 # Recode characters as factors and recode factors if "order" mode
    predictor_levels <- NULL
    if (!is.matrix(x) && !inherits(x, "Matrix") && ncol(x) > 0) {
        character_ind <- sapply(x, is.character)

        if (unordered_predictors %in% "order") {
          # Recode characters and unordered factors
            ordered_ind <- sapply(x, is.ordered)
            factor_ind <- sapply(x, is.factor)
            recode_ind <- character_ind | (factor_ind & !ordered_ind)

          # Numeric response
            if (is.factor(y)) { num_y <- as.numeric(y) } else { num_y <- y }

          # Recode each column
            x[recode_ind] <- lapply(x[recode_ind], function(x_j) {
                if (!is.factor(x_j)) x_j <- as.factor(x_j)
                levels_j <- levels(x_j)
                if (length(levels_j) > 1 && is.factor(y) && nlevels(y) > 2) {
                    levels_j <- pca_order(y=y, x=x_j)
                } else if(length(levels_j) > 1) {
                  # Order factor levels by mean response
                    means <- sapply(levels_j, function(k) {
                        mean(num_y[x_j == k])
                    })
                    levels_j <- as.character(levels_j[order(means)])
                }
              # Return reordered factor
                factor(x_j, levels=levels_j, ordered=TRUE, exclude=NULL)
            })
        } else # Recode characters only
            x[character_ind] <- lapply(x[character_ind], factor)

      # Save levels
        if (any(sapply(x, is.factor))) predictor_levels <- lapply(x, levels)
    }

    if (unordered_predictors == "partition") {
        ordered_idx <- sapply(x, is.ordered)
        factor_idx <- sapply(x, is.factor)
        names_of_unordered <- predictor_names[factor_idx & !ordered_idx]

        if (length(names_of_unordered) > 0) {
            n_level <- sapply(x[, factor_idx & !ordered_idx, drop=FALSE],
                              nlevels)
            limit_n_level <- (.Machine$sizeof.longlong * 8) - 1
            if (max(n_level) > limit_n_level)
                stop("Too many levels in unordered categorical variable ",
                     names_of_unordered[which.max(n_level)],
                     ". Only ", limit_n_level, " levels allowed on this ",
                     "system. Consider using the 'order' option.", sep="")
        } else {
            names_of_unordered <- character()
        }
    } else if (unordered_predictors == "ignore" ||
                   unordered_predictors == "order") {
        names_of_unordered <- character()
    }

    if (length(names_of_unordered) > 0 && split_rule %in% c("beta", "maxstat"))
        stop("Unordered factor splitting not implemented for 'maxstat' or ",
             "'beta' splitting rule.")

  # Class weights: NULL for no weights (all 1)
    if (identical(response_weights, numeric())) {
        if (is.factor(y)) {
            response_weights <- rep(1, nlevels(droplevels(y)))
        } else
            response_weights <- rep(1, length(unique(y)))
    } else {
        if (!(tree_type %in% "classification"))
            stop("'response_weights' argument only valid for classification ",
                 "forests.")
        if (!is.numeric(response_weights) || any(response_weights < 0))
            stop("Invalid value for 'response_weights'. Please give a vector ",
                 "of non-negative values.")
        if (length(response_weights) != length(unique(as.numeric(y))))
            stop("Number of response weights not equal to number of classes.")

        ## Reorder (C++ expects order as appearing in the data)
        response_weights <- response_weights[unique(as.numeric(y))]
    }

  # Extra trees number of random splits
    if (!is.numeric(n_random_split) || n_random_split < 1)
        stop("Invalid value for 'n_random_split', please give a positive ",
             "integer.")

    if (n_random_split > 1 && split_rule != "extratrees")
        warning("Argument 'n_random_split' ignored if 'split_rule' is not ",
                "'extratrees'.")

  # Maxstat splitting
    if (alpha <= 0 || alpha > 1)
        stop("Invalid value for 'alpha', please give a value in (0,1].")

    if (min_prop < 0 || min_prop > 0.5)
        stop("Invalid value for 'min_prop', please give a value in [0,0.5].")

    if (round(seed) != seed) warning("Rounding 'seed' to nearest integer.")
    seed <- as.integer(round(seed))

  # Num threads; default 0 => detect from system in C++.
    if (!is.numeric(n_thread) || n_thread < 0)
        stop("Invalid value for 'n_thread'.")

  # Use sparse matrix
    if (inherits(x, "dgCMatrix")) {
        sparse_x <- x
        x <- matrix(numeric())
    } else {
        sparse_x <- NULL
        if (is.data.frame(x)) x <- data.matrix(x)
        if (storage.mode(x) %in% "integer") storage.mode(x) <- "double"
    }

    y_mat <- as.matrix(as.numeric(y))

  # Call the library's train function using cpp11
    result <- cpp11_train(
      # data
        x, y_mat, sparse_x, case_weights,
      # forest
        tree_type, n_tree,
      # variables
        predictor_names, names_of_unordered,
      # sampling/bootstrapping
        replace, sample_fraction,
      # drawing candidate predictors
        n_try,
        draw_predictor_weights, names_of_always_draw,
      # splitting rules
        split_rule,
        max_depth, min_split_n_sample, min_leaf_n_sample,
      # tree-type specific
        response_weights,
      # split-rule specific
        n_random_split, alpha, min_prop,
      # miscellaneous
        seed, save_memory, n_thread,
        verbose
    )

    if (is.factor(y)) result$response_levels <- levels(y)

    if (!is.null(predictor_levels)) result$predictor_levels <- predictor_levels

    class(result) <- "literanger"

    invisible(result)

}

