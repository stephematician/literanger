
set.seed(42)
rf_iris_df <- train(data=iris, response_name="Species")
set.seed(42)
rf_iris_mat <- train(data=data.matrix(iris), response_name="Species")

test_that("literanger object returned by train", {
    expect_is(rf_iris_df, "literanger")
    expect_is(rf_iris_mat, "literanger")
})

test_that("trained forest has default number of trees", {
    expect_equal(rf_iris_df$n_tree, formals(literanger::train)$n_tree)
    expect_equal(rf_iris_mat$n_tree, formals(literanger::train)$n_tree)
})

test_that("'predictor_names' does not include the response", {
    expect_equal(setdiff(names(iris), rf_iris_df$predictor_names), "Species")
    expect_equal(setdiff(names(iris), rf_iris_mat$predictor_names), "Species")
})

test_that("out-of-bag error estimate is number", {
    expect_true(hasName(rf_iris_df, 'oob_error'))
    expect_false(is.na(rf_iris_df$oob_error))
    expect_true(hasName(rf_iris_mat, 'oob_error'))
    expect_false(is.na(rf_iris_mat$oob_error))
})

test_that("can train when 'data' argument has more than one class", {
    dat <- iris
    class(dat) <- c("data.frame", "data.table")
    expect_silent(train(data=dat, response_name="Species"))
})

test_that("can train when response has attributes other than names", {
    iris2 <- iris
    attr(iris2$Sepal.Width, "aaa") <- "bbb"
    expect_silent(train(data=iris2, response_name="Sepal.Width"))
})

test_that("get error when predictor has missing values", {
    dat <- iris
    dat[25, "Sepal.Length"] <- NA
    expect_error(train(data=dat, response_name="Species"),
                 "Missing values in the predictors: Sepal.Length.", fixed=TRUE)
    dat_num <- iris
    dat_num[25, "Species"] <- NA
    expect_error(train(data=dat_num, response_name="Sepal.Length"),
                 "Missing values in the predictors: Species.", fixed=TRUE)
})

test_that("get error when response has missing values", {
    dat <- iris
    dat[4, "Species"] <- NA
    expect_error(train(data=dat, response_name="Species"),
                 "Missing values in the response.", fixed=TRUE)
    dat_num <- iris
    dat_num[4,"Sepal.Length"] <- NA_real_
    expect_error(train(data=dat_num, response_name="Sepal.Length"),
                 "Missing values in the response.", fixed=TRUE)
})

test_that("missing value allowed in discarded column", {
    dat <- iris
    dat[1, "Sepal.Width"] <- NA
    expect_silent(train(data=dat, response_name="Species",
                        predictor_names="Sepal.Length"))
})

test_that("can train when a predictor is named 'none'", {
    dat <- data.frame(y=rbinom(100, 1, .5),
                      x=rbinom(100, 1, .5),
                      none=rbinom(100, 1, .5))
    rf_df <- train(data=dat, response_name="y")
    rf_mat <- train(data=data.matrix(dat), response_name="y")
    expect_equal(rf_df$predictor_names, c("x", "none"))
    expect_equal(rf_mat$predictor_names, c("x", "none"))
})

test_that("get error when sample fraction is outside (0,1]", {
    expect_error(train(data=iris, response_name="Species", sample_fraction=0))
    expect_error(train(data=iris, response_name="Species",
                       sample_fraction=1 * (1 + .Machine$double.eps)))
})

test_that("can use function as 'n_try' argument ", {
    rf <- train(data=iris, response_name="Species", n_try=function(n) n - 1)
    expect_equal(3, rf$n_try)
})

test_that("get error when 'n_try' function has error", {
    expect_error(
        train(data=iris, response_name="Species",
              n_try=function(n) stop("this is some error")),
        "'n_try' function evaluation resulted in an error.",
        fixed=TRUE)
})

if (requireNamespace("tibble", quietly=TRUE)) tb <- tibble::as_tibble(iris)

test_that("can train with tibble", {
    skip_if_not_installed("tibble")
    expect_silent(train(data=tb, response_name="Species"))
    expect_silent(train(data=tb, response_name="Sepal.Length"))
})

