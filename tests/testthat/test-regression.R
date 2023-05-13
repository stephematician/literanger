## Tests for random forests for regression

## Initialize the random forest for regression
iris_mat <- data.matrix(iris)
set.seed(42)
rf_df <- train(data=iris, response_name="Sepal.Length")
set.seed(42)
rf_mat <- train(data=iris_mat, response_name="Sepal.Length")

## Basic tests (for all random forests equal)
test_that("training result has correct class and size", {
    expect_is(rf_df, "literanger")
    expect_equal(length(rf_df), 14)
})

test_that("prediction is a numeric vector", {
    pred <- predict(rf_df, newdata=iris)
    expect_is(pred$values, "numeric")
    expect_null(dim(pred$values))
})

test_that("error if sample fraction is vector for regression", {
    expect_error(
        train(data=iris, response_name="Sepal.Length",
              sample_fraction=c(0.1, 0.2)),
        paste("Invalid value for 'sample_fraction'. Vector values only valid",
              "for classification forests."),
        fixed=TRUE
    )
})


test_that("trained forest has default number of trees", {
    expect_equal(rf_df$n_tree, formals(literanger::train)$n_tree)
})

test_that("trained forest has correct predictor names", {
    expect_equal(setdiff(names(iris), rf_df$predictor_names), "Sepal.Length")
})

test_that("can use matrix interface for training", {
    expect_equal(rf_mat$tree_type, "regression")
    expect_equal(setdiff(names(iris), rf_mat$predictor_names), "Sepal.Length")
})

test_that("can use matrix interface for prediction", {
    expect_silent(predict(rf_mat, newdata=iris_mat, prediction_type="bagged"))
})

test_that("can use save_memory option when training", {
    rf <- train(data=iris, response_name="Sepal.Length", save_memory=TRUE)
    expect_equal(rf$tree_type, "regression")
})

test_that("can omit response variable during prediction", {
    n <- 50
    dat1 <- data.frame(x1=runif(n), x2=runif(n), y=rbinom(n, 1, 0.5))
    rf1 <- train(data=dat1, response_name="y")

    expect_silent(predict(rf1, newdata=dat1))
    expect_silent(predict(rf1, newdata=dat1[, 1:2]))

    dat2 <- data.frame(y=rbinom(n, 1, 0.5), x1=runif(n), x2=runif(n))
    rf2 <- train(data=dat2, response_name="y")

    expect_silent(predict(rf2, newdata=dat2))
    expect_silent(predict(rf2, newdata=dat2[, 2:3]))
})

test_that("predicted values not all the same", {
    n <- 50
    dat <- data.frame(x=runif(n), y=rbinom(n, 1, 0.5))
    rf <- train(data=dat, response_name="y")

    expect_gt(diff(range(predict(rf, newdata=dat)$values)), 0)
    expect_gt(diff(range(predict(rf, newdata=dat[, 1, drop=FALSE])$values)), 0)
})


## Splitrule
test_that("default split metric is 'variance'", {
    set.seed(42)
    rf_var <- train(data=iris, response_name="Sepal.Length",
                    split_rule="variance")

    expect_equal(rf_df$split_rule, "variance")
    expect_equal(rf_var$split_rule, "variance")
    expect_equal(rf_df$oob_error, rf_var$oob_error)
})

