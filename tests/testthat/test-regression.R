## Tests for random forests for regression

## Initialize the random forest for regression
dat <- data.matrix(iris)
rg.reg <- train(data=iris, response_name="Sepal.Length")
rg.mat <- train(data=dat, response_name="Sepal.Length")

## Basic tests (for all random forests equal)
test_that("training result has correct class and size", {
    expect_is(rg.reg, "literanger")
    expect_equal(length(rg.reg), 14)
})

test_that("prediction is a numeric vector", {
    pred <- predict(rg.reg, newdata=iris)
    expect_is(pred$values, "numeric")
    expect_null(dim(pred$values))
})

test_that("trained forest has default number of trees", {
    expect_equal(rg.reg$n_tree, formals(literanger::train)$n_tree)
})

test_that("trained forest has correct predictor names", {
    expect_equal(setdiff(names(iris), rg.reg$predictor_names), "Sepal.Length")
})

test_that("can use matrix interface for training", {
    expect_equal(rg.mat$tree_type, "regression")
    expect_equal(setdiff(names(iris), rg.mat$predictor_names), "Sepal.Length")
})

test_that("can use matrix interface for prediction", {
    expect_silent(predict(rg.mat, newdata=dat, prediction_type="bagged"))
})

test_that("can use save_memory option when training", {
    rf <- train(data=iris, response_name="Sepal.Length", save_memory=T)
    expect_equal(rf$tree_type, "regression")
})

test_that("can omit response variable during prediction", {
    n <- 50

    dt <- data.frame(x1=runif(n), x2=runif(n), y=rbinom(n, 1, 0.5))
    rf <- train(data=dt, response_name="y")
    expect_silent(predict(rf, newdata=dt))
    expect_silent(predict(rf, newdata=dt[, 1:2]))

    dt2 <- data.frame(y=rbinom(n, 1, 0.5), x1=runif(n), x2=runif(n))
    rf <- train(data=dt2, response_name="y")
    expect_silent(predict(rf, newdata=dt2))
    expect_silent(predict(rf, newdata=dt2[, 2:3]))
})

test_that("results not all the same", {
    n <- 50

    dt <- data.frame(x=runif(n), y=rbinom(n, 1, 0.5))
    rf <- train(data=dt, response_name="y")
    expect_gt(diff(range(predict(rf, newdata=dt)$values)), 0)
    expect_gt(diff(range(predict(rf, newdata=dt[, 1, drop = FALSE])$values)), 0)
})


## Splitrule
test_that("'variance' is default splitrule", {
    set.seed(42)
    rf1 <- train(data=iris, response_name="Sepal.Length", n_tree=5)

    set.seed(42)
    rf2 <- train(data=iris, response_name="Sepal.Length", n_tree=5,
                 split_rule="variance")

    expect_equal(rf1$split_rule, "variance")
    expect_equal(rf2$split_rule, "variance")
    expect_equal(rf1$oob_error, rf2$oob_error)
})

