## Tests for random forests for classification

## Initialize the random forest for classification
dat <- data.matrix(iris)
rg.class <- train(data=iris, response_name="Species")
rg.mat <- train(data=dat, response_name="Species", classification=T)

## Basic tests (for all random forests equal)
test_that("training result has correct class and size", {
    expect_is(rg.class, "literanger")
    expect_equal(length(rg.class), 15)
})

test_that("prediction is a factor", {
    expect_silent(pred <- predict(rg.class, newdata=iris))
    expect_is(pred$values, "factor")
    expect_null(dim(pred$values))
})

test_that("trained forest has default number of trees", {
    expect_equal(rg.class$n_tree, formals(literanger::train)$n_tree)
})

test_that("trained forest has correct predictor names", {
    expect_equal(setdiff(names(iris), rg.class$predictor_names), "Species")
})

test_that("can use matrix interface for training", {
    expect_equal(rg.mat$tree_type, "classification")
    expect_equal(setdiff(names(iris), rg.mat$predictor_names), "Species")
})

test_that("can use matrix interface for prediction", {
    expect_silent(predict(rg.mat, newdata=dat, prediction_type="bagged"))
    # expect_output(predict(rg.mat, newdata=dat, prediction_type="inbag"),
    #               regexp=predict_out_regex)
})

test_that("can use save_memory option when training", {
    expect_silent(
        rf <- train(data=iris, response_name="Species", save_memory=T)
    )
    expect_equal(rf$tree_type, "classification")
})

test_that("can omit response variable during prediction", {
    n <- 50

    dt <- data.frame(x1=runif(n), x2=runif(n), y=factor(rbinom(n, 1, 0.5)))
    rf <- train(data=dt, response_name="y")
    expect_silent(predict(rf, newdata=dt))
    expect_silent(predict(rf, newdata=dt[, 1:2]))

    dt2 <- data.frame(y=factor(rbinom(n, 1, 0.5)), x1=runif(n), x2=runif(n))
    rf <- train(data=dt2, response_name="y")
    expect_silent(predict(rf, dt2))
    expect_silent(predict(rf, dt2[, 2:3]))
})

## Special tests for random forests for classification
test_that("prediction works for single observations", {
    expect_silent(pred <- predict(rg.class, newdata=head(iris, 1)))
    expect_equal(pred$values, iris[1,"Species"])
})

## Splitrule
test_that("error if variance split-rule", {
    expect_error(train(data=iris, response_name="Species",
                        split_rule="variance"))
})

test_that("'gini' is default split-rule", {
    set.seed(42)
    rf1 <- train(data=iris, response_name="Species")

    set.seed(42)
    rf2 <- train(data=iris, response_name="Species", split_rule="gini")

    expect_equal(rf1$split_rule, "gini")
    expect_equal(rf2$split_rule, "gini")
    expect_equal(rf1$oob_error, rf2$oob_error)
})

test_that("results using extratress split-rule differ from gini", {
    set.seed(42)
    rf1 <- train(data=iris, response_name="Species", split_rule="extratrees")

    set.seed(42)
    rf2 <- train(data=iris, response_name="Species", split_rule="gini")

    expect_equal(rf1$split_rule, "extratrees")
    expect_equal(rf2$split_rule, "gini")
    expect_false(rf1$oob_error == rf2$oob_error)
})

test_that("training with numerically near-identical splits succeeds", {
    dat <- data.frame(a = factor(1:2),
                      z = c(1.7629414498915687570246291215880773,
                            1.7629414498915689790692340466193854))
    expect_silent(train(data=dat, response_name="a", n_thread=1, n_tree=1))
})

test_that("warning given if unused factor levels in response", {
    expect_warning(
        rf <- train(data=iris[1:100, ], response_name="Species"),
        "Dropped unused factor level(s) in response variable: virginica.",
        fixed=T
    )
    pred <- predict(rf, newdata=iris)
    expect_equal(levels(pred$values), levels(iris$Species))
})

test_that("predictions with unused factor levels provided", {
    expect_warning(
        rf <- train(data=iris[51:150, ], response_name="Species"),
        "Dropped unused factor level(s) in response variable: setosa.",
        fixed=T
    )
    pred <- predict(rf, newdata=iris)
    expect_equal(sum(is.na(pred$values)), 0)
})

test_that("logical response is converted to factor", {
    dat <- iris
    dat[["Species"]] <- dat[["Species"]] == "setosa"
    rf <- train(data=dat, response_name="Species")
    pred <- predict(rf, newdata=iris)
    expect_is(pred$values, "numeric")
    expect_null(dim(pred$values))
})

