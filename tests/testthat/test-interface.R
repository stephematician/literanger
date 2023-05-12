
test_that("supports dependent variable has attributes other than names", {
    iris2 <- iris
    attr(iris2$Sepal.Width, "aaa") <- "bbb"
    expect_silent(train(data=iris2, response_name="Sepal.Width"))
})

test_that("supports dependent variable is matrix with one column", {
    iris2 <- iris
    iris2$Sepal.Width = scale(iris$Sepal.Width)
    expect_silent(train(data=iris2, response_name="Sepal.Width"))
})

test_that("result same with x/y interface, classification", {
    set.seed(300)
    rf <- train(data=iris, response_name="Species")
    pred <- predict(rf, newdata=iris)

    set.seed(300)
    rf_xy <- train(y=iris[, 5], x=iris[, -5])
    pred_xy <- predict(rf, newdata=iris[, -5])

    expect_equal(rf$oob_error, rf_xy$oob_error)
    expect_equal(pred$values, pred_xy$values)
})

test_that("result same with x/y interface, regression", {
    set.seed(300)
    rf <- train(data=iris, response_name="Sepal.Length")
    pred <- predict(rf, newdata=iris)

    set.seed(300)
    expect_silent(rf_xy <- train(y=iris[, 1], x=iris[, -1]))
    expect_silent(pred_xy <- predict(rf_xy, newdata=iris[, -1]))

    expect_equal(rf$oob_error, rf_xy$oob_error)
    expect_equal(pred$values, pred_xy$values)
})

test_that("column order does not change prediction", {
    dat <- iris[, c(sample(1:4), 5)]
    rf <- train(data=iris, response_name="Species")

    set.seed(42)
    pred1 <- predict(rf, iris)$values

    set.seed(42)
    pred2 <- predict(rf, dat)$values

    expect_equal(pred1, pred2)
})


if (requireNamespace("tibble", quietly=TRUE)) tb <- tibble::as_tibble(iris)

test_that("training works with tibbles", {
    skip_if_not_installed("tibble")
    rf1 <- train(data=tb, response_name="Species", seed=1)
    pred1tb <- predict(rf1, newdata=tb, seed=2)
    pred1 <- predict(rf1, newdata=iris, seed=2)

    rf2 <- train(data=iris, response_name="Species", seed=1)
    pred2tb <- predict(rf2, newdata=tb, seed=2)
    pred2 <- predict(rf2, newdata=iris, seed=2)

    expect_equal(rf1$oob_error, rf2$oob_error)
    expect_equal(pred1, pred1tb)
    expect_equal(pred2, pred2tb)
    expect_equal(pred1tb, pred2tb)
})

