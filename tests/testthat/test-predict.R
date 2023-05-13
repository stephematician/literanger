
test_that("prediction has acceptable accuracy for classification forest", {
    rf <- train(data=iris, response_name="Species")
    pred <- predict(rf, newdata=iris)
    expect_gt(mean(iris$Species == pred$values), 0.9)
})

test_that("can use case weights", {
    expect_silent(
        train(data=iris, response_name="Species",
              case_weights = rep(1, nrow(iris)))
    )
  ## Should only predict setosa now
    rf <- train(data=iris, response_name="Species",
                case_weights=c(rep(1, 50), rep(0, 100)))
    pred <- predict(rf, newdata=iris)$values
    expect_true(all(pred == "setosa"))
})

test_that("result is independent of column ordering", {
    dat <- iris[, c(1:2, 5, 3:4)]
    rf <- train(data=dat, response_name="Species")
    set.seed(42)
    pred1 <- predict(rf, newdata=iris)$values
    set.seed(42)
    pred2 <- predict(rf, newdata=dat)$values
    expect_equal(pred1, pred2)
    expect_gte(mean(pred2 == dat$Species), 0.9)
  ## No response column
    expect_gte(mean(predict(rf, newdata=dat[, -3])$values == dat$Species), 0.9)
})

test_that("meaningful predictions with max_depth = 1", {
    rf <- train(data=iris, response_name="Sepal.Length", max_depth=1)
    pred <- predict(rf, newdata=iris)$values
    expect_gte(min(pred), min(iris$Sepal.Length))
    expect_lte(max(pred), max(iris$Sepal.Length))
})

test_that("error if predictors have missing values", {
    rf <- train(data=iris, response_name="Species")
    dat <- iris
    dat[4, 4] <- dat[25,1] <- NA
    expect_error(
        predict(rf, newdata=dat),
        "Missing values in the predictors: Sepal.Length, Petal.Width."
    )
})

test_that("error if unknown value for prediction type", {
    rf <- train(data=iris, response_name="Species")
    expect_error(predict(rf, newdata=iris, prediction_type="foo"))
})

test_that("missing value allowed in discarded column", {
    rf <- train(data=iris, response_name="Species",
                predictor_names="Sepal.Length")
    dat <- iris
    dat[1, "Sepal.Width"] <- NA
    expect_silent(predict(rf, newdata=dat))
})
