
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
