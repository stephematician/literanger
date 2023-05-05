
test_that("classification forest with ordered predictors accepts split-rule", {
    expect_silent(rf <- train(data=iris, response_name="Species",
                              split_rule="extratrees"))
    expect_is(rf, "literanger")
    expect_lt(rf$oob_error, 0.2)
})

test_that("classification forest with ordered predictors has acceptable out-of-bag error", {
    rf <- train(data=iris, response_name="Species", split_rule="extratrees")
    expect_lt(rf$oob_error, 0.2)
})

test_that("regression forest with ordered predictors accepts split-rule", {
    expect_silent(
        rf <- train(data=iris, response_name="Sepal.Length",
                    split_rule="extratrees")
    )
    expect_is(rf, "literanger")
    expect_gt(1 - rf$oob_error / var(iris$Sepal.Length), 0.5)
})

test_that("regression forest with ordered predictors has acceptable out-of-bag error", {
    rf <- train(data=iris, response_name="Sepal.Length", split_rule="extratrees")
    expect_gt(1 - rf$oob_error / var(iris$Sepal.Length), 0.5)
})

test_that("classification forest can have large number of random splits", {
    expect_silent(
        train(data=iris, response_name="Species",
              split_rule="extratrees", n_random_split=100)
    )
})

