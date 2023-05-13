
test_that("can train forests with ordered predictors", {
    expect_silent(rf <- train(data=iris, response_name="Species",
                              split_rule="extratrees"))
    expect_silent(rf <- train(data=iris, response_name="Sepal.Length",
                              split_rule="extratrees"))
})

test_that("trained forests with ordered predictors have acceptable out-of-bag error", {
    rf_class <- train(data=iris, response_name="Species",
                      split_rule="extratrees")
    expect_lt(rf_class$oob_error, 0.2)
    rf_num <- train(data=iris, response_name="Sepal.Length",
                    split_rule="extratrees")
    expect_gt(1 - rf_num$oob_error / var(iris$Sepal.Length), 0.5)
})

test_that("forest can have large number of random splits", {
    expect_silent(train(data=iris, response_name="Species",
                        split_rule="extratrees", n_random_split=100))
    expect_silent(train(data=iris, response_name="Sepal.Length",
                        split_rule="extratrees", n_random_split=100))
})

