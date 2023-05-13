
test_that("can train a regression forest", {
    expect_silent(
        rf <- train(data=iris, response_name="Sepal.Length",
                    split_rule="maxstat")
    )
    expect_is(rf, "literanger")
})

test_that("get error when training a classification forest", {
    expect_error(
        train(data=iris, response_name="Species", split_rule = "maxstat"),
        "'arg' should be one of .gini., .extratrees., .hellinger."
    )
})

test_that("get error when alpha is outside (0,1]", {
    expect_error(train(data=iris, response_name=="Sepal.Length",
                       split_rule="maxstat", alpha=0))
    expect_error(train(data=iris, response_name=="Sepal.Length",
                       split_rule="maxstat",
                       alpha=0.5 * (1 + .Machine$double.eps)))
})

test_that("get error when minprop is outside [0,0.5]", {
    expect_error(train(data=iris, response_name=="Sepal.Length",
                       split_rule="maxstat",
                       minprop=-.Machine$double.xmin))
    expect_error(train(data=iris, response_name=="Sepal.Length",
                       split_rule="maxstat",
                       minprop=0.5 * (1 + .Machine$double.eps)))
})

test_that("trained forest has acceptable r-squared", {
    rf <- train(data=iris, response_name="Sepal.Length", split_rule="maxstat")
    expect_gt(1 - rf$oob_error / var(iris$Sepal.Length), 0.5)
})

