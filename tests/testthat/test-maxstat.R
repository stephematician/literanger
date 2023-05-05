
test_that("regression forest accepts split-rule", {
    expect_silent(
        rf <- train(data=iris, response_name="Sepal.Length",
                    split_rule="maxstat")
    )
    expect_is(rf, "literanger")
})

test_that("regression forest has acceptable r-squared", {
    rf <- train(data=iris, response_name="Sepal.Length", split_rule="maxstat")
    expect_gt(1 - rf$oob_error / var(iris$Sepal.Length), 0.5)
})


test_that("error if training a classificaiton forest with split-rule", {
    expect_error(
        train(data=iris, response_name="Species", split_rule = "maxstat"),
        "'arg' should be one of .gini., .extratrees., .hellinger."
    )
})

test_that("error if alpha or min_prop is out of range", {
    expect_error(train(data=iris, response_name=="Sepal.Length",
                       split_rule="maxstat", alpha=-1))
    expect_error(train(data=iris, response_name=="Sepal.Length",
                       split_rule="maxstat", alpha=2))
    expect_error(train(data=iris, response_name=="Sepal.Length",
                       split_rule="maxstat", minprop=-1))
    expect_error(train(data=iris, response_name=="Sepal.Length",
                       split_rule="maxstat", minprop=1))
})

