
test_that("classification forest accepts split-rule", {
    expect_silent(
        rf <- train(data=droplevels(iris[1:100,]), response_name="Species",
                    split_rule="hellinger")
    )
    expect_is(rf, "literanger")
})

test_that("classification forest has acceptable out-of-bag error", {
    rf <- train(data=droplevels(iris[1:100,]), response_name="Species",
                split_rule="hellinger")
    expect_lt(rf$oob_error, 0.3)
})

test_that("classification forest with non-factor response accepts split-rule", {
    dat <- iris[1:100, ]
    dat$Species <- as.numeric(dat$Species)
    expect_silent(
        rf <- train(data=dat, response_name="Species",
                    classification=T, split_rule="hellinger")
    )
    expect_is(rf, "literanger")
    expect_lt(rf$oob_error, 0.3)
})

test_that("error if response data has more than two classes", {
    expect_error(
        train(data=iris, response_name="Species", split_rule="hellinger"), 
        "Hellinger split metric only implemented for binary classification.",
        fixed=T
    )
})

test_that("error if training a regression forest with split-rule", {
    expect_error(
        train(data=iris, response_name="Sepal.Length",
              split_rule="hellinger"),
        "'arg' should be one of .variance., .extratrees., .maxstat., .beta."
    )
})

