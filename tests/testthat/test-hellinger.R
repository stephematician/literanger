
test_that("can train a classification forest", {
    expect_silent(
        rf <- train(data=droplevels(iris[1:100,]), response_name="Species",
                    split_rule="hellinger")
    )
    expect_is(rf, "literanger")
})

test_that("get error when training a regression forest", {
    expect_error(
        train(data=iris, response_name="Sepal.Length",
              split_rule="hellinger"),
        "'arg' should be one of .variance., .extratrees., .maxstat., .beta."
    )
})

test_that("can train a forest with non-factor response", {
    dat <- iris[1:100, ]
    dat$Species <- as.numeric(dat$Species)
    expect_silent(
        rf <- train(data=dat, response_name="Species",
                    classification=TRUE, split_rule="hellinger")
    )
    expect_is(rf, "literanger")
    expect_lt(rf$oob_error, 0.3)
})

test_that("get error when training with non-binary outcomes", {
    expect_error(
        train(data=iris, response_name="Species", split_rule="hellinger"),
        "Hellinger split metric only implemented for binary classification.",
        fixed=TRUE
    )
})

test_that("trained forest has acceptable out-of-bag error", {
    rf <- train(data=droplevels(iris[1:100,]), response_name="Species",
                split_rule="hellinger")
    expect_lt(rf$oob_error, 0.3)
})

