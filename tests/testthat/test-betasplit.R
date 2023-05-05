# Generate data with 0..1 outcome
n <- 100
p <- 4
beta <- c(0, 1, 2, 3)
x <- replicate(p, runif(n))
y <- as.vector(x %*% beta)
y <- (y-min(y))/(max(y)-min(y))
dat <- data.frame(y = y, x)

test_that("regression forest accepts beta log-likelihood metric", {
    rf <- train(data=dat, response_name="y", split_rule="beta", seed=1)
    expect_is(rf, "literanger")
})

test_that("regression forest with beta log-likelihood metric has acceptable out-of-bag error", {
    rf <- train(data=dat, response_name="y", split_rule="beta", seed=1)
    expect_lt(rf$oob_error, 0.2)
})

test_that("error if outcomes outside [0,1] combined with beta log-likelihood metric", {
    expect_error(
        train(data=iris, response_name="Species", split_rule="beta"),
        "'arg' should be one of .gini., .extratrees., .hellinger.",
    )
})

