
# Generate data with 0..1 outcome
set.seed(42)
n <- 100
p <- 4
beta <- c(0, 1, 2, 3)
x <- replicate(p, runif(n))
y <- as.vector(x %*% beta)
min_y <- as.numeric((beta < 0) %*% beta)
max_y <- as.numeric((beta > 0) %*% beta)
y <- (y - min_y) / (max_y - min_y)
dat <- data.frame(y = y, x)

test_that("can train a regression forest", {
    set.seed(42)
    rf <- train(data=dat, response_name="y", split_rule="beta")
    expect_is(rf, "literanger")
})

test_that("get error when training a classification forest", {
    expect_error(
        train(data=iris, response_name="Species", split_rule="beta"),
        "'arg' should be one of .gini., .extratrees., .hellinger."
    )
})

test_that("get error when outcomes outside (0,1)", {
    expect_error(
        train(data=iris, response_name="Sepal.Length", split_rule="beta"),
        paste("Beta log-likelihood metric applicable to regression data in the",
              "interval (0,1)."),
        fixed=TRUE
    )
})

test_that("trained forest has acceptable r-squared", {
    set.seed(42)
    rf <- train(data=dat, response_name="y", split_rule="beta")
    # expect_lt(rf$oob_error, 0.2)
    expect_gt(1 - rf$oob_error / var(dat$y), 0.65)
})

