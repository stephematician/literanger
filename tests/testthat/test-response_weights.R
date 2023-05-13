
test_that("can use response weights", {
    expect_silent(
        train(data=iris, response_name="Species",
              response_weights=c(0.5, 1, 0.1))
    )
})

test_that("get error when training regression forest", {
    expect_error(
        train(data=iris, response_name="Sepal.Length",
              response_weights=c(0.5, 1, 0.1)),
        "'response_weights' argument only valid for classification forests",
        fixed=TRUE
    )
})

test_that("prediction accuracy for minority class increases with higher weight", {
    set.seed(3)
    n <- 100
    x <- rnorm(n)
    beta0 <- -3
    beta <- 1
    y <- factor(rbinom(n, 1, plogis(beta0 + beta * x)))
    dat <- data.frame(y=y, x)
    is_minor <- dat$y == "1"

    set.seed(42)
    rf <- train(data=dat, response_name="y", min_split_n_sample=50,
                response_weights=c(1, 1))
    pred <- predict(rf, newdata=dat[is_minor,])
    accuracy_minor <- mean(pred$values == 1, na.rm=TRUE)

    set.seed(42)
    rf_wtd <- train(data=dat, response_name="y", min_split_n_sample=50,
                    response_weights=c(0.01, 0.99))
    pred_wtd <- predict(rf_wtd, newdata=dat[is_minor,])
    accuracy_minor_wtd <- mean(pred_wtd$values == 1, na.rm=TRUE)

    expect_gt(accuracy_minor_wtd, accuracy_minor)
})

test_that("get error when response weights of wrong size", {
    expect_error(
        train(data=iris, response_name="Species", response_weights=c(0.5, 1)),
        "Number of response weights not equal to number of classes.",
        fixed=TRUE
    )
})

test_that("get error when response weights are NA", {
    expect_error(
        train(data=iris, response_name="Species",
              response_weights=c(0.5, 1, NA)),
        "missing value where TRUE/FALSE needed",
        fixed=TRUE
    )
})

test_that("get error when response weights are not numeric", {
    expect_error(
        train(data=iris, response_name="Species",
              response_weights = c(0.5, 1, "a")),
        paste("Invalid value for 'response_weights'. Please give a vector of",
              "non-negative values."),
        fixed=TRUE
    )
})

