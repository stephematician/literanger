
test_that("silent when data argument has more than one class", {
    dat <- iris
    class(dat) <- c("data.frame", "data.table")
    expect_silent(train(data=dat, response_name="Species"))
})

test_that("error if sample fraction is 0 or >1", {
    expect_error(train(data=iris, response_name="Species", sample_fraction=0))
    expect_error(train(data=iris, response_name="Species", sample_fraction=2))
})

test_that("error if sample fraction is vector for regression", {
    expect_error(
        train(data=iris, response_name="Sepal.Length",
              sample_fraction=c(0.1, 0.2)), 
        paste("Invalid value for 'sample_fraction'. Vector values only valid",
              "for classification forests."),
        fixed=T
    )
})

test_that("error if sample fraction is vector of wrong size", {
    expect_error(
        train(data=iris, response_name="Species", sample_fraction=c(0.1, 0.2)), 
        paste("Invalid value for 'sample_fraction'. Expecting 3 values,",
              "provided 2."),
        fixed=T
    )
})

test_that("error if element of sample fraction vector is <0 or >1", {
    expect_error(
        train(data=iris, response_name="Species",
              sample_fraction=c(0.1, 1.1, 0.3)), 
        paste("Invalid value for 'sample_fraction'. Please give a value in",
              "(0,1] or a vector of values in [0,1]."),
        fixed=T
    )
    expect_error(
        train(data=iris, response_name="Species",
              sample_fraction = c(-3, 0.5, 0.3)), 
        paste("Invalid value for 'sample_fraction'. Please give a value in",
              "(0,1] or a vector of values in [0,1]."),
        fixed=T
    )
})

test_that("error if sum of sample fraction vector is 0", {
    expect_error(
        train(data=iris, response_name="Species", sample_fraction=rep(0, 3)), 
        paste("'sample_fraction' too small (results in zero samples)."),
        fixed=T
    )
})

test_that("error if sample without replacement and not enough samples", {
    expect_error(
        train(data=iris, response_name="Species",
              sample_fraction=c(0.2, 0.3, 0.4), replace=F), 
        paste("Not enough samples in class virginica; available: 50,",
              "requested: 60."),
        fixed=T
    )
    expect_silent(train(data=iris, response_name="Species",
                        sample_fraction=c(0.2, 0.3, 0.4), replace=T))
})

test_that("error if sample fraction and case weights", {
    expect_error(
        train(data=iris, response_name="Species",
              sample_fraction=c(0.2, 0.3, 0.4),
              case_weights=rbinom(nrow(iris), 1, 0.5)), 
        paste("Combination of 'case_weights' argument and class-wise sampling",
              "not supported."),
        fixed=T
    )
})

test_that("error if missing values in predictors", {
    dat <- iris
    dat[25, 1] <- NA
    expect_error(train(data=dat, response_name="Species"),
                 "Missing values in the predictors: Sepal.Length.", fixed=T)
    dat <- iris
    dat[4, 5] <- NA
    expect_error(train(data=dat, response_name="Species"),
                "Missing values in the response.", fixed=T)
})

test_that("missing value allowed in discarded column", {
    dat <- iris
    dat[1, "Sepal.Width"] <- NA
    expect_silent(train(data=dat, response_name="Species",
                        predictor_names="Sepal.Length"))
})

test_that("no error if variable named forest", {
    dat <- iris
    dat$forest <- rnorm(150)
    rf <- train(data=dat, response_name="Species")
    expect_silent(predict(rf, newdata=dat))
})

test_that("out-of-bag error estimate is number", {
    rf <- train(data=iris, response_name="Species")
    expect_false(is.na(rf$oob_error))
})

test_that("does not crash when variable named 'none'", {
    dat <- data.frame(y=rbinom(100, 1, .5), 
                      x=rbinom(100, 1, .5), 
                      none=rbinom(100, 1, .5))
    rf <- train(data=dat, response_name="y")
    expect_equal(rf$predictor_names, c("x", "none"))
    expect_silent(predict(rf, newdata=dat))
})

test_that("n_try function input works as expected", {
    rf <- train(data=iris, response_name="Species", n_try=function(n) n - 1)
    expect_equal(3, rf$n_try)
})

test_that("n_try function error halts the train function", {
    expect_error(train(data=iris, response_name="Species",
                       n_try=function(n) stop("this is some error")), 
                 "'n_try' function evaluation resulted in an error.", fixed=T)
})

