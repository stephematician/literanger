
test_that("classification foreset with unordered predictors accepts split-rule", {
    n <- 20
    dat <- data.frame(x=sample(c("A", "B", "C", "D"), n, replace=T), 
                      y=factor(rbinom(n, 1, 0.5)), 
                      stringsAsFactors=F)
    expect_silent(
        rf <- train(data=dat, response_name="y", min_split_n_sample=n / 2, 
                    split_rule="extratrees",
                    unordered_predictors="partition")
    )
    expect_is(rf, "literanger")
})

test_that("regression forest with unordered predictors accepts split-rule", {
    expect_silent(
        rf <- train(data=iris, response_name="Sepal.Length",
                    split_rule="extratrees",
                    unordered_predictors="partition")
    )
    expect_is(rf, "literanger")
})

test_that("regression forest with unordered predictors has acceptable out-of-bag error", {
    rf <- train(data=iris, response_name="Sepal.Length",
                split_rule="extratrees", unordered_predictors="partition")
    expect_gt(1 - rf$oob_error / var(iris$Sepal.Length), 0.5)
})
