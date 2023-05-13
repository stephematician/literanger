dat_class <- {
    set.seed(2)
    n <- 20
    data.frame(x=sample(c("A", "B", "C", "D"), n, replace=TRUE),
               y=factor(rbinom(n, 1, 0.5)),
               stringsAsFactors=FALSE)
}
dat_num <- {
    set.seed(2)
    n <- 20
    data.frame(x=sample(c("A", "B", "C", "D"), n, replace=TRUE),
               y=rnorm(n, 1, 0.5),
               stringsAsFactors=FALSE)
}

dat_class_ftor <- dat_class
dat_class_ftor$x <- factor(dat_class_ftor$x, ordered=FALSE)

dat_num_ftor <- dat_num
dat_num_ftor$x <- factor(dat_num_ftor$x, ordered=FALSE)

test_that("can train a forest with partitioning", {
    expect_silent(rf_class <- train(data=dat_class, response_name="y",
                                    unordered_predictors="partition"))
    expect_true(length(rf_class$names_of_unordered) > 0)
    expect_silent(rf_num <- train(data=dat_num, response_name="y",
                                  unordered_predictors="partition"))
    expect_true(length(rf_num$names_of_unordered) > 0)
})

test_that("can train a forest with re-ordering via PCA score", {
    expect_silent(rf_class <- train(data=dat_class, response_name="y",
                                    unordered_predictors="order"))
    expect_true(length(rf_class$names_of_unordered) == 0)
    expect_silent(rf_num <- train(data=dat_num, response_name="y",
                                  unordered_predictors="order"))
    expect_true(length(rf_num$names_of_unordered) == 0)
})

test_that("get error when too many levels in factor for partitioning", {
    n <- 100
    dat <- data.frame(x=factor(1:100, ordered=FALSE), y=rbinom(n, 1, 0.5))
    expect_error(
        train(data=dat, response_name="y", unordered_predictors="partition"),
        "Too many levels in unordered categorical variable x",
        fixed=TRUE
    )
})

test_that("same out-of-bag error for character or unordered factor with partitioning", {
    set.seed(2)
    rf_class_char <- train(data=dat_class, response_name="y",
                           unordered_predictors="partition")
    set.seed(2)
    rf_class_ftor <- train(data=dat_class_ftor, response_name="y",
                           unordered_predictors="partition")
    expect_equal(rf_class_char$oob_error, rf_class_ftor$oob_error)

    set.seed(2)
    rf_num_char <- train(data=dat_num, response_name="y",
                         unordered_predictors="partition")
    set.seed(2)
    rf_num_ftor <- train(data=dat_num_ftor, response_name="y",
                           unordered_predictors="partition")
    expect_equal(rf_num_char$oob_error, rf_num_ftor$oob_error)
})

test_that("same out-of-bag error for character or unordered factor re-ordered by PCA score", {
    set.seed(2)
    rf_class_char <- train(data=dat_class, response_name="y",
                           unordered_predictors="order")
    set.seed(2)
    rf_class_ftor <- train(data=dat_class_ftor, response_name="y",
                            unordered_predictors="order")
    expect_equal(rf_class_char$oob_error, rf_class_ftor$oob_error)

    set.seed(2)
    rf_num_char <- train(data=dat_num, response_name="y",
                         unordered_predictors="order")
    set.seed(2)
    rf_num_ftor <- train(data=dat_num_ftor, response_name="y",
                           unordered_predictors="order")
    expect_equal(rf_num_char$oob_error, rf_num_ftor$oob_error)
})

test_that("can train forest when single-level predictor is re-ordered by PCA score", {
    n <- 20
    dat_class_one <- data.frame(x=sample(c("A"), n, replace=TRUE),
                                y=factor(sample(c("A", "B", "C", "D"),
                                         n, replace=TRUE)),
                                stringsAsFactors=FALSE)
    expect_silent(train(data=dat_class_one, response_name="y",
                        unordered_predictors="order"))

    dat_num_one <- data.frame(x=sample(c("A"), n, replace=TRUE), y=rnorm(n),
                              stringsAsFactors=FALSE)
    expect_silent(train(data=dat_num_one, response_name="y",
                        unordered_predictors="order"))
})

test_that("result same when training forests if no unordered factors", {
    set.seed(100)
    rf_class_ignr <- train(data=iris, response_name="Species",
                           unordered_predictors="ignore")
    set.seed(100)
    rf_class_ordr <- train(data=iris, response_name="Species",
                           unordered_predictors="order")
    set.seed(100)
    rf_class_part <- train(data=iris, response_name="Species",
                           unordered_predictors="partition")
    expect_equal(rf_class_ignr$oob_error, rf_class_ordr$oob_error)
    expect_equal(rf_class_ordr$oob_error, rf_class_part$oob_error)

    is.factor <- sapply(iris, is.factor)
    set.seed(100)
    rf_num_ignr <- train(data=iris[!is.factor], response_name="Sepal.Length",
                         unordered_predictors="ignore")
    set.seed(100)
    rf_num_ordr <- train(data=iris[!is.factor], response_name="Sepal.Length",
                         unordered_predictors="order")
    set.seed(100)
    rf_num_part <- train(data=iris[!is.factor], response_name="Sepal.Length",
                         unordered_predictors="partition")
    expect_equal(rf_num_ignr$oob_error, rf_num_ordr$oob_error)
    expect_equal(rf_num_ordr$oob_error, rf_num_part$oob_error)
})

test_that("can train forests with unordered predictors and 'extratrees'", {
    expect_silent(rf <- train(data=dat_class, response_name="y",
                              split_rule="extratrees",
                              unordered_predictors="partition"))
    expect_silent(rf <- train(data=iris, response_name="Sepal.Length",
                              split_rule="extratrees",
                              unordered_predictors="partition"))
})

test_that("trained forests with unordered predictors and 'extratrees' has acceptable out-of-bag error", {
    set.seed(42)
    rf_class <- train(data=iris, response_name="Species",
                      split_rule="extratrees", unordered_predictors="partition")
    expect_lt(rf_class$oob_error, 0.2)
    set.seed(42)
    rf_num <- train(data=iris, response_name="Sepal.Length",
                    split_rule="extratrees", unordered_predictors="partition")
    expect_gt(1 - rf_num$oob_error / var(iris$Sepal.Length), 0.5)
})

test_that("maximally selected rank statistics metric fails with partitioning", {
    expect_error(
        train(data=iris, response_name="Sepal.Length", split_rule="maxstat",
              unordered_predictors="partition"),
        paste("Unordered factor splitting not implemented for 'maxstat' or",
              "'beta' splitting rule."),
        fixed=TRUE
    )
})

test_that("can predict unobserved levels given unordered predictors", {
    set.seed(1)
    n <- 20
    dat_train <- data.frame(x1=sample(c("A", "B", "C"), n, replace=TRUE),
                            x2=sample(c("A", "B", "C"), n, replace=TRUE),
                            y=rbinom(n, 1, 0.5),
                            stringsAsFactors=FALSE)

    dat_test <- data.frame(x1=sample(c("A", "B", "C", "D"), n, replace=TRUE),
                           x2=sample(c("A", "B", "C", "D"), n, replace=TRUE),
                           stringsAsFactors=FALSE)

    rf_ignr <- train(data=dat_train, response_name="y",
                     unordered_predictors="ignore")
    expect_warning(
        predict(rf_ignr, newdata=dat_test),
        "Predictor levels found that were not present during training",
        fixed=TRUE
    )

    rf_part <- train(data=dat_train, response_name="y",
                     unordered_predictors="partition")
    expect_warning(
        predict(rf_part, newdata=dat_test),
        "Predictor levels found that were not present during training",
        fixed=TRUE
    )

    rf_ordr <- train(data=dat_train, response_name="y",
                     unordered_predictors="order")
    expect_warning(
        predict(rf_ordr, newdata=dat_test),
        "Predictor levels found that were not present during training",
        fixed=TRUE
    )
})

test_that("warning for maximally selected rank statistics with re-ordering via PCA score", {
    expect_warning(
        train(data=iris, response_name="Sepal.Length", split_rule="maxstat",
              unordered_predictors="order"),
        paste("The 'order' mode for unordered factor handling with the",
              "'maxstat' splitrule is experimental."),
        fixed=TRUE
    )
})

test_that("can train on NA factor levels when using re-ordering via PCA score", {
    dat <- data.frame(x=addNA(factor(c("a", "a", NA, NA, "b", "b"))),
                     y=c(1, 2, 3, 4, 5, 6))
    expect_silent(
        train(data=dat, response_name="y", unordered_predictors="order")
    )
})

test_that("can use re-ordering via PCA score when numerics in data", {
    n <- 20

    dat_bin <- data.frame(x1=sample(c("A", "B", "C"), n, replace=TRUE),
                          x2=sample(1:3, n, replace=TRUE),
                          y=factor(sample(c("A", "B"), n, replace=TRUE)),
                          stringsAsFactors=FALSE)
    expect_silent(
        rf_bin <- train(data=dat_bin, response_name="y",
                        unordered_predictors="order")
    )
    expect_silent(predict(rf_bin, newdata=dat_bin))

  # Multiclass classification
    dat_fac <- dat_bin
    dat_fac$y <- factor(sample(c("A", "B", "C", "D"), n, replace=TRUE))
    expect_silent(
        rf_fac <- train(data=dat_fac, response_name="y",
                        unordered_predictors="order")
    )
    expect_silent(predict(rf_fac, newdata=dat_fac))

  # Regression
    dat_num <- dat_class
    dat_num$y <- rnorm(n)
    expect_silent(
        rf_num <- train(data=dat_num, response_name="y",
                        unordered_predictors="order")
    )
    expect_silent(predict(rf_num, newdata=dat_num))
})

test_that("can use partitioning with a large number of levels", {
    n <- 43
    dat <- data.frame(x=factor(1:n, ordered=FALSE),  y=rbinom(n, 1, 0.5))

    expect_silent(
        rf <- train(data=dat, response_name="y", split_rule="extratrees")
    )
  #  max_split <- max(sapply(1:rf$n_tree, function(i) {
  #      max(log2(rf$forest$split.values[[i]]), na.rm=TRUE)
  #  }))
  #  expect_lte(max_split, n)
})

