dt.class <- {
    set.seed(2)
    n <- 20
    data.frame(x=sample(c("A", "B", "C", "D"), n, replace=T),
               y=factor(rbinom(n, 1, 0.5)),
               stringsAsFactors=F)
}
dt.num <- {
    set.seed(2)
    n <- 20
    data.frame(x=sample(c("A", "B", "C", "D"), n, replace=T),
               y=rnorm(n, 1, 0.5),
               stringsAsFactors=F)
}

test_that("can train a classification forest with partitioning", {
    rf <- train(data=dt.class, response_name="y",
                unordered_predictors="partition")
    expect_true(length(rf$names_of_unordered) > 0)
})

test_that("can train a classification forest with re-ordering via PCA score", {
    rf <- train(data=dt.class, response_name="y",
                unordered_predictors="order")
    expect_true(length(rf$names_of_unordered) == 0)
})

test_that("can train a regression forest with partitioning", {
    rf <- train(data=dt.num, response_name="y",
                unordered_predictors="partition")
    expect_true(length(rf$names_of_unordered) > 0)
})

test_that("can train a regression forest with re-ordering via PCA score", {
    rf <- train(data=dt.num, response_name="y",
                unordered_predictors="order")
    expect_true(length(rf$names_of_unordered) == 0)
})

test_that("error if too many levels in factor for partitioning", {
    n <- 100
    dt <- data.frame(x=factor(1:100, ordered=F), y=rbinom(n, 1, 0.5))

    expect_error(
        train(data=dt, response_name="y", unordered_predictors="partition"),
        "Too many levels in unordered categorical variable x",
        fixed=T
    )
})

test_that("same out-of-bag error for character or unordered factor with partitioning", {
    set.seed(2)
    rf.class.char <- train(data=dt.class, response_name="y",
                           unordered_predictors="partition")
    dt.class.factor <- dt.class
    dt.class.factor$x <- factor(dt.class.factor$x, ordered=FALSE)
    set.seed(2)
    rf.class.factor <- train(data=dt.class.factor, response_name="y",
                            unordered_predictors="partition")

    expect_equal(rf.class.char$oob_error, rf.class.factor$oob_error)

    set.seed(2)
    rf.num.char <- train(data=dt.num, response_name="y",
                         unordered_predictors="partition")
    dt.num.factor <- dt.num
    dt.num.factor$x <- factor(dt.num.factor$x, ordered=FALSE)
    set.seed(2)
    rf.num.factor <- train(data=dt.num.factor, response_name="y",
                           unordered_predictors="partition")

    expect_equal(rf.num.char$oob_error, rf.num.factor$oob_error)
})

test_that("same out-of-bag error for character or unordered factor re-ordered by PCA score", {
    set.seed(2)
    rf.class.char <- train(data=dt.class, response_name="y",
                           unordered_predictors="order")
    dt.class.factor <- dt.class
    dt.class.factor$x <- factor(dt.class.factor$x, ordered=FALSE)
    set.seed(2)
    rf.class.factor <- train(data=dt.class.factor, response_name="y",
                            unordered_predictors="order")

    expect_equal(rf.class.char$oob_error, rf.class.factor$oob_error)

    set.seed(2)
    rf.num.char <- train(data=dt.num, response_name="y",
                         unordered_predictors="order")
    dt.num.factor <- dt.num
    dt.num.factor$x <- factor(dt.num.factor$x, ordered=FALSE)
    set.seed(2)
    rf.num.factor <- train(data=dt.num.factor, response_name="y",
                           unordered_predictors="order")

    expect_equal(rf.num.char$oob_error, rf.num.factor$oob_error)
})

test_that("can train classification forest when single-level predictor is re-ordered by PCA score", {
    n <- 20
    dt.bin <- data.frame(x=sample(c("A"), n, replace=T),
                         y=factor(sample(c("A", "B"), n, replace=T)),
                         stringsAsFactors=F)
    expect_silent(
        train(data=dt.bin, response_name="y", unordered_predictors="order")
    )
    dt.factor <- data.frame(x=sample(c("A"), n, replace=T),
                            y=factor(sample(c("A", "B", "C", "D"), n, replace=T)),
                            stringsAsFactors=F)
    expect_silent(
        train(data=dt.factor, response_name="y", unordered_predictors="order")
    )
})

test_that("can train regression forest when single-level predictor is re-ordered by PCA score", {
    n <- 20
    dt.bin <- data.frame(x=sample(c("A"), n, replace=T), y=rnorm(n),
                         stringsAsFactors=F)
    expect_silent(
        train(data=dt.bin, response_name="y", unordered_predictors="order")
    )
    dt.factor <- data.frame(x=sample(c("A"), n, replace=T), y=rnorm(n),
                            stringsAsFactors=F)
    expect_silent(
        train(data=dt.factor, response_name="y", unordered_predictors="order")
    )
})

test_that("same result when training forests if no unordered factors", {
    set.seed(100)
    rf1.class <- train(data=iris, response_name="Species",
                       unordered_predictors="ignore")
    set.seed(100)
    rf2.class <- train(data=iris, response_name="Species",
                       unordered_predictors="order")
    set.seed(100)
    rf3.class <- train(data=iris, response_name="Species",
                       unordered_predictors="partition")

    expect_equal(rf1.class$oob_error, rf2.class$oob_error)
    expect_equal(rf1.class$oob_error, rf3.class$oob_error)
    #expect_equal(rf1$confusion.matrix, rf2$confusion.matrix)
    #expect_equal(rf1$confusion.matrix, rf3$confusion.matrix)
    is.factor <- sapply(iris, is.factor)
    set.seed(100)
    rf1.num <- train(data=iris[!is.factor], response_name="Sepal.Length",
                     unordered_predictors="ignore")
    set.seed(100)
    rf2.num <- train(data=iris[!is.factor], response_name="Sepal.Length",
                     unordered_predictors="order")
    set.seed(100)
    rf3.num <- train(data=iris[!is.factor], response_name="Sepal.Length",
                     unordered_predictors="partition")

    expect_equal(rf1.num$oob_error, rf2.num$oob_error)
    expect_equal(rf1.num$oob_error, rf3.num$oob_error)
})

test_that("can train classification forest with unordered predictors and 'extratrees'", {
    expect_silent(
        rf <- train(data=dt.class, response_name="y", split_rule="extratrees",
                    unordered_predictors="partition")
    )
    expect_is(rf, "literanger")
})

test_that("can train regression forest when single-level unordered predictors are ordered by PCA score", {
  # Regression
    dt.num <- data.frame(x=sample(c("A"), n, replace=T),
                         y=rnorm(n),
                         stringsAsFactors=F)
    expect_silent(
        train(data=dt.num, response_name="y", unordered_predictors="order")
    )
})

test_that("can train regression forest with unordered predictors and 'extratrees'", {
    expect_silent(
        rf <- train(data=iris, response_name="Sepal.Length",
                    split_rule="extratrees",
                    unordered_predictors="partition")
    )
    expect_is(rf, "literanger")
})

test_that("trained regression forest with unordered predictors and 'extratrees' has acceptable out-of-bag error", {
    rf <- train(data=iris, response_name="Sepal.Length",
                split_rule="extratrees", unordered_predictors="partition")
    expect_gt(1 - rf$oob_error / var(iris$Sepal.Length), 0.5)
})

test_that("maximally selected rank statistics metric fails with partitioning", {
    expect_error(
        train(data=iris, response_name="Sepal.Length", split_rule="maxstat",
              unordered_predictors="partition"),
        paste("Unordered factor splitting not implemented for 'maxstat' or",
              "'beta' splitting rule."),
        fixed=T
    )
})

test_that("warning for maximally selected rank statistics with PCA ordering", {
    expect_warning(
        train(data=iris, response_name="Sepal.Length", split_rule="maxstat",
              unordered_predictors="order"),
        paste("The 'order' mode for unordered factor handling with the",
              "'maxstat' splitrule is experimental."),
        fixed=T
    )
})

test_that("can predict unobserved levels given unordered predictors", {
    set.seed(1)
    n <- 20
    dt.train <- data.frame(x1=sample(c("A", "B", "C"), n, replace=T),
                           x2=sample(c("A", "B", "C"), n, replace=T),
                           y=rbinom(n, 1, 0.5),
                           stringsAsFactors=F)

    dt.test <- data.frame(x1=sample(c("A", "B", "C", "D"), n, replace=T),
                          x2=sample(c("A", "B", "C", "D"), n, replace=T),
                          stringsAsFactors=F)

    rf.ignore <- train(data=dt.train, response_name="y",
                       unordered_predictors="ignore")
    expect_warning(
        predict(rf.ignore, newdata=dt.test),
        "Predictor levels found that were not present during training", fixed=T
    )

    rf.partition <- train(data=dt.train, response_name="y",
                          unordered_predictors="partition")
    expect_warning(
        predict(rf.partition, newdata=dt.test),
        "Predictor levels found that were not present during training", fixed=T
    )

    rf.order <- train(data=dt.train, response_name="y",
                      unordered_predictors="order")
    expect_warning(
        predict(rf.order, newdata=dt.test),
        "Predictor levels found that were not present during training", fixed=T
    )
})

test_that("can train on NA factor levels when using re-ordering via PCA score", {
    dt <- data.frame(x=addNA(factor(c("a", "a", NA, NA, "b", "b"))),
                     y=c(1, 2, 3, 4, 5, 6))
    expect_silent(
        train(data=dt, response_name="y", unordered_predictors="order")
    )
})

test_that("can use re-ordering via PCA score when numerics in data", {
    n <- 20

    dt.bin <- data.frame(x1=sample(c("A", "B", "C"), n, replace=T),
                         x2=sample(1:3, n, replace=T),
                         y=factor(sample(c("A", "B"), n, replace=T)),
                         stringsAsFactors=F)
    expect_silent(
        rf.bin <- train(data=dt.bin, response_name="y",
                        unordered_predictors="order")
    )
    expect_silent(predict(rf.bin, newdata=dt.bin))

  # Multiclass classification
    dt.fac <- dt.bin
    dt.fac$y <- factor(sample(c("A", "B", "C", "D"), n, replace=T))
    expect_silent(
        rf.fac <- train(data=dt.fac, response_name="y",
                        unordered_predictors="order")
    )
    expect_silent(predict(rf.fac, newdata=dt.fac))

  # Regression
    dt.num <- dt.class
    dt.num$y <- rnorm(n)
    expect_silent(
        rf.num <- train(data=dt.num, response_name="y",
                        unordered_predictors="order")
    )
    expect_silent(predict(rf.num, newdata=dt.num))
})

test_that("can use partitioning with a large number of levels", {
    n <- 43
    dt <- data.frame(x=factor(1:n, ordered = FALSE),  y=rbinom(n, 1, 0.5))

    expect_silent(
        rf <- train(data=dt, response_name="y", split_rule = "extratrees")
    )
  #  max_split <- max(sapply(1:rf$n_tree, function(i) {
  #      max(log2(rf$forest$split.values[[i]]), na.rm=T)
  #  }))
  #  expect_lte(max_split, n)
})

