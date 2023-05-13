
set.seed(42)
rf_class_df <- train(data=iris, response_name="Species")

test_that("can use case weights", {
    expect_silent(
        train(data=iris, response_name="Species",
              case_weights = rep(1, nrow(iris)))
    )
  ## Should only predict setosa now
    rf <- train(data=iris, response_name="Species",
                case_weights=c(rep(1, 50), rep(0, 100)))
    pred <- predict(rf, newdata=iris)$values
    expect_true(all(pred == "setosa"))
})

test_that("result is independent of column ordering", {
    dat <- iris[, c(1:2, 5, 3:4)]
    rf <- train(data=dat, response_name="Species")
    set.seed(42)
    pred_iris <- predict(rf, newdata=iris)$values
    set.seed(42)
    pred_shufl <- predict(rf, newdata=dat)$values
    expect_equal(pred_iris, pred_shufl)
    expect_gte(mean(pred_shufl == dat$Species), 0.9)
  ## No response column
    expect_gte(mean(predict(rf, newdata=dat[, -3])$values == dat$Species), 0.9)
})

test_that("meaningful predictions with max_depth = 1", {
    rf_class <- train(data=iris, response_name="Species", max_depth=1)
    pred_class <- predict(rf_class, newdata=iris)$values
    expect_true(all(pred_class %in% levels(iris$Species)))

    rf_num <- train(data=iris, response_name="Sepal.Length", max_depth=1)
    pred_num <- predict(rf_num, newdata=iris)$values
    expect_gte(min(pred_num), min(iris$Sepal.Length))
    expect_lte(max(pred_num), max(iris$Sepal.Length))
})

test_that("error if predictors have missing values", {
    dat <- iris
    dat[4, 4] <- dat[25,1] <- NA
    expect_error(
        predict(rf_class_df, newdata=dat),
        "Missing values in the predictors: Sepal.Length, Petal.Width."
    )
})

test_that("error if unknown value for prediction type", {
    expect_error(predict(rf_class_df, newdata=iris, prediction_type="foo"))
})

test_that("missing value allowed in discarded column", {
    rf <- train(data=iris, response_name="Species",
                predictor_names="Sepal.Length")
    dat <- iris
    dat[1, "Sepal.Width"] <- NA
    expect_silent(predict(rf, newdata=dat))
})

test_that("can omit response variable during prediction", {
    n <- 50
    dat_ftor <- data.frame(x1=runif(n), x2=runif(n),
                           y=factor(rbinom(n, 1, 0.5)))
    rf_ftor <- train(data=dat_ftor, response_name="y")
    expect_silent(predict(rf_ftor, newdata=dat_ftor))
    expect_silent(predict(rf_ftor, newdata=dat_ftor[, 1:2]))
})

test_that("can predict when a predictor is named 'forest'", {
    dat <- iris
    dat$forest <- rnorm(150)
    rf <- train(data=dat, response_name="Species")
    expect_silent(predict(rf, newdata=dat))
})

test_that("can predict when a predictor is named 'none'", {
    dat <- data.frame(y=rbinom(100, 1, .5),
                      x=rbinom(100, 1, .5),
                      none=rbinom(100, 1, .5))
    rf <- train(data=dat, response_name="y")
    expect_silent(predict(rf, newdata=dat))
})

if (requireNamespace("tibble", quietly=TRUE)) tb <- tibble::as_tibble(iris)

test_that("same result for prediction with tibble vs data.frame", {
    skip_if_not_installed("tibble")
    set.seed(42)
    rf_tb <- train(data=tb, response_name="Species")
    set.seed(57)
    pred_tb_tb <- predict(rf_tb, newdata=tb)
    set.seed(57)
    pred_tb_df <- predict(rf_tb, newdata=iris)

    set.seed(42)
    rf_df <- train(data=iris, response_name="Species")
    set.seed(57)
    pred_df_tb <- predict(rf_df, newdata=tb)
    set.seed(57)
    pred_df_df <- predict(rf_df, newdata=iris)

    expect_equal(rf_df$oob_error, rf_tb$oob_error)
    expect_equal(pred_df_df, pred_df_tb)
    expect_equal(pred_df_tb, pred_tb_tb)
    expect_equal(pred_tb_tb, pred_tb_df)
})

test_that("result of prediction has correct class", {
    pred <- predict(rf_class_df, newdata=iris)
    expect_is(pred, "literanger_prediction")
})

test_that("can get 'inbag' prediction", {
    expect_silent(predict(rf_class_df, newdata=iris, prediction_type="inbag"))
})

test_that("can get 'nodes' prediction", {
    expect_silent(
        pred_nodes <- predict(rf_class_df, newdata=iris,
                              prediction_type="nodes")
    )
    expect_is(pred_nodes$nodes, "matrix")
    expect_equal(storage.mode(pred_nodes$nodes), "integer")
})

test_that("default prediction type is 'bagged'", {
    set.seed(42)
    pred_df <- predict(rf_class_df, newdata=iris)
    set.seed(42)
    pred_bag_df <- predict(rf_class_df, newdata=iris, prediction_type="bagged")
    expect_equal(pred_df, pred_bag_df)
})

