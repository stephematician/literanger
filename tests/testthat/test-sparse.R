if (requireNamespace("Matrix", quietly=T))
    iris_sparse <- Matrix::Matrix(data.matrix(iris), sparse=T)

## 0/1 sparse data
n <- 100
p <- 5
x <- replicate(p, rbinom(n, 1, .1))
y <- rbinom(n, 1, .5)
dat <- data.frame(y = y, x)
dat_matrix <- data.matrix(dat)
if (requireNamespace("Matrix", quietly=T))
    dat_sparse <- Matrix::Matrix(dat_matrix, sparse=T)

test_that("out-of-bag error same using sparse or non-sparse data for iris classification", {
    skip_if_not_installed("Matrix")
    set.seed(56)
    rf1 <- train(data=iris_sparse, response_name="Species", classification=T)
    set.seed(56)
    rf2 <- train(data=iris, response_name="Species")
  
    expect_equal(rf1$oob_error, rf2$oob_error)
})

test_that("prediction same using sparse data for iris classification", {
    skip_if_not_installed("Matrix")
    set.seed(56)
    rf1 <- train(data=iris_sparse, response_name="Species", classification=T)
    set.seed(56)
    rf2 <- train(data=iris, response_name="Species")
    pred1 <- factor(predict(rf1, newdata=iris_sparse)$values,
                    levels=rf2$response_values,
                    labels=rf2$response_levels)
    pred2 <- predict(rf2, newdata=iris)$values
    expect_equal(pred1, pred2)
})

test_that("out-of-bag error same with sparse data for binary classification", {
    skip_if_not_installed("Matrix")
    set.seed(56)
    rf1 <- train(data=dat_sparse, response_name="y", classification=T)
    set.seed(56)
    rf2 <- train(data=dat, response_name="y", classification=T)
  
    expect_equal(rf1$oob_error, rf2$oob_error)
})

test_that("prediction same with sparse data for binary classification", {
    skip_if_not_installed("Matrix")
    set.seed(56)
    rf1 <- train(data=dat_sparse, response_name="y", classification=T)
    set.seed(56)
    rf2 <- train(data=dat, response_name="y", classification=T)

    pred1 <- predict(rf1, newdata=dat_sparse, seed=123)$values
    pred2 <- predict(rf2, newdata=dat, seed=123)$values
    expect_equal(pred1, pred2)
})

test_that("out-of-bag error same with sparse data for 0/1 regression", {
    skip_if_not_installed("Matrix")
    set.seed(56)
    rf1 <- train(data=dat_sparse, response_name="y")
    set.seed(56)
    rf2 <- train(data=dat, response_name="y")
  
    expect_equal(rf1$oob_error, rf2$oob_error)
})

test_that("prediction same with sparse data for 0/1 regression", {
    skip_if_not_installed("Matrix")
    set.seed(56)
    rf1 <- train(data=dat_sparse, response_name="y")
    set.seed(56)
    rf2 <- train(data=dat, response_name="y")

    pred1 <- predict(rf1, newdata=dat_sparse, seed=123)$values
    pred2 <- predict(rf2, newdata=dat, seed=123)$values
    expect_equal(pred1, pred2)
})

test_that("prediction same if training or testing data is sparse", {
    skip_if_not_installed("Matrix")
    idx <- sample(nrow(iris), 2/3*nrow(iris))
    dat_train <- iris[idx, ]
    dat_test <- iris[-idx, ]
    train_sparse <- Matrix::Matrix(data.matrix(dat_train), sparse=T)
    test_sparse <- Matrix::Matrix(data.matrix(dat_test), sparse=T)
  # to convert from sparse-data response to original response
    unique_species <- unique(iris$Species)
    species_map <- setNames(unique_species, as.integer(unique_species))
  
    set.seed(42)
    rf1 <- train(data=dat_train, response_name="Species")
    pred1 <- predict(rf1, newdata=dat_test, seed=123)$values
    pred1_sparse <- predict(rf1, newdata=test_sparse, seed=123)$values

    set.seed(42)
    rf2 <- train(data=train_sparse, response_name="Species", classification=T)
    pred2 <- factor(
        unname(species_map[predict(rf2, newdata=dat_test, seed=123)$values]),
        levels=levels(iris$Species)
    )
    pred2_sparse <- factor(
        unname(species_map[predict(rf2, newdata=test_sparse, seed=123)$values]),
        levels=levels(iris$Species)
    )

    expect_equal(pred1, pred1_sparse)
    expect_equal(pred2, pred2_sparse)
    expect_equal(pred1, pred2)
})

