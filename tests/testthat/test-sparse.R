if (requireNamespace("Matrix", quietly=TRUE))
    iris_sparse <- Matrix::Matrix(data.matrix(iris), sparse=TRUE)

test_that("out-of-bag error same using sparse or non-sparse data", {
    skip_if_not_installed("Matrix")
    set.seed(56)
    rf_class_sparse <- train(data=iris_sparse, response_name="Species",
                             classification=TRUE)
    set.seed(56)
    rf_class_df <- train(data=iris, response_name="Species")

    expect_equal(rf_class_df$oob_error, rf_class_sparse$oob_error)

    skip_if_not_installed("Matrix")
    set.seed(56)
    rf_num_sparse <- train(data=iris_sparse, response_name="Sepal.Length")
    set.seed(56)
    rf_num_df <- train(data=iris, response_name="Sepal.Length")

    expect_equal(rf_num_df$oob_error, rf_num_sparse$oob_error)

})

test_that("prediction same using sparse or non-sparse data", {
    skip_if_not_installed("Matrix")
    set.seed(56)
    rf_class_df <- train(data=iris, response_name="Species")
    pred_class_df <- predict(rf_class_df, newdata=iris)$values
    set.seed(56)
    rf_class_sparse <- train(data=iris_sparse, response_name="Species",
                             classification=TRUE)
    pred_class_sparse <- predict(rf_class_sparse, newdata=iris_sparse)$values
    pred_class_sparse <- factor(pred_class_sparse,
                                levels=rf_class_df$response_values,
                                labels=rf_class_df$response_levels)

    expect_equal(pred_class_df, pred_class_sparse)

    set.seed(56)
    rf_num_sparse <- train(data=iris_sparse, response_name="Sepal.Length")
    pred_num_sparse <- predict(rf_num_sparse, newdata=iris_sparse)$values
    set.seed(56)
    rf_num_df <- train(data=iris, response_name="Sepal.Length")
    pred_num_df <- predict(rf_num_df, newdata=iris)$values

    expect_equal(pred_num_df, pred_num_sparse)

})

test_that("prediction same if one of training or testing data is sparse", {
    skip_if_not_installed("Matrix")
    idx <- sample(nrow(iris), 2/3*nrow(iris))
    train_df <- iris[idx, ]
    test_df <- iris[-idx, ]
    train_sparse <- Matrix::Matrix(data.matrix(train_df), sparse=TRUE)
    test_sparse <- Matrix::Matrix(data.matrix(test_df), sparse=TRUE)
  # to convert from sparse-data response to original response
    unique_species <- unique(iris$Species)
    species_map <- setNames(unique_species, as.integer(unique_species))

    set.seed(42)
    rf_df <- train(data=train_df, response_name="Species")
    set.seed(57)
    pred_df_df <- predict(rf_df, newdata=test_df)$values
    set.seed(57)
    pred_df_sparse <- predict(rf_df, newdata=test_sparse)$values

    set.seed(42)
    rf_sparse <- train(data=train_sparse, response_name="Species",
                       classification=TRUE)
    set.seed(57)
    pred_sparse_df <- factor(
        unname(species_map[predict(rf_sparse, newdata=test_sparse)$values]),
        levels=levels(iris$Species)
    )
    set.seed(57)
    pred_sparse_sparse <- factor(
        unname(species_map[predict(rf_sparse, newdata=test_sparse)$values]),
        levels=levels(iris$Species)
    )

    expect_equal(pred_df_df, pred_df_sparse)
    expect_equal(pred_df_sparse, pred_sparse_df)
    expect_equal(pred_sparse_df, pred_sparse_sparse)
})

