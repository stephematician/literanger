# classification...
test_that("predict.all for classification returns numeric matrix of size trees x n", {
  rf <- ranger(Species ~ ., iris, num.trees = 5, write.forest = TRUE)
  pred <- predict(rf, iris, predict.all = TRUE)
  expect_is(pred$predictions, "matrix")
  expect_equal(dim(pred$predictions), 
              c(nrow(iris), rf$num.trees))
})
# TODO: need to adjust this to "check that it is actally bagging"
test_that("Majority vote of predict.all for classification is equal to forest prediction", {
  rf <- ranger(Species ~ ., iris, num.trees = 5, write.forest = TRUE)
  pred_forest <- predict(rf, iris, predict.all = FALSE)
  pred_trees <- predict(rf, iris, predict.all = TRUE)
  ## Majority vote, NA for ties
  pred_num <- apply(pred_trees$predictions, 1, function(x) {
    res <- which(tabulate(x) == max(tabulate(x)))
    if (length(res) == 1) {
      res
    } else {
      NA
    }
  })
  pred <- integer.to.factor(pred_num, rf$forest$levels)
  idx <- !is.na(pred)
  expect_equal(pred[idx], pred_forest$predictions[idx])
})
# TODO: confusion matrix tests...
test_that("confusion matrix is of right dimension", {
  expect_equal(dim(rg.class$confusion.matrix), 
               rep(nlevels(iris$Species), 2))
})

test_that("confusion matrix has right dimnames", {
  expect_equal(dimnames(rg.class$confusion.matrix),
               list(true = levels(iris$Species), predicted = levels(iris$Species)))
})

test_that("confusion matrix rows are the true classes", {
  expect_equal(as.numeric(rowSums(rg.class$confusion.matrix)), 
               as.numeric(table(iris$Species)))
})

test_that("confusion matrix rows are the true classes if using case weights", {
  rf <- ranger(Species ~ ., data = iris, num.trees = 5, 
               case.weights = c(rep(100, 5), rep(5, 145)))
  expect_equal(as.numeric(rowSums(rf$confusion.matrix)), 
               as.numeric(table(iris$Species)))
})

# TODO: inbag stuff?

test_that("Number of samples is right sample fraction, replace=FALSE, default", {
  rf <- ranger(Species ~ ., iris, num.trees = 5, keep.inbag = TRUE, replace = FALSE)
  num.inbag <- sapply(rf$inbag.counts, function(x) {
    sum(x > 0)
  })
  sample.fraction <- mean(num.inbag/nrow(iris))
  
  expect_gt(sample.fraction, 0.6)
  expect_lt(sample.fraction, 0.7)
})

test_that("Number of samples is right sample fraction, replace=FALSE, 0.3", {
  rf <- ranger(Species ~ ., iris, num.trees = 5, keep.inbag = TRUE, replace = FALSE, sample.fraction = 0.3)
  num.inbag <- sapply(rf$inbag.counts, function(x) {
    sum(x > 0)
  })
  sample.fraction <- mean(num.inbag/nrow(iris))
  
  expect_gt(sample.fraction, 0.25)
  expect_lt(sample.fraction, 0.35)
})
test_that("Number of samples is right sample fraction, replace=TRUE, default", {
  rf <- ranger(Species ~ ., iris, num.trees = 5, keep.inbag = TRUE, replace = TRUE)
  num.inbag <- sapply(rf$inbag.counts, function(x) {
    sum(x > 0)
  })
  
  sample.fraction <- mean(num.inbag/nrow(iris))
  expected.sample.fraction <- 1-exp(-1)
  
  expect_gt(sample.fraction, expected.sample.fraction-0.05)
  expect_lt(sample.fraction, expected.sample.fraction+0.05)
})

test_that("Number of samples is right sample fraction, replace=TRUE, 0.5", {
  rf <- ranger(Species ~ ., iris, num.trees = 5, keep.inbag = TRUE, replace = TRUE, sample.fraction = 0.5)
  num.inbag <- sapply(rf$inbag.counts, function(x) {
    sum(x > 0)
  })
  
  sample.fraction <- mean(num.inbag/nrow(iris))
  expected.sample.fraction <- 1-exp(-0.5)
  
  expect_gt(sample.fraction, expected.sample.fraction-0.05)
  expect_lt(sample.fraction, expected.sample.fraction+0.05)
})

test_that("Number of samples is right sample fraction, replace=FALSE, 0.3, weighted", {
  rf <- ranger(Species ~ ., iris, num.trees = 5, keep.inbag = TRUE, replace = FALSE, sample.fraction = 0.3, case.weights = runif(nrow(iris)))
  num.inbag <- sapply(rf$inbag.counts, function(x) {
    sum(x > 0)
  })
  sample.fraction <- mean(num.inbag/nrow(iris))
  
  expect_gt(sample.fraction, 0.25)
  expect_lt(sample.fraction, 0.35)
})

test_that("Number of samples is right sample fraction, replace=TRUE, 0.5, weighted", {
  rf <- ranger(Species ~ ., iris, num.trees = 5, keep.inbag = TRUE, replace = TRUE, sample.fraction = 0.5, case.weights = runif(nrow(iris)))
  num.inbag <- sapply(rf$inbag.counts, function(x) {
    sum(x > 0)
  })
  
  sample.fraction <- mean(num.inbag/nrow(iris))
  expected.sample.fraction <- 1-exp(-0.5)
  
  expect_gt(sample.fraction, expected.sample.fraction-0.05)
  expect_lt(sample.fraction, expected.sample.fraction+0.05)
})


# TODO: classification tests (inbag?)
test_that("Inbag counts match sample fraction, classification", {
  ## With replacement
  rf <- ranger(Species ~ ., iris, num.trees = 5, sample.fraction = c(0.2, 0.3, 0.4), 
               replace = TRUE, keep.inbag = TRUE)
  inbag <- do.call(cbind, rf$inbag.counts)
  expect_equal(unique(colSums(inbag[iris$Species == "setosa", ])), 30)
  expect_equal(unique(colSums(inbag[iris$Species == "versicolor", ])), 45)
  expect_equal(unique(colSums(inbag[iris$Species == "virginica", ])), 60)
  
  ## Without replacement
  rf <- ranger(Species ~ ., iris, num.trees = 5, sample.fraction = c(0.1, 0.2, 0.3), 
               replace = FALSE, keep.inbag = TRUE)
  inbag <- do.call(cbind, rf$inbag.counts)
  expect_equal(unique(colSums(inbag[iris$Species == "setosa", ])), 15)
  expect_equal(unique(colSums(inbag[iris$Species == "versicolor", ])), 30)
  expect_equal(unique(colSums(inbag[iris$Species == "virginica", ])), 45)
  
  ## Different order, without replacement
  dat <- iris[c(51:100, 101:150, 1:50), ]
  rf <- ranger(Species ~ ., dat, num.trees = 5, sample.fraction = c(0.1, 0.2, 0.3), 
               replace = FALSE, keep.inbag = TRUE)
  inbag <- do.call(cbind, rf$inbag.counts)
  expect_equal(unique(colSums(inbag[dat$Species == "setosa", ])), 15)
  expect_equal(unique(colSums(inbag[dat$Species == "versicolor", ])), 30)
  expect_equal(unique(colSums(inbag[dat$Species == "virginica", ])), 45)
})

test_that("OOB error is correct for 1 tree, classification", {
  n <- 50
  dat <- data.frame(y = factor(rbinom(n, 1, .5)), x = rnorm(n))
  rf <- ranger(y ~ ., dat, num.trees = 1)
  expect_equal(rf$prediction.error, mean(rf$predictions != dat$y, na.rm = TRUE))
})

test_that("OOB error is correct for 1 tree, regression", {
  n <- 50
  dat <- data.frame(y = rbinom(n, 1, .5), x = rnorm(n))
  rf <- ranger(y ~ ., dat, num.trees = 1)
  expect_equal(rf$prediction.error, mean((dat$y - rf$predictions)^2, na.rm = TRUE))
})

test_that("Split points are at (A+B)/2 for numeric features, regression variance splitting", {
  dat <- data.frame(y = rbinom(100, 1, .5), x = rbinom(100, 1, .5))
  rf <- ranger(y ~ x, dat, num.trees = 10)
  split_points <- sapply(1:rf$num.trees, function(i) {
    res <- treeInfo(rf, i)$splitval
    res[!is.na(res)]
  })
  expect_equal(split_points, rep(0.5, rf$num.trees))
})

test_that("Split points are at (A+B)/2 for numeric features, regression maxstat splitting", {
  dat <- data.frame(y = rbinom(100, 1, .5), x = rbinom(100, 1, .5))
  rf <- ranger(y ~ x, dat, num.trees = 10, splitrule = "maxstat", alpha = 1)
  split_points <- sapply(1:rf$num.trees, function(i) {
    res <- treeInfo(rf, i)$splitval
    res[!is.na(res)]
  })
  expect_equal(split_points, rep(0.5, rf$num.trees))
})

test_that("Split points are at (A+B)/2 for numeric features, classification", {
  dat <- data.frame(y = factor(rbinom(100, 1, .5)), x = rbinom(100, 1, .5))
  rf <- ranger(y ~ x, dat, num.trees = 10)
  split_points <- sapply(1:rf$num.trees, function(i) {
    res <- treeInfo(rf, i)$splitval
    res[!is.na(res)]
  })
  expect_equal(split_points, rep(0.5, rf$num.trees))
})

test_that("Tree depth creates trees of correct size", {
  # Recursive function to get tree depth
  depth <- function(rf, tree, i) {
    left <- rf$forest$child.nodeIDs[[tree]][[1]][i] + 1
    right <- rf$forest$child.nodeIDs[[tree]][[2]][i] + 1
    if (left <= 1) {
      0
    } else {
      1 + max(c(depth(rf, tree, left), depth(rf, tree, right)))
    }
  }
  forest_depth <- function(rf) {
    sapply(1:rf$num.trees, depth, rf = rf, i = 1)
  }
  
  # Depth 1
  rf <- ranger(Species ~ ., iris, num.trees = 5, max.depth = 1)
  expect_true(all(forest_depth(rf) <= 1))
  
  # Depth 4
  rf <- ranger(Species ~ ., iris, num.trees = 5, max.depth = 4)
  expect_true(all(forest_depth(rf) <= 4))
  
  # Random depth (deeper trees)
  max.depth <- round(runif(1, 1, 20))
  dat <- data.frame(y = runif(100, 0, 1), x = runif(100, 0, 1))
  rf <- ranger(y ~ ., dat, num.trees = 5, min.node.size = 1, max.depth = max.depth)
  expect_true(all(forest_depth(rf) <= max.depth))
})

test_that("Tree depth 0 equivalent to unlimited", {
  set.seed(200)
  rf1 <- ranger(Species ~ ., iris, num.trees = 5, max.depth = 0)
  
  set.seed(200)
  rf2 <- ranger(Species ~ ., iris, num.trees = 5)
  
  expect_equal(sapply(rf1$forest$split.varIDs, length), 
               sapply(rf2$forest$split.varIDs, length))
})

# TODO: prediction num tress
test_that("If num.trees set, these number is used for predictions", {
  rf <- ranger(Species ~ ., iris, num.trees = 5, write.forest = TRUE)
  pred <- predict(rf, iris, predict.all = TRUE, num.trees = 3)
  expect_equal(pred$num.trees, 3)
  expect_equal(dim(pred$predictions), c(nrow(iris), 3))
})

test_that("If num.trees not set, all trees are used for prediction", {
  rf <- ranger(Species ~ ., iris, num.trees = 5, write.forest = TRUE)
  pred <- predict(rf, iris, predict.all = TRUE)
  expect_equal(pred$num.trees, 5)
  expect_equal(dim(pred$predictions), c(nrow(iris), 5))
})

test_that("Warning if predicting with corrected impurity importance", {
  rf <- ranger(Species ~ ., iris, num.trees = 5, importance = "impurity_corrected")
  expect_warning(predict(rf, iris))
})











# TODO: variable importance (one day)

test_that("maxstat impurity importance is positive", {
  rf <- ranger(Surv(time, status) ~ ., veteran, num.trees = 5, 
               splitrule = "maxstat", importance = "impurity")
  expect_gt(mean(rf$variable.importance), 0)
  
  rf <- ranger(Sepal.Length ~ ., iris, num.trees = 5, 
               splitrule = "maxstat", importance = "impurity")
  expect_gt(mean(rf$variable.importance), 0)
})

test_that("maxstat corrected impurity importance is positive (on average)", {
  rf <- ranger(Surv(time, status) ~ ., veteran, num.trees = 50, 
               splitrule = "maxstat", importance = "impurity_corrected")
  expect_gt(mean(rf$variable.importance), 0)
  
  rf <- ranger(Sepal.Length ~ ., iris, num.trees = 5, 
               splitrule = "maxstat", importance = "impurity_corrected")
  expect_gt(mean(rf$variable.importance), 0)
})


test_that("Corrected importance working for sparse data", {
  rf <- ranger(data = dat_sparse, dependent.variable.name = "y", classification = TRUE, 
               num.trees = 5, importance = "impurity_corrected")
  expect_equal(names(rf$variable.importance), colnames(dat_sparse)[-1])
})





# TODO: some tests for regression output


test_that("Mean of predict.all for regression is equal to forest prediction", {
  rf <- ranger(Petal.Width ~ ., iris, num.trees = 5, write.forest = TRUE)
  pred_forest <- predict(rf, iris, predict.all = FALSE)
  pred_trees <- predict(rf, iris, predict.all = TRUE)
  expect_equal(rowMeans(pred_trees$predictions), pred_forest$predictions)
})








####### Split weights

test_that("draw predictor weights can be zero or one", {
    weights <- replicate(formals(ranger)$n_tree,
                         sample(c(0, 0, 1, 1)), simplify=F)
    rf <- ranger(data=iris, response_name="Species", 
                 draw_predictor_weights=weights)
    # selected_correctly <- sapply(1:rf$num.trees, function(i) {
    #     all(treeInfo(rf, i)[,"splitvarID"] %in% c(which(weights[[i]] > 0) - 1, NA))
    # })
    # expect_true(all(selected_correctly))
})

test_that("Tree-wise split select weights work with 0s", {
  num.trees <- 5
  weights <- replicate(num.trees, sample(c(0, 0, 0.5, 0.5)), simplify = FALSE)
  rf <- ranger(Species ~ ., iris, mtry = 2, num.trees = num.trees, 
               split.select.weights = weights)
  selected_correctly <- sapply(1:num.trees, function(i) {
    all(treeInfo(rf, i)[,"splitvarID"] %in% c(which(weights[[i]] > 0) - 1, NA))
  })
  expect_true(all(selected_correctly))
})