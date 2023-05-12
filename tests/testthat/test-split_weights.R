
test_that("weights for drawing predictors (shared over forest) needs input of correct length", {
    expect_silent(
        train(data=iris, response_name="Species",
              draw_predictor_weights=c(0.1, 0.2, 0.3, 0.4))
    )
    expect_error(
        train(data=iris, response_name="Species",
              draw_predictor_weights=c(0.1, 0.2, 0.3)),
        paste("Size of 'draw_predictor_weights' (numeric) not equal to number",
              "of predictors."),
        fixed=T
    )
})

test_that("can have tree-wise weights for drawing predictors", {
    weights <- replicate(formals(literanger::train)$n_tree, runif(ncol(iris)-1),
                       simplify=F)
    expect_silent(train(data=iris, response_name="Species",
                        draw_predictor_weights =weights))

    extra_weights <- c(weights, list(runif(ncol(iris)-1)))
    expect_error(
        train(data=iris, response_name="Species",
              draw_predictor_weights=extra_weights),
        "Size of 'draw_predictor_weights' (list) not equal to number of trees.",
        fixed=T
    )
})

test_that("can provide names of predictors that are always candidates for split", {
    expect_silent(
        train(data=iris, response_name="Species", n_try=2,
              names_of_always_draw=c("Petal.Length", "Petal.Width"))
    )
})

test_that("can mix draw predictor weights and names that are always candidates", {
    iris_var <- setdiff(names(iris), 'Species')
    n_var <- length(iris_var)
    last_var <- iris_var[n_var]
    with_last_zero <- c(rep(1, n_var - 1), 0)
    expect_silent(
        train(data=iris, response_name="Species", n_try=n_var - 1,
              names_of_always_draw=last_var,
              draw_predictor_weights=with_last_zero)
    )
})

