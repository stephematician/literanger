
test_that("get same prediction twice via set.seed for pseudo-rng seeding", {
    set.seed(42)
    ind <-1:150 %in% sample(150, 100)

    set.seed(2)
    rf1 <- train(data=iris[ind, ], response_name="Species")
    pred1 <- predict(rf1, newdata=iris[!ind, ])

    set.seed(2)
    rf2 <- train(data=iris[ind, ], response_name="Species")
    pred2 <- predict(rf2, newdata=iris[!ind, ])

    expect_equal(pred1$predictions, pred2$predictions)
})

test_that("get same prediction twice via 'seed' argument for pseudo-rng seeding", {
    set.seed(42)
    ind <-1:150 %in% sample(150, 100)

    rf1 <- train(data=iris[ind, ], response_name="Species", seed=1)
    pred1 <- predict(rf1, newdata=iris[!ind, ], seed=2)

    rf2 <- train(data=iris[ind, ], response_name="Species", seed=1)
    pred2 <- predict(rf2, newdata=iris[!ind, ], seed=2)

    expect_equal(pred1$predictions, pred2$predictions)
})


