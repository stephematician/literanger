## Initialize the random forests
ind = 1:150 %in% sample(150, 100)

set.seed(2)
mod1 <- train(data=iris[ind, ], response_name="Species")
pred1 <- predict(mod1, newdata=iris[!ind, ])

set.seed(2)
mod2 <- train(data=iris[ind, ], response_name="Species")
pred2 <- predict(mod2, newdata=iris[!ind, ])

## Tests
test_that("same result with same seed", {
  expect_equal(pred1$predictions, pred2$predictions)
})

