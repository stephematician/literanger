## Tests for character data

## Initialize random forests
dat <- iris
dat$Test <- paste0("AA", as.character(1:nrow(dat)))

## Tests
test_that("can train and predict if data has character vector", {
    expect_silent(rf <- train(data=dat, response_name="Species"))
    expect_silent(predict(rf, newdata=dat))
})

