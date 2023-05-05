#include <testthat.h>

/* call declaration */
#include "cpp11_train.decl.h"

// Initialize a unit test context. This is similar to how you
// might begin an R test file with 'context()', expect the
// associated context should be wrapped in braced.
context("cpp11_train") {

  // The format for specifying tests is similar to that of
  // testthat's R functions. Use 'test_that()' to define a
  // unit test, and use 'expect_true()' and 'expect_false()'
  // to test the desired conditions.
  test_that("two plus two equals four") {
    //expect_true(twoPlusTwo() == 4);
    expect_true(1 == 1);
  }

}

