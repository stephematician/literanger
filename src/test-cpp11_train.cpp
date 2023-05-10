/*-------------------------------------------------------------------------------
 * This file is part of 'literanger'. literanger was adapted from the 'ranger'
 * package for R Statistical Software <https://www.r-project.org>. ranger was
 * authored by Marvin N Wright with the GNU General Public License version 3.
 * The adaptation was performed by Stephen Wade in 2023. literanger carries the
 * same license, terms, and
 * permissions as ranger.
 * 
 * literanger is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * literanger is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with literanger. If not, see <http://www.gnu.org/licenses/>.
 *
 * Written by:
 *
 *   Stephen Wade
 *   Cancer Council New South Wales
 *   Woolloomooloo NSW 2011
 *   Australia
 *-------------------------------------------------------------------------------
 */
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

