literanger: A fast implementation of random forests for multiple imputation
===========================================================================

![R-CMD-Check](https://github.com/stephematician/literanger/actions/workflows/check-standard.yaml/badge.svg)

_Stephen Wade_

🚧 Under construction 🚧

`literanger` is an adaption of the [`ranger`][ranger_cran] R package for
training and predicting from random forest models within multiple imputation
algorithms. `ranger` is a fast implementation of random forests
([Breiman, 2001][brieman2001_doi]) or recursive partitioning, particularly
suited for high dimensional data ([Wright et al, 2017][wright2017_doi]).
`literanger` enables random forests to be embedded in the fully conditional
specification framework for multiple imputation known as 'Multiple Imputation
via Chained Equations' (Van Buuren 2007).

Implementations of multiple imputation with random forests include:

1.  [`mice`][mice_cran] which uses random forests to predict in a similar fashion
    to [Doove et al, 2014][doove2014_doi], i.e. for each observation, a draw is
    taken from the sample of all values that belong to the the terminal node of
    a randomly drawn tree.
2.  [`miceRanger`][miceranger_cran] and [`missRanger`][missranger_cran] which
    use predictive mean matching.

This package enables a minor variation on `mice`'s use of random forests.
The prediction can be drawn from the in-bag samples in the terminal node for
each _missing_ data point. Thus, the computational effort during prediction then
scales with the number of missing values, rather than with the product of the
size of the whole dataset and the number of trees (as in `mice`).

A more general advantage of this package is re-cycling of the trained forest
object and the separation of the (training) data from the forest, see `ranger` [issue #304](https://github.com/imbs-hl/ranger/issues/304).

A multiple imputation algorithm using this package is under development: called
[`smirf`][smirf_github].

[mice_cran]: https://cran.r-project.org/package=MICE
[miceranger_cran]: https://cran.r-project.org/package=miceRanger
[missranger_cran]: https://cran.r-project.org/package=missRanger
[ranger_cran]: https://cran.r-project.org/package=ranger
[smirf_github]: https://github.com/stephematician/smirf


## Example

```r
require(literanger)

train.idx <- sample(nrow(iris), 2/3 * nrow(iris))
iris.train <- iris[ train.idx, ]
iris.test  <- iris[-train.idx, ]
rg.iris <- train(data=iris.train, response_name="Species")
pred.iris.bagged <- predict(rg.iris, newdata=iris.test,
                            prediction_type="bagged")
pred.iris.inbag  <- predict(rg.iris, newdata=iris.test,
                            prediction_type="inbag")
# compare bagged vs actual test values
table(test.iris$Species, pred.iris.bagged$values)
# compare bagged prediction vs in-bag draw
table(pred.iris.bagged$values, pred.iris.inbag$values)
```


## Installation

Installation is easy using [`devtools`][devtools_cran]:

```r
library(devtools)
install_github('stephematician/literanger')
```

The [`cpp11`][cpp11_cran] package is also required, available on CRAN:

```r
install.packages('cpp11')
```

[cpp11_cran]: https://cran.r-project.org/package=cpp11
[devtools_cran]: https://cran.r-project.org/package=devtools


## To-do

Not exhaustive:

-   prediction type: terminal nodes for every tree (e.g. for mice algo);
-   ~~finish documentation, e.g. this README~~;
-   prepare CRAN submission;
-   implement variable importance measures;
-   probability and survival forests.


## References

Breiman, L. (2001). Random forests. _Machine learning_, 45, pp. 5-32.
[doi:10.1023/A:1010933404324](https://doi.org/10.1023/A:1010933404324).

Doove, L.L., Van Buuren, S. and Dusseldorp, E., 2014. Recursive partitioning for
missing data imputation in the presence of interaction effects. _Computational
Statistics & Data Analysis, 72_, pp. 92-104.
[doi:10.1016/j.csda.2013.10.025](https://doi.org/10.1016/j.csda.2013.10.025).

Wright, M. N. and Ziegler, A., 2017. ranger: A fast implementation of random
forests for high dimensional data in C++ and R. _Journal of Statistical
Software, 77_(i01), pp. 1-17.
[doi:10.18637/jss.v077.i01](https://doi.org/10.18637/jss.v077.i01).

[brieman2001_doi]: https://doi.org/10.1023/A:1010933404324
[doove2014_doi]: https://doi.org/10.1016/j.csda.2013.10.025
[wright2017_doi]: https://doi.org/10.18637/jss.v077.i01

