literanger: A fast implementation of random forests for multiple imputation
===========================================================================

_Stephen Wade_

🚧 Under construction 🚧


See [`ranger`][ranger_cran] (Wright and Ziegler, 2017).

MICE: See [`mice`][mice_cran] (Doove et al, 2014).



[mice_cran]: https://cran.r-project.org/package=MICE
[ranger_cran]: https://cran.r-project.org/package=ranger


## Example

```r
require(literanger)
# TODO: 
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


### Alternatives

For those looking for a predictive mean matching approach to random forests in
multiple imputation:

1.  [`miceRanger`][miceranger_cran].
2.  [`missRanger`][missranger_cran].

[miceranger_cran]: https://cran.r-project.org/package=miceRanger
[missranger_cran]: https://cran.r-project.org/package=missRanger


## To-do

Not exhaustive:

-   finish documentation, e.g. this README;
-   prepare CRAN submission;
-   implement variable importance measures;
-   probability and survival forests?


## References

Doove, L.L., Van Buuren, S. and Dusseldorp, E., 2014. Recursive partitioning for
missing data imputation in the presence of interaction effects. _Computational
Statistics & Data Analysis, 72_, pp. 92-104.
[doi.10.1016/j.csda.2013.10.025](https://dx.doi.org/10.1016/j.csda.2013.10.025)

Wright, M. N. and Ziegler, A., 2017. ranger: A fast implementation of random
forests for high dimensional data in C++ and R. _Journal of Statistical
Software, 77_(i01), pp. 1-17.
[doi.10.18637/jss.v077.i01](https://dx.doi.org/10.18637/jss.v077.i01).

