# Changelog - literanger 

[![Common Changelog](https://common-changelog.org/badge.svg)](https://common-changelog.org)

## [0.0.3]() - 2023-xx-xx



## [0.0.2](https://github.com/stephematician/literanger/releases/tag/v0.0.1) - 2023-07-11

Update to pass CRAN's ASAN check

### Changed

-   Improve performance of node splitting ([`d3f6424`](https://github.com/stephematician/literanger/commit/d3f64245))

### Added

-   Add re-entrant log gamma to speed up beta splitting rule
    ([`d7f058d`](https://github.com/stephematician/literanger/commit/d7f058dd))
-   Minor fixes to documentation ([`91b6c6d`](https://github.com/stephematician/literanger/commit/91b6c6d),
    [`0f62d02`](https://github.com/stephematician/literanger/commit/0f62d027))

### Fixed

-   Fix potential illegal access and incorrect unweighted sampling without
    replacement ([`b6df5d9`](https://github.com/stephematician/literanger/commit/b6df5d9))


## [0.0.1](https://github.com/stephematician/literanger/releases/tag/v0.0.1) - 2023-06-25

_First release_

A refactoring and adaptation of the ranger package
<https://github.com/imbs-hl/ranger> for random forests. Has faster prediction
mode intended for embedding into the multiple imputation algorithm proposed by
Doove et al in:

Doove, L. L., Van Buuren, S., & Dusseldorp, E. (2014). Recursive partitioning
for missing data imputation in the presence of interaction effects.
_Computational statistics & data analysis_, 72, 92-104.

### Added

-   Fit classification and regression trees
-   Prediction via most frequent value or mean
-   Get predictions as terminal node identifiers in each tree or as a random
    draw from inbag values in a random tree

