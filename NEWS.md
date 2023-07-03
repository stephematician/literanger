literanger NEWS
===============

# version 0.0.2

Performance enhancements
-   Faster (correct) test for number of candidate values in node splitting.
-   Remove lock on log gamma (beta splitting rule).

Bug-fixes
-   Fix container overrun and incorrect (unweighted) sampling without
    replacement.


# version 0.0.1

This is the initial release of literanger, a refactoring and adaptation of the
ranger package <https://github.com/imbs-hl/ranger> for random forests. The
purpose of this update was to refactor the prediction code to enable efficient
prediction when embedded into the multiple imputation algorithm proposed by
Doove et al in:

Doove, L. L., Van Buuren, S., & Dusseldorp, E. (2014). Recursive partitioning
for missing data imputation in the presence of interaction effects.
_Computational statistics & data analysis_, 72, 92-104.

Currently supports:
-   Classification and regression trees/forests.
-   Prediction types:
    -   Conventional 'bagged' prediction (most frequent value or mean).
    -   Terminal node identifiers for all trees.
    -   Prediction given by drawing a tree for each prediction and then drawing
        an in-bag value from the terminal node.

