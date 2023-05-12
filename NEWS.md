# literanger 0.0.0.9000

To be released as 0.0.1

This is the initial release of literanger, a refactoring and adaptation of the
ranger package <https://github.com/imbs-hl/ranger> for random forests. The
purpose of this update was to refactor the prediction code to enable efficient
prediction when embedded into the multiple imputation algorithm proposed by
Doove et al in:

Doove, L. L., Van Buuren, S., & Dusseldorp, E. (2014). Recursive partitioning
for missing data imputation in the presence of interaction effects.
_Computational statistics & data analysis_, 72, 92-104.

Currently supports:

-   Classification and regression trees/forests
-   Bagged prediction or prediction given by drawing from a random tree and a
    random in-bag value.





