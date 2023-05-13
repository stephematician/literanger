# literanger submission

-   This is a new release (fourth re-submission).
-   Removed 'C++11' system requirement from previous submission.
-   Updated case in title and replace T/F with TRUE/FALSE.
-   Passes r-lib/actions/examples/check-standard.yaml
-   Below are R CMD CHECK results from win-builder and my machine (Ubuntu 22.04).


# `win-builder` R CMD CHECK results 

```
❯ checking CRAN incoming feasibility ... [10s] NOTE
  Maintainer: 'Stephen Wade <stephematician@gmail.com>'
  
  New submission
  
  Possibly misspelled words in DESCRIPTION:
    Buuren (15:62)
    Doove (20:17)
    al (12:75, 20:26)
    et (12:72, 20:23)

0 errors ✔ | 0 warnings ✔ | 1 note ✖
```

-   The misspelled words above are names.


# Ubuntu 22.04 (my machine) R CMD CHECK results

```
❯ checking installed package size ... NOTE
    installed size is  8.0Mb
    sub-directories of 1Mb or more:
      libs   7.9Mb

0 errors ✔ | 0 warnings ✔ | 1 note ✖
```

-   The NOTE is due to linking to the Eigen3 library, which I use to
    support sparse matrices.

