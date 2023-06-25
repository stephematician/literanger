# literanger submission

This is a new release (fifth re-submission).

-   Removed Eigen dependency (eliminates header downloads).

Older changes compared to earlier submissions:

-   Removed 'C++11' system requirement from previous submission.
-   Updated case in title and replace T/F with TRUE/FALSE.
-   Passes r-lib/actions/examples/check-standard.yaml
-   Below are R CMD CHECK results from win-builder and my machine (Ubuntu 22.04).


# `win-builder` R CMD CHECK results 

```
❯ checking CRAN incoming feasibility ... [11s] NOTE
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
Duration: 48.3s

❯ checking installed package size ... NOTE
    installed size is  7.7Mb
    sub-directories of 1Mb or more:
      libs   7.5Mb

0 errors ✔ | 0 warnings ✔ | 1 note ✖
```

-   The NOTE is likely due to linking cpp11 and other C++11 standard libraries.

