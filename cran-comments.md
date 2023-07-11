# literanger submission

Submit version 0.0.2: bug fix and performance improvements.

-   Passes r-lib/actions/examples/check-standard.yaml.
-   Fixes container overflow identified in:
    https://www.stats.ox.ac.uk/pub/bdr/memtests/clang-ASAN/literanger/00check.log


# `win-builder` R CMD CHECK results 

```
  Maintainer: 'Stephen Wade <stephematician@gmail.com>'
  
  New submission
  
  Package was archived on CRAN
  
  Possibly misspelled words in DESCRIPTION:
    Buuren (15:62)
    Doove (20:17)
    al (12:75, 20:26)
    et (12:72, 20:23)
  
  CRAN repository db overrides:
    X-CRAN-Comment: Archived on 2023-07-10 as issues were not corrected
      in time.

0 errors ✔ | 0 warnings ✔ | 1 note ✖
```

-   The misspelled words above are names.


# Ubuntu 22.04 (my machine) R CMD CHECK results

```
── R CMD check results ─────────────────────────────────── literanger 0.0.2 ────
Duration: 52.8s

❯ checking installed package size ... NOTE
    installed size is  7.7Mb
    sub-directories of 1Mb or more:
      libs   7.6Mb

0 errors ✔ | 0 warnings ✔ | 1 note ✖
```

-   The NOTE is likely due to linking cpp11 and other C++11 standard libraries.

