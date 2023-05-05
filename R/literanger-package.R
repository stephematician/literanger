#' \pkg{literanger}: Random forests for multiple imputation algorithms.
#'
#' 'literanger' is an adaption of the 'ranger' R package for training and
#' predicting from random forest models within multiple imputation problems.
#' 'ranger' is a fast implementation of random forests (Breiman, 2001) or
#' recursive partitioning, particularly suited for high dimensional data
#' (Wright et al, 2017a). literanger enables random forests to be embedded in
#' the fully conditional specification framework for multiple imputation known
#' as  'Multiple Imputation via Chained Equations' (Van Buuren 2007).
#' Specifically, literanger records terminal nodes' in-bag data as needed by
#' the algorithm described in Doove et al (2014). Alternatively, the usual
#' bagged prediction can be provided as in the imputation algorithm called
#' 'missForest' (Stekhoven et al, 2014).
#'
#' Classification and regression forests are implemented as in the original
#' Random Forest (Breiman, 2001) or using extremely randomized trees (Geurts et
#' al, 2006). "data.frame", "matrix", and sparse matrices ("dgCmatrix") are
#' supported.
#'
#' Split selection may be based on improvement in metrics such as the variance,
#' Gini impurity, beta log-likelihood (Weinhold et al, 2019), Hellinger distance
#' (Cieslak et al, 2012) or maximally selected statistics (Wright et al, 2017b).
#'
#' See <https://github.com/imbs-hl/ranger> for the development version of ranger
#' or <https://github.com/stephematician/literanger> for development version of
#' this package.
#'
#' For alternative approaches to multiple imputation that employ random forests,
#' see 'missRanger' <https://cran.r-project.org/package=missRanger> and
#' 'miceRanger' <https://cran.r-project.org/package=miceRanger>, which use
#' predictive mean matching combined with the original 'ranger' algorithm.
#'
#' This package was adapted from the "ranger" package for R Statistical
#' Software. The C++ core is provided under the same license terms as the
#' original C++ core in the 'ranger' package, namely the MIT license
#' <https://www.r-project.org/Licenses/MIT>. The wrappers in this package around
#' the core are licensed under the same terms of the corresponding components in
#' the 'ranger' R package, namely the GPL3 license
#' <https://www.r-project.org/Licenses/GPL-3>, <http:://www.gnu.org/licenses>.
#'
#' License statement for C++ core of ranger:
#'
#'    Copyright (c) [2014-2018] [Marvin N. Wright]
#'    
#'    This software may be modified and distributed under the terms of the MIT
#'    license.
#'    
#'    Please note that the C++ core of ranger is distributed under MIT license
#'    and the R package "ranger" under GPL3 license.
#'
#' License statement for the additions to the core in the 'ranger' R package:
#'
#'    Ranger is free software: you can redistribute it and/or modify it under
#'    the terms of the GNU General Public License as published by the Free
#'    Software Foundation, either version 3 of the License, or (at your option)
#'    any later version.
#'    
#'    Ranger is distributed in the hope that it will be useful, but WITHOUT ANY
#'    WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
#'    FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
#'    details.
#'    
#'    You should have received a copy of the GNU General Public License along
#'    with Ranger. If not, see <http://www.gnu.org/licenses/>.
#'    
#'    Written by:
#'    Marvin N. Wright
#'    Institut für Medizinische Biometrie und Statistik
#'    Universität zu Lübeck
#'    Ratzeburger Allee 160
#'    23562 Lübeck
#'    http://www.imbs-luebeck.d
#'
#' @references
#'
#' -   Cieslak, D. A., Hoens, T. R., Chawla, N. V., & Kegelmeyer, W. P. (2012).
#'     Hellinger distance decision trees are robust and skew-insensitive. _Data
#'     Mining and Knowledge Discovery_, 24, 136-158.
#'     \doi{10.1007/s10618-011-0222-1}.
#' -   Doove, L. L., Van Buuren, S., & Dusseldorp, E. (2014). Recursive
#'     partitioning for missing data imputation in the presence of interaction
#'     effects. _Computational Statistics & Data Analysis_, 72, 92-104.
#'     \doi{10.1016/j.csda.2013.10.025}.
#' -   Geurts, P., Ernst, D., & Wehenkel, L. (2006). Extremely randomized trees.
#'     _Machine Learning_, 63, 3-42. \doi{10.1007/s10994-006-6226-1}.
#' -   Stekhoven, D.J. and Buehlmann, P. (2012). MissForest--non-parametric
#'     missing value imputation for mixed-type data. _Bioinformatics_, 28(1),
#'     112-118. \doi{10.1093/bioinformatics/btr597}.
#' -   Van Buuren, S. (2007). Multiple imputation of discrete and continuous
#'     data by fully conditional specification. _Statistical Methods in Medical
#'     Research, 16(3), 219-242. \doi{10.1177/0962280206074463}.
#' -   Weinhold, L., Schmid, M., Wright, M. N., & Berger, M. (2019). A random
#'     forest approach for modeling bounded outcomes. _arXiv preprint_,
#'     arXiv:1901.06211. \doi{10.48550/arXiv.1901.06211}.
#' -   Wright, M. N., & Ziegler, A. (2017a). ranger: A Fast Implementation of
#'     Random Forests for High Dimensional Data in C++ and R. _Journal of
#'     Statistical Software_, 77, 1-17. \doi{10.18637/jss.v077.i01}.
#' -   Wright, M. N., Dankowski, T., & Ziegler, A. (2017b). Unbiased split
#'     variable selection for random survival forests using maximally selected
#'     rank statistics. _Statistics in medicine_, 36(8), 1272-1284.
#'     \doi{10.1002/sim.7212}.
#'
#' @author Stephen Wade <stephematician@gmail.com>, Marvin N Wright (original
#' 'ranger' package)
#'
#' @useDynLib literanger, .registration = TRUE
#'
#' @docType package
#' @name literanger-package
#' @md
NULL

# A recommended practise; unload this package's dynamic libraries,
# see http://r-pkgs.had.co.nz/src.html
.onUnload <- function (lib_path)
    library.dynam.unload("literanger", lib_path)

# This dummy function definition is included with the package to ensure that
# 'tools::package_native_routine_registration_skeleton()' generates the required
# registration info for the 'run_testthat_tests' symbol.
(function() {
  .Call("run_testthat_tests", FALSE, PACKAGE="literanger")
})
