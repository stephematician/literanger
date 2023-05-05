# -------------------------------------------------------------------------------
# This has been adapted from Ranger.
#
# Ranger is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Ranger is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Ranger. If not, see <http://www.gnu.org/licenses/>.
#
# Written by:
#
#   Marvin N. Wright
# Institut fuer Medizinische Biometrie und Statistik
# Universitaet zu Luebeck
# Ratzeburger Allee 160
# 23562 Luebeck
# Germany
#
# http://www.imbs-luebeck.de
# -------------------------------------------------------------------------------
#
# The adaptation was performed by Stephen Wade. The changes include eliminating
# features that were not required for the multiple imputation procedure.

# Order factor levels with PCA approach
# Reference: Coppersmith, D., Hong, S.J. & Hosking, J.R. (1999) Partitioning
# Nominal Attributes in Decision Trees. Data Min Knowl Discov 3:197.
# \doi{10.1023/A:1009869804967}.
#' @importFrom stats cov.wt prcomp
pca_order <- function(y, x) {
    x <- droplevels(x)
    if (nlevels(x) < 2) return(as.character(levels(x)))

  # Create contingency table of the nominal outcome with the nominal covariate
    N <- table(x, droplevels(y))

  # PCA of weighted covariance matrix of class probabilites
    P <- N/rowSums(N)
    S <- stats::cov.wt(P, wt = rowSums(N))$cov
    pc1 <- stats::prcomp(S, rank. = 1)$rotation
    score <- P %*% pc1

  # Return ordered factor levels
    as.character(levels(x)[order(score)])
}

