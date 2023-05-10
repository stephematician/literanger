# -------------------------------------------------------------------------------
# This file is part of 'literanger'. literanger was adapted from the 'ranger'
# package for R statistical software. ranger was authored by Marvin N Wright
# with the GNU General Public License version 3. The adaptation was performed by
# Stephen Wade in 2023. literanger carries the same license, terms, and
# permissions as ranger.
#
# literanger is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# literanger is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with literanger. If not, see <http://www.gnu.org/licenses/>.
#
# Written by:
#
#   Stephen Wade
#   Cancer Council New South Wales
#   Woolloomooloo NSW 2011
#   Australia
# -------------------------------------------------------------------------------

# Order factor levels with PCA approach
# Reference: Coppersmith, D., Hong, S.J. & Hosking, J.R. (1999) Partitioning
# Nominal Attributes in Decision Trees. Data Min Knowl Discov 3:197.
# \doi{10.1023/A:1009869804967}.
#' @author Marvin N Wright
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

