\name{dimVAR}
\alias{dimVAR}
\title{
Assessing variances under factor solutions
}
\description{
\code{dimVAR} is a support function that assesses the proportion of and cumulative variances for a 
range of factor solutions.
}
\usage{
dimVAR(R, maxdim, graph = TRUE, verbose = TRUE)
}
\arguments{
  \item{R}{
  (Regularized) correlation \code{matrix}.
  }
  \item{maxdim}{
  A \code{numeric} integer or \code{integer} indicating the maximum factor dimension to be assessed.
  }
  \item{graph}{
  A \code{logical} indicating if the results should be visualized.
  }
  \item{verbose}{
  A \code{logical} indicating if the function should run silently.\cr
  Runs silently when \code{verbose = FALSE}.
  }
}
\details{
To assess a factor solution under \eqn{m} factors one might look at the proportion of explained variance.
The \code{dimVAR} function calculates the proportion of variance explained by any factor as well as 
the proportion of variance explained by all factors for each factor solution raging from 1 to \code{maxdim}.
Qualitatively, we want the proportion of variance explained by all factors to be appreciable (rules of thumb would say in excess of 70\%).
Moreover, one would want the proportion of variance explained by the \eqn{k}th factor in relation to the \eqn{(k - 1)}th factor to be appreciable and the proportion of variance of the \eqn{(k + 1)}th factor in relation to the \eqn{k}th factor to be negligible.

When \code{graph = TRUE} also a graph is returned visualizing the total cumulative variance against the dimension of the factor solution.
Hence, it plots the total cumulative variances against the respective factor solutions ranging from 1 to \code{maxdim}.
The point at which the graph flattens out is indicative of a formative number of latent factors.
}
\value{
Returns an object of class \code{list}.
\item{$CumVar}{Contains a \code{numeric} vector with the cumulative variances explained for each factor solution from 1 to \code{maxdim}.}
\item{$varianceTables}{This slot is itself a \code{list}. It contains, for each factor solution, a matrix with the sum of squares (SS), proportion variance (PV), and cumulative variance (CV) for each respective latent feature. Say one wants to access the variance table for a solution under 5-factors. Then one can call \code{$varianceTables$`dimension = 5`}. Similar calls are made to retrieve the variance table for other factor solutions.}
}
\references{
Ledermann, W. (1937). On the rank of the reduced correlational matrix in multiple factor analysis. Psychometrika, 2:85--93.

Peeters, C.F.W. \emph{et al.} (2019). Stable prediction with radiomics data.
}
\author{
Carel F.W. Peeters <cf.peeters@vumc.nl>
}
\note{
\itemize{
\item The argument \code{maxdim} cannot exceed the Ledermann-bound (Ledermann, 1937): \eqn{\lfloor [2p + 1 - (8p + 1)^{1/2}]/2\rfloor}, where \eqn{p} indicates the observed-feature dimension.
Usually, one wants to set \code{maxdim} much lower than this bound.

\item The tabulations in the \code{$varianceTables} slot are based on unrotated maxmimum likelihood factor solutions. Note that the total cumulative variance does not depend on the choice of (orthogonal) rotation.
}
}
\seealso{
\code{\link{dimGB}}, \code{\link{FAsim}}, \code{\link{mlFA}}, \code{\link{SMC}}
}
\examples{
## Simulate some high-dimensional data according to the factor model
simDAT <- FAsim(p = 50, m = 5, n = 40)

## Regularize the correlation matrix
RegR <- regcor(simDAT$data)

## Assess proportion and cumulative variances for a range of factor solutions
## Inspect, for example, the variance table for the 5-factor solution
V <- dimVAR(RegR$optCor, maxdim = 20)
V$varianceTables$`dimension = 5`
}
\concept{Cumulative variance}
