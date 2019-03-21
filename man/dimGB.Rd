\name{dimGB}
\alias{dimGB}
\title{
Assess the latent dimensionality using Guttman bounds
}
\description{
\code{dimGB} is a function that calculates the first, second, and third Guttman (lower-)bounds to the dimensionality of the latent vector.
These can be used to choose the number of latent factors.
}
\usage{
dimGB(R, graph = TRUE, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{R}{
  (Regularized) correlation \code{matrix}.
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
The communality in factor analysis refers to the amount of variance (of feature \eqn{j}) explained by the latent features.
The correlation of any feature with itself can then be decomposed into common variance (the communality) and unique variance.
This implies that unity (1) minus the unique variance for feature \eqn{j} equals the communality for feature \eqn{j}.
From the matrix perspective one can then construct a reduced correlation matrix: the correlation matrix with communalities in the diagonal.
This reduced correlation matrix is, by the assumptions on the factor model, Gramian and of rank \eqn{m}, with \eqn{m} indicating the intrinsic dimensionality of the latent vector.
The dimension of the latent vector (i.e., the number of common factors) can then be assessed by evaluating the rank of the sample correlation matrix in which the diagonal elements are replaced with appropriate communality estimates.

In our case, which is often high-dimensional, we use the regularized correlation matrix as our sample-representation of the population correlation matrix.
The diagonal elements are then replaced with Guttman's lower-bound estimates for the communalities (Guttman, 1956).
Guttman (1956) gives 3 (ordered) lower-bound estimates.
The first estimate is the most conservative, using 0 as a lower-bound estimate of the communalities.
From this perspective, every positive eigenvalue of the reduced sample correlation matrix is indicative of a latent factor whose contribution to variance-explanation is above and beyond mere unique variance. 
The decisonal approach would then be to retain all such factors.
See Peeters \emph{et al.} (2019) for additional detail.

The Guttman approach has historically been used as a lower-bound estimate of the latent dimensionality.
We consider the decisional approach stated above to give an upper-bound.
Peeters \emph{et al.} (2019) contains an extensive simulation study showing that in high-dimensional situations this decisional approach provides a reliable upper-bound.
The choice of the number of factors can be further assessed with the \code{\link{SMC}} and \code{\link{dimVAR}} functions.
Assessments provided by these latter functions may inform if the result of the decisional rule above should be accepted or be treated as an upper-bound.

When \code{graph = TRUE} the Guttman bounds are visualized. 
It plots the consecutive eigenvalues for each of the reduced correlation matrices.
The number of positive eigenvalues for each respective reduced correlation matrix then corresponds to each of the respective Guttman bounds.
The visualization may be of limited value when the feature-dimension gets (very) large.
}
\value{
The function returns an object of class \code{table}.
The entries correspond to the first, second, and third Guttman bounds.
}
\references{
Guttman, L. (1956). Best possible systematic estimates of communalities. Psychometrika, 21:273--285.

Peeters, C.F.W. \emph{et al.} (2019). Stable prediction with radiomics data.
}
\author{
Carel F.W. Peeters <cf.peeters@vumc.nl>
}
\note{
\itemize{
\item Again, from a historical perspective, the decisional rule would have been used as a lower-bound to the question of the number of latent common factors. In high-dimensional situations we recommend to use it as an upper-bound.

\item Other functions for factor analytic dimensionality assessment are \code{\link{dimIC}} and \code{\link{dimLRT}}. In high-dimensional situations usage of \code{dimGB} is recommended over these other functions.
}
}
\seealso{
\code{\link{SMC}}, \code{\link{dimVAR}}, \code{\link{FAsim}}
}
\examples{
## Simulate some data according to a factor model with 5 latent factors
## $cormatrix gives the correlation matrix on the generated data
simDAT <- FAsim(p = 50, m = 5, n = 100)
simDAT$cormatrix

## Evaluate the Guttman bounds
## First Guttman bound indicates to retain 5 latent factors
GB <- dimGB(simDAT$cormatrix)
print(GB)
}
\concept{Guttman bounds}
