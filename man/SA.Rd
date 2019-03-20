\name{SA}
\alias{SA}
\title{
Calculate the KMO measure of feature-sampling adequacy
}
\description{
\code{SA} is a function that calculates the Kaiser-Meyer-Olkin (KMO) measure of sampling adequacy.
}
\usage{
SA(R)
}
\arguments{
  \item{R}{
  (Regularized) covariance or correlation \code{matrix}.
  }
}
\details{
The \code{SA} function calculates the Kaiser-Meyer-Olkin (KMO) measure of feature-sampling adequacy (Kaiser & Rice, 1974).
It provides a practical option for the assessment of factorability.
Factorability refers to the assessment of the ability to identify coherent common latent factors from a given correlation matrix.
In common factor analysis the observed features are assumed to be independent \emph{given} the common latent features.
Under this crucial model assumption, the inverse of the population correlation matrix is diagonal.
Hence, to assess factorability one could assess if the inverse of the sample correlation matrix is near-diagonal.
The KMO index provides for such an assessment by "comparing the sizes of the off-diagonal entries of the regularized correlation matrix to the sizes of the off-diagonal entries of its scaled inverse" (Peeters \emph{et al.}, 2019). 
It takes values in \eqn{[0,1]} and larger values are preferred.
A KMO index between .9 and 1 would be considered to be indicative of great factorability.
For rules of thumb regarding interpretation of KMO index value, see Kaiser (1970).
The \code{SA} function calculates an overall KMO index as well as the KMO index per observed feature.

The intended usage of the \code{SA} function is to assess if performing a factor analysis on a given (regularized) correlation matrix can be considered appropriate. As such, it succeeds usage of the \code{\link{regcor}} function (for high-dimensional and/or strongly collinear settings) and precedes 
usage of the \code{\link{dimGB}} and \code{\link{mlFA}} functions.
}
\value{
The function returns an object of class \code{list}:
\item{$KMO}{A \code{numeric} scalar representing the overall KMO index.}
\item{$KMOfeature}{A \code{numeric} vector giving the KMO index per feature.}
}
\references{
Kaiser, H.F. (1970). A second-generation little jiffy. Psychometrika, 35:401--415.

Kaiser, H.F., & Rice., J. (1974). Little jiffy, mark IV. Educational and Pscyhological Measurement,
34:111--117.

Peeters, C.F.W. \emph{et al.} (2019). Stable prediction with radiomics data.
}
\author{
Carel F.W. Peeters <cf.peeters@vumc.nl>
}
\note{
The input \code{matrix} \code{R} should be nonsingular for the KMO to be computed.
When \code{R} is singular one may regularize it using the \code{\link{regcor}} function.
}

\seealso{
\code{\link{regcor}}, \code{\link{dimGB}}, \code{\link{mlFA}}
}
\examples{
## Generate some (high-dimensional) data
p = 25
n = 10
set.seed(333)
X = matrix(rnorm(n*p), nrow = n, ncol = p)
colnames(X)[1:25] = letters[1:25]

## Obtain regularized correlation matrix
RegR <- regcor(X, fold = 5, verbose = TRUE)

## Assess factorability through the KMO index
factorable <- SA(RegR$optCor)
factorable$KMO
factorable$KMOfeature
}
\concept{KMO-index}
