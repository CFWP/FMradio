\name{SMC}
\alias{SMC}
\title{
Compare squared multiple correlations with model-based communalities
}
\description{
\code{SMC} is a function that compares the best lower-bound estimates to the communalities with the model-based communalities implied by a factor solution of dimension \eqn{m}.
}
\usage{
SMC(R, LM)
}
\arguments{
  \item{R}{
  (Regularized) correlation \code{matrix}.
 }
  \item{LM}{
  (Rotated) factor loadings \code{matrix}.
 }
}
\details{
This function can be used to qualitatively assess the choice of dimensionality (as well as the fit) in the \eqn{m}-factor model.
This is done using the concept of communalities.
The communality refers to the amount of variance of feature \eqn{j} explained by the latent features.
It is then of interest to compare lower-bound estimates of the (population) communalities to the extracted communalities under the \eqn{m}-factor model.

Guttman (1956) gave the best possible lower-bound estimates to the communalities, which can essentially be considered squared multiple correlations: the proportion of variance in feature \eqn{j} that is explained by the remaining \eqn{p - 1} features.
To assess a factor model, these might be compared to the retrieved estimated communalities under the \eqn{m}-factor model.
When the chosen latent dimensionality is sufficient then one would expect that, for almost all features, the retrieved communality approximately equals or exceeds its corresponding lower-bound estimate.
If this is not the case then one might have extracted too few factors.
}
\value{
The function returns a \code{matrix}.
The first column (labeled 'SMC') contains the lower-bound estimates to the communalities.
The second column (labeled 'Communalities') contains the retrieved estimated communalities under the \eqn{m}-factor model.
}
\references{
Guttman, L. (1956). Best possible systematic estimates of communalities. Psychometrika, 21:273--285.

Peeters, C.F.W. \emph{et al.} (2019). Stable prediction with radiomics data.
\href{https://arxiv.org/abs/1903.11696}{arXiv:1903.11696 [stat.ML]}.
}
\author{
Carel F.W. Peeters <carel.peeters@wur.nl>
}
\note{
Note that the choice of orthogonal rotation does not affect the model-implied communality estimates.
}
\seealso{
\code{\link{dimGB}}, \code{\link{FAsim}}, \code{\link{mlFA}}, \code{\link{dimVAR}}
}
\examples{
## Simulate some high-dimensional data according to the factor model
simDAT <- FAsim(p = 50, m = 5, n = 40)

## Regularize the correlation matrix
RegR <- regcor(simDAT$data)

## Fit 5-factor model to the regularized correlation matrix
fit <- mlFA(RegR$optCor, m = 5)

## Compare lower-bound estimates to communalities with model-implied ones
C <- SMC(RegR$optCor, fit$Loadings)
print(C)
}
\concept{Squared multiple correlation}
