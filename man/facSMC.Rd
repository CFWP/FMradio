\name{facSMC}
\alias{facSMC}
\title{
Evaluate the determinacy of factor scores
}
\description{
\code{facSMC} is a function with which one may evaluate the determinacy of factor scores.
}
\usage{
facSMC(R, LM)
}
\arguments{
  \item{R}{
  (Regularized) correlation \code{matrix}.
  }
  \item{LM}{
  A (rotated) loadings \code{matrix}. Usually the \code{$Loadings}-slot object from the \code{\link{mlFA}} function output.
  }
}
\details{
The \code{facSMC} function calculates the squared multiple correlations between the observed features and the common latent factors.
The closer to unity, the lesser the problem of factor-score indeterminacy and the better one is able to uniquely determine the  factor scores.
In practice, a squared multiple correlation equalling or exceeding .9 would be considered adequate.
See Mulaik (2010, Chapter 13) and Peeters \emph{et al.} (2019, Supplementary Materials) for further details.
}
\value{
The function returns a \code{numeric} vector indicating, for each factor, the squared multiple correlation between the observed features and the common latent factor.
}
\references{
Mulaik, S.A. (2010). Foundations of Factor Analysis. Boca Raton: Chapman & Hall/CRC, 2nd
edition.

Peeters, C.F.W. \emph{et al.} (2019). Stable prediction with radiomics data.
\href{https://arxiv.org/abs/1903.11696}{arXiv:1903.11696 [stat.ML]}.
}
\author{
Carel F.W. Peeters <cf.peeters@vumc.nl>
}
\note{
Note that the computations assume an orthogonal factor model.
Hence, only orthogonal rotations of the loadings matrix should be used (or no rotation at all).
}
\seealso{
\code{\link{facScore}}
}
\examples{
## Simulate some data according to a factor model with 5 latent factors
## Simulate high-dimensional situation in the sense that p > n
## $cormatrix gives the correlation matrix on the generated data
simDAT <- FAsim(p = 50, m = 5, n = 40, loadingvalue = .9)
simDAT$cormatrix

## Regularize the correlation matrix
RegR <- regcor(simDAT$data)

## Evaluate the Guttman bounds
## First Guttman bound indicates to retain 5 latent factors
GB <- dimGB(RegR$optCor)
print(GB)

## Produce ML factor solution under 5 factors
## Print loadings structure of this solution
fit <- mlFA(RegR$optCor, 5)
print(fit$Loadings, digits = 2, cutoff = .3, sort = TRUE)

## Obtain factor-scores
scores <- facScore(scale(simDAT$data), fit$Loadings, fit$Uniqueness)
print(scores)

## Evaluate determinacy of factor scores
fd <- facSMC(RegR$optCor, fit$Loadings)
print(fd)
}
\concept{Factor determinacy}
