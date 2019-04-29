\name{facScore}
\alias{facScore}
\title{
Compute factor scores
}
\description{
\code{facScore} is a function that computes factor scores, the score each person/object attains on each latent factor.
}
\usage{
facScore(X, LM, UM, type = "thomson")
}
\arguments{
  \item{X}{
  A (scaled and possibly subsetted) data \code{matrix}.
  }
  \item{LM}{
  A (rotated) loadings \code{matrix}. Usually the \code{$Loadings}-slot object from the \code{\link{mlFA}} function output.
  }
  \item{UM}{
  A diagonal uniquenesses \code{matrix}. Usually the \code{$Uniqueness}-slot object from the \code{\link{mlFA}} function output.
  }
  \item{type}{
  A \code{character} indicating the type of factor score to calculate.\cr
  Must be one of: "thomson", "bartlett", "anderson".
  }
}
\details{
Once a factor model is fitted one may desire an estimate of the score each object/individual would obtain on each of the latent factors.
Such scores are referred to as factor scores.
The \code{facScore} function provides several types of factor score estimates.
The default are Thomson-type scores (Thomson, 1939).
These may be viewed as (empirical) Bayesian-type scores.
Bartlett-type scores (Bartlett, 1937) are unbiased but less efficient in terms of mean-squared error.
Under the orthogonal model the latent factors are orthogonal in
the population and, hence, the Thomson and Bartlett-type factor scores will be near orthogonal
in the sample. 
Anderson and Rubin (1956) constructed an alternative estimator for the factor
scores that enforces their orthogonality in the sample.
}
\value{
The function returns an object of class \code{data.frame}.
Observations are represented in the rows.
Each column represent a latent factor.
}
\references{
Anderson, T.W., & Rubin, H. (1956). Statistical inference in factor analysis. In Proceedings
of the Third Berkeley Symposium on Mathematical Statistics and Probability, volume 5:
Contributions to Econometrics, Industrial Research, and Psychometry, pages 111--150.
Berkeley, CA: University of California Press.

Bartlett, M.S. (1937). The statistical conception of mental factors. British Journal of Psychology,
28:97--104.

Peeters, C.F.W. \emph{et al.} (2019). Stable prediction with radiomics data.
\href{https://arxiv.org/abs/1903.11696}{arXiv:1903.11696 [stat.ML]}.

Thomson, G. (1939). The Factorial Analysis of Human Ability. London: University of Londen
Press.
}
\author{
Carel F.W. Peeters <cf.peeters@vumc.nl>
}
\note{
The input data (argument \code{X}) are assumed to be scaled (or at least centered).
The \code{UM} matrix is assumed to be positive definite.
The \code{LM} matrix is assumed to be of full column rank.
}
\seealso{
\code{\link{dimGB}}, \code{\link{mlFA}}, \code{\link{facSMC}}
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
}
\concept{Factor scores}
