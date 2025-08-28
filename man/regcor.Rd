\name{regcor}
\alias{regcor}
\title{
Regularized correlation matrix estimation
}
\description{
\code{regcor} is a function that determines the optimal penalty value and, subsequently, the optimal Ledoit-Wolf type regularized correlation matrix using K-fold cross validation of the negative log-likelihood.
}
\usage{
regcor(X, fold = 5, verbose = TRUE)
}
\arguments{
  \item{X}{
  A (possibly centered and scaled and possibly subsetted) data \code{matrix}.
  }
  \item{fold}{
  A \code{numeric} integer or \code{integer} indicating the number of folds to use in cross-validation.
  }
  \item{verbose}{
  A \code{logical} indicating if function should run silently.\cr
  Runs silently when \code{verbose = FALSE}.
  }
}
\details{
This function estimates a Ledoit-Wolf-type (Ledoit & Wolf, 2004) regularized correlation matrix.
The optimal penalty-value is determined internally by \emph{K}-fold cross-validation of the of the negative log-likelihood function.
The procedure is efficient as it makes use of the Brent root-finding procedure (Brent, 1971).
The value at which the \emph{K}-fold cross-validated negative log-likelihood score is minimized is deemed optimal. 
The function employs the Brent algorithm as implemented in the \href{https://stat.ethz.ch/R-manual/R-devel/library/stats/html/optim.html}{optim} function.
It outputs the optimal value for the penalty parameter and the regularized correlation matrix under this optimal penalty value.
See Peeters \emph{et al.} (2019) for further details.

The optimal penalty-value can be used to assess the conditioning of the estimated regularized correlation matrix using, for example, a condition number plot (Peeters, van de Wiel, van Wieringen, 2016). 
The regularized correlation matrix under the optimal penalty can serve as the input to functions that assess factorability (\code{\link{SA}}), evaluate optimal choices of the latent common factor dimensionality (e.g., \code{\link{dimGB}}), and perform maximum likelihood factor analysis (\code{\link{mlFA}}).
}
\value{
The function returns an object of class \code{list}:
\item{$optPen}{A \code{numeric} scalar representing the optimal value for the penalty parameter.}
\item{$optCor}{A \code{matrix} representing the regularized correlation matrix under the optimal penalty-value.}
}
\references{
Brent, R.P. (1971). An Algorithm with Guaranteed Convergence for Finding a Zero of a Function. Computer Journal 14: 422--425.

Ledoit, O, & Wolf, M. (2004). A well-conditioned estimator for large-dimensional covariance matrices. Journal of Multivariate Analysis, 88:365--411.

Peeters, C.F.W. \emph{et al.} (2019). Stable prediction with radiomics data.
\href{https://arxiv.org/abs/1903.11696}{arXiv:1903.11696 [stat.ML]}.

Peeters, C.F.W., van de Wiel, M.A., & van Wieringen, W.N. (2016). The spectral condition number plot for regularization parameter determination, arXiv:1608.04123v1 [stat.CO].
}
\author{
Carel F.W. Peeters <carel.peeters@wur.nl>
}
\note{
Note that, for argument \code{X}, the observations are expected to be in the rows and the features are expected to be in the columns.
}
\seealso{
\code{\link{RF}}, \code{\link{subSet}}, \code{\link{SA}}, \code{\link{dimGB}}, \code{\link{mlFA}}
}
\examples{
## Generate some (high-dimensional) data
## Get correlation matrix
p = 25
n = 10
set.seed(333)
X = matrix(rnorm(n*p), nrow = n, ncol = p)
colnames(X)[1:25] = letters[1:25]
R <- cor(X)

## Redundancy visualization, at threshold value .9
radioHeat(R, diag = FALSE, threshold = TRUE, threshvalue = .9)

## Redundancy-filtering of correlation matrix
Rfilter <- RF(R, t = .9)
dim(Rfilter)

## Subsetting data
DataSubset <- subSet(X, Rfilter)
dim(DataSubset)

## Obtain regularized correlation matrix
RegR <- regcor(DataSubset, fold = 5, verbose = TRUE)
RegR$optPen  ## optimal penalty-value
}
\concept{regularized correlation}
