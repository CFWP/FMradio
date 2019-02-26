\name{regcor}
\alias{regcor}
\title{
Regularized correlation matrix estimation
}
\description{
Function determining the optimal penalty value and, subsequently, the optimal Ledoit-Wolf type regularized correlation matrix using K-fold cross validation of the negative log-likelihood.
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
  A \code{logical} indicating if function should run silently.
  Runs silently when \code{verbose = FALSE}.
  }
}
\details{
Regularized correlation matrix. (LW-type: ref to LW.)
How do we determine the penalty value.

also input into dimGB and SA functions.
These, together with this output is Input into mlFA.

  # - The estimator is a Ledoit-Wolf type estimator
  # - The K-fold cv procedure makes use of the Brent algorithm
}
\value{
The function returns an object of class \code{list}:
\item{$optPen}{A \code{numeric} scalar representing the optimal value for the penalty parameter.}
\item{$optCor}{A \code{matrix} representing the regularized correlation matrix under the optimal penalty-value.}
}
\references{
Peeters, C.F.W. \emph{et al.} (2019). Stable prediction with radiomics data.
}
\author{
Carel F.W. Peeters <cf.peeters@vumc.nl>
}
\note{
Note that, for argument \code{X}, the observations are expected to be in the rows and the features are expected to be in the columns.
}
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
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
