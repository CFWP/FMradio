\name{subSet}
\alias{subSet}
\title{
Subset a data matrix or expression set
}
\description{
\code{subSet} is a convenience function that subsets a data \code{matrix} or an \code{ExpressionSet} object. 
}
\usage{
subSet(X, Rf)
}
\arguments{
  \item{X}{
  A data \code{matrix} or an \code{ExpressionSet} object. 
  }
  \item{Rf}{
  A filtered (correlation) \code{matrix} (as returned by the \code{\link{RF}} function).
  }
}
\details{
The \code{subSet} convenience function may directly follow usage of the \code{\link{RF}} in the sense that the latters return-value can be used as the \code{Rf} argument.
It then subsets a data \code{matrix} or an \code{ExpressionSet} object to those features retained by the redundancy-filtering.
The function returns a subsetted \code{matrix} or \code{ExpressionSet} (depending on the class of the \code{X} argument).
The subsetted data can then be used for penalty-parameter selection and regularized correlation matrix estimation provided by the \code{\link{regcor}} function.
}
\value{
Returns a subsetted data \code{matrix} or \code{ExpressionSet}.
}
\references{
Peeters, C.F.W. \emph{et al.} (2019). Stable prediction with radiomics data.
}
\author{
Carel F.W. Peeters <cf.peeters@vumc.nl>
}
\note{
If argument \code{X} is a \code{matrix}, the observations are expected to be in the rows and the features are expected to be in the columns.
}
\seealso{
\code{\link{regcor}}
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
}
\concept{data subsetting}
