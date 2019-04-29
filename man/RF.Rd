\name{RF}
\alias{RF}
\title{
Redundancy filtering of a square (correlation) matrix
}
\description{
\code{RF} is a function that performs redundancy filtering (RF) of a square (correlation) matrix.
}
\usage{
RF(R, t = .95)
}
\arguments{
  \item{R}{
  Square (correlation) \code{matrix}.
  }
  \item{t}{
  A scalar \code{numeric} indicating the absolute value for thresholding.
  }
}
\details{
Radiomic features can be very strongly correlated.
The sample correlation matrix on extracted radiomic features will then often display strong collinearity.
The collinearity may be so strong as to imply redundant information, in the sense that some entries will approach perfect (negative) correlation. 
Hence, one may wish to perform redundancy-filtering on the raw sample correlation matrix in such situations.

The \code{RF} function uses an Algorithm from Peeters \emph{et al.} (2019) to remove the minimal number of redundant features under absolute marginal correlation threshold \code{t}.
We recommend setting \eqn{\mathrm{t} \in [.9,.95]}.
Details of the algorithm can be found in Peeters \emph{et al.} (2019).

The function returns a redundancy-filtered correlation \code{matrix}.
This return output may subsequently be used in the \code{\link{subSet}} function.
This is a convenience function that subsets a dataset to the features retained after redundancy-filtering.
}
\value{
Returns a redundancy-filtered \code{matrix}.
}
\references{
Peeters, C.F.W. \emph{et al.} (2019). Stable prediction with radiomics data.
\href{https://arxiv.org/abs/1903.11696}{arXiv:1903.11696 [stat.ML]}.
}
\author{
Carel F.W. Peeters <cf.peeters@vumc.nl>
}
\note{
\itemize{
\item While geared towards the redundancy filtering of correlation matrices, the function is quite general, in the sense that it can be used to filter any square \code{matrix}.
\item When the input matrix \code{R} is a correlation matrix, then argument \code{t} should satisfy \eqn{-1 < \mathrm{t} < 1}, for the return matrix to be sensical for further analysis.
}
}
\seealso{
\code{\link{subSet}}, \code{\link{regcor}}
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
}
\concept{redundancy-filtering}
