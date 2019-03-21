\name{mlFA}
\alias{mlFA}
\title{
Maximum likelihood factor analysis
}
\description{
\code{mlFA} is a function that performs a maximum likelihood factor analysis.
}
\usage{
mlFA(R, m)
}
\arguments{
  \item{R}{
  (Regularized) correlation \code{matrix}.
  }
  \item{m}{
  A \code{numeric} integer or \code{integer} indicating the latent dimension of the factor solution (i.e., the number of factors).
  }
}
\details{
This function is basically a wrapper around the 
\href{https://stat.ethz.ch/R-manual/R-devel/library/stats/html/factanal.html}{factanal} function from the \href{https://stat.ethz.ch/R-manual/R-devel/library/stats/html/stats-package.html}{stats} package.
Its purpose is to produce a factor solution of the chosen dimension (argument \code{m}) by a maximum likelihood estimation procedure (Joreskog, 1967).
The wrapper ensures that the model is fitted under the same circumstances under which latent dimensionality is assessed with functions such as \code{\link{dimLRT}} and \code{\link{dimIC}}.
The function produces a Varimax rotated (Kaiser, 1958) factor solution.
The output can be used to produce factor scores by the \code{\link{facScore}} function.
}
\value{
The function returns an object of class \code{list}:
\item{$Loadings}{A matrix of class \code{loadings} representing the loadings matrix in which in which each element \eqn{\lambda_{jk}} is the loading of the \eqn{j}th feature on the \eqn{k}th latent factor}
\item{$Uniqueness}{A \code{matrix} representing the diagonal matrix carrying the unique variances.}
\item{$rotmatrix}{A \code{matrix} representing the Varimax rotation matrix.}
}
\references{
Joreskog, K.G (1967). Some contributions to maximum likelihood factor analysis. Psychometrika,
32:443--482.

Kaiser, H.F. (1958). The varimax criterion for analytic rotation in factor analysis. Psychometrika,
23:187--200.

Peeters, C.F.W. \emph{et al.} (2019). Stable prediction with radiomics data.
}
\author{
Carel F.W. Peeters <cf.peeters@vumc.nl>
}
\note{
\itemize{
\item Note that the order of the features in the \code{$Loadings} and \code{$Uniqueness} slots of the output is determined by the order of the features for the input argument \code{R}. As the \code{$Loadings} slot gives an object of class "loadings" it can be subjected to the \code{print} function, which sorts the output to emphasize the loadings structure when calling \code{sort = TRUE}.

\item Note that the maximum likelihood procedure is stable when a regularized correlation matrix is used as the input for argument \code{R}.

\item In high-dimensional situations usage of \code{\link{dimGB}} on the regularized correlation matrix is recommended to determine the value for argument \code{m}.
}
}
\seealso{
\code{\link{dimGB}}, \code{\link{facScore}}
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
}
\concept{Factor analysis}
