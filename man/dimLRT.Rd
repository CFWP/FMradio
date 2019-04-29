\name{dimLRT}
\alias{dimLRT}
\title{
Assess the latent dimensionality using a likelihood ratio test
}
\description{
\code{dimLRT} is a function that evaluates a likelihood ratio test on the factor model.
It can be used to choose the number of latent factors.
}
\usage{
dimLRT(R, X, maxdim, rankDOF = TRUE, graph = TRUE, 
       alpha = .05, Bartlett = FALSE, verbose = TRUE)
}
\arguments{
  \item{R}{
  (Regularized) correlation \code{matrix}.
  }
  \item{X}{
  A (possibly centered and scaled and possibly subsetted) data \code{matrix}.
  }
  \item{maxdim}{
  A \code{numeric} integer or \code{integer} indicating the maximum factor dimension to be assessed.
  }
  \item{rankDOF}{
  A \code{logical} indicating if the degrees of freedom should be based on the rank of the raw correlation matrix.
  }
  \item{graph}{
  A \code{logical} indicating if the results should be visualized.
  }
  \item{alpha}{
  A \code{numeric} scalar representing the alpha level. Only used when \code{graph = TRUE}.
  }
  \item{Bartlett}{
  A \code{logical} indicating if the Bartlett correction should be applied.
  }
  \item{verbose}{
  A \code{logical} indicating if the function should run silently.\cr
  Runs silently when \code{verbose = FALSE}.
  }
}
\details{
The most formal approach to factor analytic dimensionality assessment is through likelihood ratio (LR) testing.
The basic idea is to test the \eqn{m}-factor model against the saturated model.
The corresponding LR criterion then converges, under the standard correlation matrix and corresponding parameter estimates under \eqn{m}-factors, to \eqn{(n - 1)} times a certain discrepancy function evaluated at the maximum-likelihood-parameters under the \eqn{m}-factor model.
This quantity is approximately \eqn{\chi^{2}}-distributed under certain regularity conditions (Amemiya & Anderson, 1990).
The general strategy would then be to sequentially test solutions of increasing dimensionality \eqn{m = 1, \ldots, \mbox{maxdim}} until the null hypothesis (stating that the \eqn{m}-factor model holds) is \emph{not} rejected at Type-I error level \code{alpha}.

The degrees of freedom for the LRT under the \eqn{m}-factor model equals the number of parameters in the saturated model (i.e., the unstructured sample correlation) minus the number of freely estimable parameters in the \eqn{m}-factor model.
Note that the general stategy above makes use of asymptotic results.
In our setting, however, the observation dimension (\eqn{n}) is usually small relative to the feature dimension (\eqn{p}).
Hence, the standard test will in a sense overestimate the degrees of freedom.
One simple option dealing with this observation would be to adapt the degrees of freedom to incorporate the rank deficiency of \code{R}.
This road is taken when \code{rankDOF = TRUE}.
Bartlett (1950) proposed a correction factor when the sample size is small to make the test statistic behave more \eqn{\chi^{2}}-like.
This correction factor is used when \code{Bartlett = TRUE}.

When \code{graph = TRUE} the LRT results are visualized. 
The graph plots the LRT \eqn{p}-values against the consecutive dimensions of the factor solution.
A horizontal line is plotted at the value provided in the \code{alpha} argument.

Unless the number of observations is much larger than the number of features, the LRT is not recommended for inference in general.
In Peeters \emph{et al.} (2019) the LRT was assessed in a comparative setting inviolving high-dimensional factor models.
}
\value{
The function returns an object of class \code{data.frame}.
The first column represents the assessed dimensions running from 1 to \code{maxdim}.
The second column represents the observed values of the LRT statistic.
The third column represents the corresponding \eqn{p}-values.
}
\references{
Amemiya, Y., & Anderson, T.W. (1990). Asymptotic chi-square tests for a large class of factor
analysis models. The Annals of Statistics, 18:1453--1463.

Bartlett, M.S. (1950). Tests of significance in factor analysis. British Journal of Psychology (Statistics Section), 3:77--85.

Ledermann, W. (1937). On the rank of the reduced correlational matrix in multiple factor analysis. Psychometrika, 2:85--93.

Peeters, C.F.W. \emph{et al.} (2019). Stable prediction with radiomics data.
\href{https://arxiv.org/abs/1903.11696}{arXiv:1903.11696 [stat.ML]}.
}
\author{
Carel F.W. Peeters <cf.peeters@vumc.nl>, Caroline Ubelhor
}
\note{
\itemize{
\item Note that, for argument \code{X}, the observations are expected to be in the rows and the features are expected to be in the columns.

\item The argument \code{maxdim} cannot exceed the Ledermann-bound (Ledermann, 1937): \eqn{\lfloor [2p + 1 - (8p + 1)^{1/2}]/2\rfloor}, where \eqn{p} indicates the observed-feature dimension.
Usually, one wants to set \code{maxdim} much lower than this bound.

\item note that, if \eqn{p > n}, then the maximum rank of the raw correlation matrix is \eqn{n - 1}. In this case there is an alternative Ledermann-bound when \code{rankDOF = TRUE}. The number of information points in the correlation matrix is then given as \eqn{n\times (n-1)/2} and this number must exceed \eqn{p\times \mbox{maxdim} + p - (\mbox{maxdim} \times (\mbox{maxdim} - 1))/2}, putting more restrictions on \code{maxdim}.

\item Other functions for factor analytic dimensionality assessment are \code{\link{dimGB}} and \code{\link{dimIC}}. In high-dimensional situations usage of \code{\link{dimGB}} on the regularized correlation matrix is recommended.
}
}
\seealso{
\code{\link{dimGB}}, \code{\link{FAsim}}
}
\examples{
## Simulate some data according to the factor model
## $cormatrix gives the correlation matrix on the generated data
simDAT <- FAsim(p = 50, m = 5, n = 500)
simDAT$cormatrix

## Calculate the LRT for models of factor dimension 1 to 20
LRT <- dimLRT(simDAT$cormatrix, simDAT$data, maxdim = 20, rankDOF = FALSE)
print(LRT)
}
\concept{Factor analytic likelihood ratio testing}
