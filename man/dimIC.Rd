\name{dimIC}
\alias{dimIC}
\title{
Assess the latent dimensionality using information criteria
}
\description{
A function that calculates either the AIC or the BIC on the factor model.
These can be used to choose the number of latent factors.
}
\usage{
dimIC(R, n, maxdim, Type = "BIC", graph = TRUE, verbose = TRUE)
}
\arguments{
  \item{R}{
  (Regularized) correlation \code{matrix}.
  }
  \item{n}{
  A \code{numeric} scalar representing the sample size.
  }
  \item{maxdim}{
  A \code{numeric} integer or \code{integer} indicating the maximum factor dimension to be assessed.
  }
  \item{Type}{
  A \code{character} indicating the type of IC to be calculated.\cr
  Must be one of: "AIC", "BIC".
  }
  \item{graph}{
  A \code{logical} indicating if the results should be visualized.
  }
  \item{verbose}{
  A \code{logical} indicating if the function should run silently.\cr
  Runs silently when \code{verbose = FALSE}.
  }
}
\details{
Information criteria (IC) are often used in selecting the number of latent factor to retain.
IC aim to balance model fit with model complexity. 
They evaluate (minus 2 times) the maximized value of the (model-dependent) likelihood function weighed with a penalty function that is dependent on the free parameters in the model. 
Different penalizations define the different types of IC.
The strategy would be to determine IC scores for a range of consecutive values of the latent factor dimension.
This function then determines scores for factor solutions ranging from 1 to \code{maxdim} latent factors.
The solution with the lowest IC score is deemed optimal.
The function allows for the calculation of either the Akaike information criterion (AIC; Akaike, 1973) or the Bayesian information criterion (BIC; Schwarz, 1978).
Also see the Supplementary Material of Peeters \emph{et al.} (2019) for additional detail.

When \code{graph = TRUE} the IC scores are visualized. 
The graph plots the IC score against the consecutive dimensions of the factor solution.
}
\value{
The function returns an object of class \code{data.frame}.
The first column represents the assessed dimensions running from 1 to \code{maxdim}.
The second column represents the corresponding values of the chosen information criterion.
}
\references{
Akaike, H. (1973). Information theory and an extension of the maximum likelihood principle. In: B. N. Petrov and F. Csaki (Eds.) Second International Symposium on Information Theory, pages 267-281. Budapest: Akademiai Kaido.

Ledermann, W. (1937). On the rank of the reduced correlational matrix in multiple factor analysis. Psychometrika, 2:85-93.

Peeters, C.F.W. \emph{et al.} (2019). Stable prediction with radiomics data.

Schwarz, G.E. (1978). Estimating the dimension of a model. Annals of Statistics, 6:461-464.
}
\author{
Carel F.W. Peeters <cf.peeters@vumc.nl>
}
\note{
\itemize{
\item The argument \code{maxdim} cannot exceed the Ledermann-bound (Ledermann, 1937): \eqn{\lfloor [2p + 1 - (8p + 1)^{1/2}]/2\rfloor}, where \eqn{p} indicates the observed-feature dimension.
Usually, one wants to set \code{maxdim} much lower than this bound.

\item Other functions for factor analytic dimensionality assessment are \code{\link{dimGB}} and \code{\link{dimLRT}}. In high-dimensional situations usage of \code{\link{dimGB}} on the regularized correlation matrix is recommended.
}
}
\seealso{
\code{\link{dimGB}}, \code{\link{FAsim}}
}
\examples{
## Simulate some data according to the factor model
## $cormatrix gives the correlation matrix on the generated data
simDAT <- FAsim(p = 50, m = 5, n = 100)
simDAT$cormatrix

## Calculate the AIC for models of factor dimension 1 to 20
AIC <- dimIC(simDAT$cormatrix, n = 100, Type = "AIC", maxdim = 20)
print(AIC)

## Calculate the BIC for models of factor dimension 1 to 20
BIC <- dimIC(simDAT$cormatrix, n = 100, Type = "BIC", maxdim = 20)
print(BIC)
}
\concept{Information criteria}
