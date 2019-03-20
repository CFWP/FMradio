\name{scree}
\alias{scree}
\title{
Scree plot for an \eqn{m}-factor solution
}
\description{
\code{scree} is a function that produces a scree plot under an \eqn{m}-factor solution.
}
\usage{
scree(R, m, type = "ML", vertical = TRUE, verbose = TRUE)
}
\arguments{
  \item{R}{
  (Regularized) correlation \code{matrix}.
  }
  \item{m}{
  A \code{numeric} integer or \code{integer} indicating the dimension of the factor solution.
  }
  \item{type}{
  A \code{character} indicating the type of eigenvalues to be represented in the scree plot.\cr
  Must be one of: "ML", "reduced".
  }
  \item{vertical}{
  A \code{logical} indicating if  the value of argument \code{m} should be indicated in the scree plot by a vertical line.
  }
  \item{verbose}{
  A \code{logical} indicating if the function should run silently.\cr
  Runs silently when \code{verbose = FALSE}.
  }
}
\details{
Which types of eigenvalues plotted.
How should it be used.

1 before the elbow

  # - If type = ML, the eigenvalues under the m-factor solution are returned
  # - If type = reduced, the eigenvalues of the reduced correlation matrix
  #   under the m-factor solution are returned
  
Are as much eigenvalues as features.
}
\value{
The function returns an object of class \code{data.frame}.
The first column represents the assessed latent dimensions running from 1 to \code{maxdim}.
The second column represents the observed values of the LRT statistic.
}
\references{
Cattell, R.B. (1966). The Scree Test For The Number Of Factors. Multivariate Behavioral Research. 1: 245--276.

Ledermann, W. (1937). On the rank of the reduced correlational matrix in multiple factor analysis. Psychometrika, 2:85--93.

Peeters, C.F.W. \emph{et al.} (2019). Stable prediction with radiomics data.
}
\author{
Carel F.W. Peeters <cf.peeters@vumc.nl>
}
\note{
Different from usual scree (Cattell, 1966).
Indeed, usual scree operates on raw correlation matrix.
Here, under a ML factor model.

\itemize{
\item The argument \code{m} cannot exceed the Ledermann-bound (Ledermann, 1937): \eqn{\lfloor [2p + 1 - (8p + 1)^{1/2}]/2\rfloor}, where \eqn{p} indicates the observed-feature dimension.
}
}
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## Simulate some high-dimensional data according to the factor model
simDAT <- FAsim(p = 50, m = 5, n = 40)

## Regularize the correlation matrix
RegR <- regcor(simDAT$data)

## Produce scree plots
MLeig  <- scree(RegR$optCor, m = 5, type = "ML"); MLeig
REDeig <- scree(RegR$optCor, m = 5, type = "reduced"); REDeig
}
\concept{Scree plot}



  
  # Return object
  EVt <- as.data.frame(cbind(dims, EVs))
  colnames(EVt) <- c("latent dimension", "Eigenvalue")
  
  # Return
  return(EVt)
}