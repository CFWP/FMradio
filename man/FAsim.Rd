\name{FAsim}
\alias{FAsim}
\title{
Simulate data according to the common factor analytic model
}
\description{
\code{FAsim} is a function that enables the simulation of data according to the common factor analytic model.
}
\usage{
FAsim(p, m, n, simplestructure = TRUE, balanced = TRUE,
      loadingfix = TRUE, loadingnegative = TRUE,
      loadingvalue = .8, loadingvaluelow = .2, numloadings,
      loadinglowerH = .7, loadingupperH = .9, 
      loadinglowerL = .1, loadingupperL = .3)
}
\arguments{
  \item{p}{
  A \code{numeric} integer or \code{integer} indicating the number of observed features.
 }
  \item{m}{
  A \code{numeric} integer or \code{integer} indicating the latent dimension of the factor solution (i.e., the number of factors).
 }
  \item{n}{
  A \code{numeric} integer or \code{integer} indicating the number of samples.
 }
  \item{simplestructure}{
  A \code{logical} indicating if the generating factor structure should be factorially pure.
 }
  \item{balanced}{
  A \code{logical} indicating if the high (i.e., qualitatively 'significant') loadings should be divided evenly over the respective factors.
 }
  \item{loadingfix}{
  A \code{logical} indicating if the loadings should have a fixed value.
 }
  \item{loadingnegative}{
  A \code{logical} indicating if, next to positive, also negative loadings should be present.
 }
  \item{loadingvalue}{
  A \code{numeric} indicating the value for high (i.e., qualitatively 'significant') loadings. Used when \code{loadingfix = TRUE}.
 }
  \item{loadingvaluelow}{
  A \code{numeric} indicating the value for low loadings. Used when \code{loadingfix = TRUE} & \code{simplestructure = FALSE}.
 }
  \item{numloadings}{
  A \code{numeric} vector with length equalling argument \code{m}, indicating the number of high (i.e., qualitatively 'significant') loadings per factor.\cr 
  Used when \code{balanced = FALSE}.
 }
  \item{loadinglowerH}{
  A \code{numeric} indicating the lower-bound of high (i.e., qualitatively 'significant') loadings. Used when \code{loadingfix = FALSE}.
 }
  \item{loadingupperH}{
  A \code{numeric} indicating the upper-bound of high (i.e., qualitatively 'significant') loadings. Used when \code{loadingfix = FALSE}.
 }
  \item{loadinglowerL}{
  A \code{numeric} indicating the lower-bound of low (i.e., qualitatively 'non-significant') loadings. Used when \code{loadingfix = FALSE} & \code{simplestructure = FALSE}.
 }
  \item{loadingupperL}{
  A \code{numeric} indicating the upper-bound of low (i.e., qualitatively 'non-significant') loadings. Used when \code{loadingfix = FALSE} & \code{simplestructure = FALSE}.
 }
}
\details{
\code{FAsim} provides certain flexibility when generating data according to an orthogonal common factor-analytic model.
It can produce data according to, for example, (i) factorially pure loadings structures, (ii) loadings-structures with only positive entries or both positive and negative loadings, (iii) loadings-structures with fixed values or varying values, (iv) balanced and unbalanced loadings-structures.
}
\value{
The function returns an object of class \code{list}:
\item{$data}{A standardized data \code{matrix} of size \eqn{n \times p}.}
\item{$loadings}{Loadings \code{matrix} of size \eqn{p \times m} on which the data-generation was based.}
\item{$Uniqueness}{A \code{numeric} vector of size \eqn{p} representing the uniquenesses on which the data-generation was based.}
\item{$cormatrix}{A \eqn{p \times p} correlation \code{matrix} based on the generated data in slot \code{$data}.}
}
\references{
Ledermann, W. (1937). On the rank of the reduced correlational matrix in multiple factor analysis. Psychometrika, 2:85--93.

Peeters, C.F.W. \emph{et al.} (2019). Stable prediction with radiomics data.
\href{https://arxiv.org/abs/1903.11696}{arXiv:1903.11696 [stat.ML]}.
}
\author{
Carel F.W. Peeters <cf.peeters@vumc.nl>
}
\note{
\itemize{
\item A uniform distribution is assumed when generating draws between \code{loadinglowerH} and \code{loadingupperH}.

\item A uniform distribution is assumed when generating draws between \code{loadinglowerL} and \code{loadingupperL}.

\item The argument \code{m} cannot exceed the Ledermann-bound (Ledermann, 1937): \eqn{\lfloor [2p + 1 - (8p + 1)^{1/2}]/2\rfloor}, where \eqn{p} indicates the observed-feature dimension.
}
}
\seealso{
\code{\link{dimGB}}, \code{\link{mlFA}}, \code{\link{facScore}}
}
\examples{
## Simulate some data according to a factor model with 3 latent factors
## Balanced and factorially pure loadings structure
simDAT <- FAsim(p = 24, m = 3, n = 40, loadingvalue = .9)
simDAT$loadings

## Simulate some data according to a factor model with 3 latent factors
## Unbalanced and factorially pure loadings structure
simDAT <- FAsim(p = 24, m = 3, n = 40, loadingvalue = .9,
                balanced = FALSE, numloadings = c(10,10,4))
simDAT$loadings

## Simulate some data according to a factor model with 3 latent factors
## Unbalanced and factorially non-pure loadings structure
simDAT <- FAsim(p = 24, m = 3, n = 40, loadingvalue = .9,
                balanced = FALSE, numloadings = c(10,10,4),
                simplestructure = FALSE)
simDAT$loadings

## Simulate some data according to a factor model with 3 latent factors
## Unbalanced and factorially non-pure loadings structure
## Non-fixed high and low loadings
simDAT <- FAsim(p = 24, m = 3, n = 40, loadingvalue = .9,
                balanced = FALSE, numloadings = c(10,10,4),
                simplestructure = FALSE, loadingfix = FALSE)
simDAT$loadings
}
\concept{Data simulation}
