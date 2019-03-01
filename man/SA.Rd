\name{SA}
\alias{SA}
\title{
Calculate the KMO measure of sampling adequacy
}
\description{
Function that calculates the Kaiser-Meyer-Olkin (KMO) measure of sampling adequacy.
}
\usage{
SA(R)
}
\arguments{
  \item{R}{
  (Regularized) covariance or correlation \code{matrix}.
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Carel F.W. Peeters <cf.peeters@vumc.nl>
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##
}
\concept{KMO-index}


SA <- function(R){
  ##############################################################################
  # Calculate Kaiser-Meyer-Olkin measure of sampling adequacy
  # R > (regularized) covariance or correlation matrix
  #
  # NOTES:
  # - On interpretation
  # - Kaiser KMO above .9
  # - About factorability
  ##############################################################################  
  
  # Dependencies:
  # require("base")
  
  # Checks
  if (!is.matrix(R)){
    stop("Input (R) should be a matrix")
  }
  if (nrow(R) != ncol(R)){
    stop("Input (R) should be square matrix")
  }
  
  # Preliminaries
  Ri      <- solve(R)
  Sc      <- diag(sqrt((1/diag(Ri))))
  P       <- Sc %*% Ri %*% Sc
  diag(R) <- 0
  diag(P) <- 0
  
  # Calculate overall KMO index
  KMO <- sum(R^2)/(sum(R^2) + sum(P^2))
  
  # Calculate KMO index per feature
  KMOfeature <- colSums(R^2)/(colSums(R^2) + colSums(P^2))
  
  # Return
  return(list(KMO = KMO, KMOfeature = KMOfeature))
}