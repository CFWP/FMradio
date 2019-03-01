\name{SA}
\alias{SA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
SA(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
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
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x) 
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line


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