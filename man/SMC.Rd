\name{SMC}
\alias{SMC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
SMC(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
}
\details{
Requires run/fit of mlFA.
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

SMC <- function(R, LM){
  ##############################################################################
  # Juxtaposing SMCs and estimates of communalities
  # R  > (regularized) covariance or correlation matrix
  # LM > (rotated) loadings matrix
  #
  # NOTES:
  # - Squared multiple correlations (SMCs) refer to the proportion of variance 
  #   in feature j that is explained by the remaining p - 1 features
  # - Note that the choice of orthogonal rotation does not affect the 
  #   communality estimates
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
  if (!is.matrix(LM)){
    stop("Input (LM) should be a matrix")    
  } 
  if (ncol(R) != nrow(LM)){
    stop("Column (or row) dimension input (R) incompatible with row-dimension input (LM)")    
  }
  
  # Calculate SMCs and communalities
  SMC           <- 1 - 1/diag(solve(R))
  Communalities <- diag((LM) %*% t(LM))
  
  # Return
  return(cbind(SMC, Communalities))
}
