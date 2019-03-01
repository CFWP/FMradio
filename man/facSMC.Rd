\name{facSMC}
\alias{facSMC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
facSMC(x)
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

facSMC <- function(R, LM){
  ##############################################################################
  # Calculate squared multiple correlations for predicting the common factors
  # R  > (regularized) correlation matrix
  # LM > (rotated) loadings matrix
  #
  # NOTES:
  # - Calculates squared multiple correlations between the observed features
  #   and the common factors
  # - The closer to unity, the better one is able to determine the factors scores
  #   from the observed features
  # - Hence, the closer to unity, the lesser the problem of factor-score 
  #   indeterminacy and the better one is able to uniquely determine the 
  #   factor scores
  # - Note that the computations assume an orthogonal factor model
  # - Hence, only orthogonal rotations of the loadings matrix should be used 
  #   (or no rotation at all)
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
    stop("Column-dimension input (R) incompatible with row-dimension input (LM)")    
  }
  
  # Determine squared multiple correlations for predicting factors
  return(diag(t(LM) %*% solve(R) %*% LM))
}