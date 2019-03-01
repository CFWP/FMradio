\name{mlFA}
\alias{mlFA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
mlFA(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
}
\details{
Wrapper around factanal.
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

mlFA <- function(R, m){
  ##############################################################################
  # Fits an m-factor model through maximum likelihood estimation
  # Produces a Varimax-rotated solution
  # R > (regularized) covariance or correlation matrix
  # m > dimension of factor solution
  # 
  # NOTES:
  # - this function is basically a wrapper around the factanal (stats) function
  # - its purpose is to produce a factor solution of the chosen dimension by
  #   maximum likelihood
  # - the wrapper ensures the model is fitted under the same circumstances
  #   under which latent dimensionality is assessed
  # - produces a Varimax-rotated solution
  # - its output can be used to produce factor scores by the facScore function
  # - note that these cannot be obtained using the factanal function as the
  #   fit is based on the (regularized) correlation matrix
  # - the function assumes that the chosen dimension m is reasonable
  # - the order of the ouput is the order of the features as given in the R 
  #   object
  # - as the 'loadings-part' of the output is of class "loadings" (factanal)
  #   it can be subjected to the print function, sorting the output to 
  #   emphasize the simple structure: e.g.,
  #   print(fit$Loadings, digits = 2, cutoff = .3, sort = TRUE)
  ##############################################################################  
  
  # Dependencies:
  # require("base")
  # require("stats")
  
  # Preliminaries for checks
  p    <- ncol(R)
  mmax <- floor((2*p + 1 - sqrt(8*p + 1))/2)
  
  # Checks
  if (!is.matrix(R)){
    stop("Input (R) should be a matrix")
  }
  if (nrow(R) != ncol(R)){
    stop("Input (R) should be square matrix")
  }
  if (class(m) != "numeric" & class(m) != "integer"){
    stop("Input (m) is of wrong class")
  }
  if (length(m) != 1){
    stop("Length input (m) must be one")
  }
  if (m <= 1){
    stop("Input (m) cannot be lower than 1")
  }
  if (m > mmax){
    stop("Input (m) is too high")
  }
  
  # Wrapper
  fit <- factanal(factors = m, covmat = R, rotation = "varimax")
  
  # Return
  rotmatrix <- fit$rotmat
  rownames(rotmatrix) <- colnames(fit$loadings)
  colnames(rotmatrix) <- rownames(rotmatrix)
  return(list(Loadings = fit$loadings, 
              Uniqueness = diag(fit$uniquenesses),
              rotmatrix = rotmatrix))
}