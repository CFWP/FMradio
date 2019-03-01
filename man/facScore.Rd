\name{facScore}
\alias{facScore}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
facScore(x)
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


facScore <- function(X, LM, UM, type = "thomson"){
  ##############################################################################
  # Finds factor scores from data and given factor-solution
  # Default is Thomson-type factor scores
  # X    > data matrix, observations in rows
  # LM   > (rotated) loadings matrix
  # UM   > diagonal uniquenesses matrix
  # type > character indicating the type of factor score to calculate
  #
  # NOTES:
  # - The input data are assumed to be scaled (or at least centered)
  # - The UM matrix is assumed to be positive definite
  # - The LM matrix is assumed to be of full column rank
  # - The factor-scores obtained the default way are Thomson-type scores
  # - These will be near-orthogonal
  # - These are close to Bartlett-type scores, which are unbiased
  # - Anderson-type scores are orthogonal
  # - The return-object is a dataframe
  ##############################################################################  
  
  # Dependencies:
  # require("base")
  # require("expm")
  
  # Checks
  if (!is.matrix(X)){
    stop("Input (X) should be a matrix")    
  }
  if (!is.matrix(LM)){
    stop("Input (LM) should be a matrix")    
  } 
  if (ncol(X) != nrow(LM)){
    stop("Column-dimension input (X) incompatible with row-dimension input (LM)")    
  }
  if (!is.matrix(UM)){
    stop("Input (UM) should be a matrix")    
  } 
  if (nrow(UM) != ncol(UM)){
    stop("Input (UM) should be square matrix")
  }
  if (ncol(X) != nrow(UM)){
    stop("Column-dimensions input (X) incompatible with dimension input (UM)")    
  }
  if (!(type %in% c("thomson", "bartlett", "anderson"))){
    stop("Input (type) should be one of {'thomson', 'bartlett', 'anderson'}")
  }

  # Needed
  UMi <- diag(1/diag(UM))
  m   <- diag(dim(LM)[2])
  
  # Obtain factor scores
  if (type == "thomson"){
    scores <- as.data.frame(X %*% UMi %*% LM %*% solve(m + t(LM) %*% UMi %*% LM))
  }
  if (type == "bartlett"){
    scores <- as.data.frame(X %*% UMi %*% LM %*% solve(t(LM) %*% UMi %*% LM))
  }
  if (type == "anderson"){
    G <- t(LM) %*% UMi %*% LM
    scores <- as.data.frame(X %*% UMi %*% LM %*% sqrtm(solve((G %*% (m + G)))))
  }
  
  # Return
  return(scores)
}