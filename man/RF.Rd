\name{RF}
\alias{RF}
\title{
Redundancy filtering of  a square (correlation) matrix
}
\description{
Function that performs redundancy filtering (RF) of a square (correlation) matrix.
}
\usage{
RF(R, t = .9)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
}
\details{
Correlation matrices can have collinearity.
This collinearity can be so strong as to border on redundancy.
Hence, redundant features can be removed.

Uses algorithm developed in Peeters et al.
See SM of Peeters.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~ Own paper should go here.
}
\author{
Carel F.W. Peeters <cf.peeters@vumc.nl>
}
\note{
\itemize{
\item While geared towards the reundancy filtering of correlation matrices, the function is quite general, in the sense that it can be use to filter any square \code{matrix}.
\item When the input matrix \code{R} is a correlation matrix, then argument \code{t} should satisfy \eqn{-1 < \mathrm{t} < 1}, for the return matrix to be sensical for further analysis.
}
}
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~ SubSet convenience function
}
\examples{
EXAMPLE where heatmap is used to check the redundancy filtering
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line





RF <- function(R, t = .9){
  ##############################################################################
  # Performs redundancy filtering (RF) of square (correlation) matrix
  # R > (correlation) matrix
  # t > absolute value for thresholding
  #
  # NOTES:
  # - When the input matrix R is a correlation matrix, then t should satisfy
  #   -1 < t < 1, for the return matrix to be sensical for further analysis
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
  if (class(t) != "numeric"){
    stop("Input (t) is of wrong class")
  }
  if (length(t) != 1){
    stop("Length input (t) must be one")
  }
  
  # Needed
  GO = TRUE
  HM <- numeric()
  
  # Loop
  while(GO){
    for(j in 1:dim(R)[1]){
      HM[j] <- sum(abs(R[j,, drop = FALSE]) >= t)
    }
    Mx <- which(HM == max(HM))
    if(max(HM) < 2){
      GO = FALSE
    } else {
      R  <- R[-Mx[1],-Mx[1]]
      HM <- numeric()
    }
  }
  
  # Return
  return(R)
}

