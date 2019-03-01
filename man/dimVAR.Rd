\name{dimVAR}
\alias{dimVAR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
dimVAR(x)
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

dimVAR <- function(R, maxdim, graph = TRUE, verbose = TRUE){
  ##############################################################################
  # Support function assessing proportion and cumulative variances for a 
  # range of factor solutions
  # R        > (regularized) covariance or correlation matrix
  # maxdim   > maximum number of latent factors to be assessed
  # graph    > logical indicating if output should also be visualized
  # verbose  > logical indicating if function should run silently
  #
  # NOTES:
  # - maxdim cannot exceed the Ledermann-bound
  # - usually, one wants to set maxdim (much) lower than the Ledermann-bound
  # - factor solutions are assessed ranging from 1 to maxdim
  # - SS, proportion variance and cumulative variance are tabulated for each
  #   factor solution
  # - these are based on the unrotated solution
  # - note that the cumulative variance does not depend on the choice of 
  #   (orthogonal) rotation
  # - when graph = TRUE, a graph is returned visualizing the total cumulative
  #   variance against the dimension of the factor solution
  ##############################################################################  
  
  # Dependencies:
  # require("base")
  # require("stats")
  # require("graphics")
  
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
  if (class(maxdim) != "numeric" & class(maxdim) != "integer"){
    stop("Input (maxdim) is of wrong class")
  }
  if (length(maxdim) != 1){
    stop("Length input (maxdim) must be one")
  }
  if (maxdim <= 1){
    stop("Input (maxdim) cannot be lower than 1")
  }
  if (maxdim > mmax){
    stop("Input (maxdim) is too high")
  }
  if (class(graph) != "logical"){ 
    stop("Input (graph) is of wrong class") 
  }
  if (class(verbose) != "logical"){ 
    stop("Input (verbose) is of wrong class") 
  }
  
  # Preliminaries
  dim <- seq(from = 1, to = maxdim, by = 1)
  L   <- list()
  mCV <- numeric()
  
  # Calculate SS, proportion variance, cumulative variance
  for (j in 1:length(dim)){
    m      <- dim[j]
    fit    <- factanal(factors = m, covmat = R, rotation = "none")
    SS     <- diag(t(fit$loadings)%*%(fit$loadings))
    PV     <- SS/p
    CV     <- cumsum(PV)
    mCV[j] <- max(CV)
    L[[j]] <- rbind(SS,PV,CV)
    names(L)[j] <- paste("dimension =", m)
    if (verbose){cat(paste("Calculating statistics for dimension m = ", 
                           m, " done\n", sep = ""))}
  }
  
  # Visualization
  if (graph){
    if (verbose){cat("Visualizing...", "\n")}
    MT <- seq(from = 1, to = maxdim, by = 2)
    plot(dim, mCV, axes = FALSE, type = "l", 
         col = "red", xlab = "dimension of factor solution", 
         ylab = "total cumulative variance")
    axis(2, ylim = c(0, max(mCV)), col = "black", lwd = 1)
    axis(1, xlim = c(1,maxdim), col = "black", lwd = 1, tick = TRUE, at = MT)
  }
  
  # Return
  return(list(CumVar = mCV, varianceTables = L))
}
