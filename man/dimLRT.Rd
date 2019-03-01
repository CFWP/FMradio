\name{dimLRT}
\alias{dimLRT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
dimLRT(x)
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

dimLRT <- function(R, X, maxdim, rankDOF = TRUE, graph = TRUE, 
                   alpha = .05, Bartlett = FALSE, verbose = TRUE){
  ##############################################################################
  # Performs dimensionality assessment using likelihood ratio testing
  # R        > (regularized) covariance or correlation matrix
  # X        > centered and scaled (subsetted) data matrix, observations in rows
  # maxdim   > maximum number of latent factors to be assessed
  # rankDOF  > logical indicating if the degrees of freedom should be based on
  #            the rank of the raw correlation matrix
  # graph    > logical indicating if output should also be visualized
  # alpha    > numeric giving alpha level, only used when graph = TRUE
  # Bartlett > logical indicating if Bartlett correction should be applied
  # verbose  > logical indicating if function should run silently
  #
  # NOTES:
  # - maxdim cannot exceed the Ledermann-bound
  # - usually, one wants to set maxdim (much) lower than the Ledermann-bound
  # - note that, if p > n, the the maximum rank of the raw correlation matrix
  #   is n - 1. In this case there is an alternative Ledermann-bound when 
  #   rankDOF = TRUE. The number of information points in the correlation matrix
  #   is then given as (n-1)*n/2 and this number must exceed 
  #   p*maxdim + p - (maxdim*(maxdim - 1))/2, putting more restrictions on 
  #   maxdim
  ##############################################################################  
  
  # Dependencies:
  # require("base")
  # require("stats")
  # require("graphics")
  
  # Checks
  if (!is.matrix(R)){
    stop("Input (R) should be a matrix")
  }
  if (nrow(R) != ncol(R)){
    stop("Input (R) should be square matrix")
  }
  if (!is.matrix(X)){
    stop("Input (X) should be a matrix")
  }
  
  # Preliminaries for further checks
  p    <- ncol(R)
  n    <- dim(X)[1]
  pr   <- .rank(cor(X))
  mmax <- floor((2*p + 1 - sqrt(8*p + 1))/2)
  
  # (Conditional) Checks
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
  if (class(rankDOF) != "logical"){ 
    stop("Input (rankDOF) is of wrong class") 
  }
  if (class(graph) != "logical"){ 
    stop("Input (graph) is of wrong class") 
  }
  if (graph){
    if (class(alpha) != "numeric"){
      stop("Input (alpha) is of wrong class") 
    }
    if (length(alpha) != 1){
      stop("Length input (alpha) must be one")
    }
    if (alpha < 0){
      stop("Input (alpha) should be strictly positive") 
    }
    if (alpha > 1){
      stop("Input (alpha) cannot exceed unity") 
    }
  }
  if (class(Bartlett) != "logical"){ 
    stop("Input (Bartlett) is of wrong class") 
  }
  if (class(verbose) != "logical"){ 
    stop("Input (verbose) is of wrong class") 
  }
  if (p >= n & rankDOF == TRUE){
    if ((pr*(pr + 1))/2 - (p*maxdim + p - (maxdim*(maxdim - 1))/2) < 1){
      stop("Input (maxdim) implies negative degrees of freedom for LRT")
    }
  }
  
  # Preliminaries
  LRT  <- numeric()
  pval <- numeric()
  dim  <- seq(from = 1, to = maxdim, by = 1)
  
  # Calculate LRTs and p-values
  for (j in 1:length(dim)){
    m          <- dim[j]
    fit        <- factanal(factors = m, covmat = R, rotation = "none")
    loadings   <- fit$loadings[1:p,]
    Uniqueness <- diag(fit$uniquenesses)
    if (Bartlett){
      LRT[j] <- (n - 1 - (2 * p + 5)/6 - (2 * m)/3) * .DF(R, loadings, Uniqueness)
    } else {
      LRT[j] <- (n - 1) * .DF(R, loadings, Uniqueness)
    }
    if (rankDOF){
      dof <- (pr*(pr + 1))/2 - (p*m + p - (m*(m - 1))/2)
    } else {
      dof <- ((p - m)^2 - (p + m))/2
    }
    pval[j] <- pchisq(LRT[j], df = dof, lower.tail = FALSE)
    if (verbose){cat(paste("Calculating likelihood ratio test for dimension m = ", 
                           m, " done\n", sep = ""))}
  }
  
  # Visualization
  if (graph){
    if (verbose){cat("Visualizing...", "\n")}
    MT <- seq(from = 1, to = maxdim, by = 2)
    plot(dim, pval, axes = FALSE, type = "l", 
         col = "red", xlab = "dimension of factor solution", 
         ylab = "p-value LR test")
    axis(2, ylim = c(0, 1), col = "black", lwd = 1)
    axis(1, xlim = c(1,maxdim), col = "black", lwd = 1, tick = TRUE, at = MT)
    abline(h = alpha, col = "blue")
  }
  
  # Return object
  LRTt <- as.data.frame(cbind(dim, LRT, pval))
  colnames(LRTt) <- c("latent dimension", "Statistic", "p.value")
  
  ## Return
  return(LRTt)
}